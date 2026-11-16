#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${CNV_VALIDATOR_ENV:-cnv-vaf-validator}"
HEARTBEAT_SECONDS="${CNV_INSTALL_HEARTBEAT:-15}"
INSTALL_LOG="${CNV_INSTALL_LOG:-$ROOT_DIR/install.log}"
PACKED_ENV_ARCHIVE="${CNV_PACKED_ENV_ARCHIVE:-$ROOT_DIR/dependencies/cnv-vaf-validator.tar.gz}"
ENV_INSTALL_RETRIES="${CNV_INSTALL_RETRIES_BEFORE_FALLBACK:-0}"
INSTALL_MODE="local"
PACKED_FALLBACK_USED=false
PACKED_PREFIX=""

usage() {
  cat <<'EOF'
Usage: bash install.sh [--local | --online]

Installation modes:
  --local   Install directly from dependencies/cnv-vaf-validator.tar.gz.
            This is the default mode and does not download packages.
  --online  First try to create/update the environment from environment.yml.
            If that fails, fall back to the bundled local archive.

Other options:
  -h, --help  Show this help message and exit.

Environment variables:
  CNV_VALIDATOR_ENV                    Target environment name.
  CNV_PACKED_ENV_ARCHIVE               Alternative conda-pack archive path.
  CNV_INSTALL_RETRIES_BEFORE_FALLBACK  Online retries before local fallback.
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --local)
      INSTALL_MODE="local"
      ;;
    --online)
      INSTALL_MODE="online"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Error: unknown argument: %s\n\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

# Keep package downloads reliable on hosts where several simultaneous
# conda-forge transfers can stall. Explicit user settings still take priority.
export MAMBA_DOWNLOAD_THREADS="${MAMBA_DOWNLOAD_THREADS:-1}"
export MAMBA_REMOTE_MAX_RETRIES="${MAMBA_REMOTE_MAX_RETRIES:-10}"
export MAMBA_REMOTE_BACKOFF_FACTOR="${MAMBA_REMOTE_BACKOFF_FACTOR:-2}"

TOTAL_STEPS=6
CURRENT_STEP=0
CURRENT_STEP_NAME="initialization"

if [[ ! "$HEARTBEAT_SECONDS" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Error: CNV_INSTALL_HEARTBEAT must be a positive integer.\n' >&2
  exit 2
fi
if [[ ! "$ENV_INSTALL_RETRIES" =~ ^[0-9]+$ ]]; then
  printf 'Error: CNV_INSTALL_RETRIES_BEFORE_FALLBACK must be a non-negative integer.\n' >&2
  exit 2
fi

# Start a fresh log on every run. Set CNV_INSTALL_LOG to save it elsewhere.
: > "$INSTALL_LOG"

color_enabled=false
if [[ -t 1 ]] && command -v tput >/dev/null 2>&1 && [[ "$(tput colors 2>/dev/null || echo 0)" -ge 8 ]]; then
  color_enabled=true
fi

if $color_enabled; then
  BOLD="$(tput bold)"
  RESET="$(tput sgr0)"
  GREEN="$(tput setaf 2)"
  YELLOW="$(tput setaf 3)"
  RED="$(tput setaf 1)"
  CYAN="$(tput setaf 6)"
else
  BOLD=""
  RESET=""
  GREEN=""
  YELLOW=""
  RED=""
  CYAN=""
fi

log() {
  printf '%s\n' "$*" | tee -a "$INSTALL_LOG"
}

render_progress() {
  local state="${1:-RUNNING}"
  local width=30
  local percent=$(( CURRENT_STEP * 100 / TOTAL_STEPS ))
  local filled=$(( CURRENT_STEP * width / TOTAL_STEPS ))
  local empty=$(( width - filled ))
  local bar_fill bar_empty
  printf -v bar_fill '%*s' "$filled" ''
  printf -v bar_empty '%*s' "$empty" ''
  bar_fill=${bar_fill// /#}
  bar_empty=${bar_empty// /-}
  printf '%s[%s%s] %3d%%  Step %d/%d: %s (%s)%s\n' \
    "$CYAN" "$bar_fill" "$bar_empty" "$percent" \
    "$CURRENT_STEP" "$TOTAL_STEPS" "$CURRENT_STEP_NAME" "$state" "$RESET" \
    | tee -a "$INSTALL_LOG"
}

begin_step() {
  CURRENT_STEP=$((CURRENT_STEP + 1))
  CURRENT_STEP_NAME="$1"
  log ""
  render_progress "STARTED"
}

finish_step() {
  render_progress "DONE"
}

print_command() {
  printf '%sCommand:%s ' "$BOLD" "$RESET" | tee -a "$INSTALL_LOG"
  printf '%q ' "$@" | tee -a "$INSTALL_LOG"
  printf '\n' | tee -a "$INSTALL_LOG"
}

print_process_tree() {
  local pid="$1"
  log "${YELLOW}Still running after ${2}s. Active process tree:${RESET}"
  if command -v pstree >/dev/null 2>&1; then
    pstree -ap "$pid" 2>/dev/null | tee -a "$INSTALL_LOG" || true
  else
    ps -o pid,ppid,etime,%cpu,%mem,stat,cmd -p "$pid" --ppid "$pid" 2>/dev/null \
      | tee -a "$INSTALL_LOG" || true
  fi
}

run_visible() {
  print_command "$@"

  # Run through tee so command output remains visible and is also recorded.
  (
    set -o pipefail
    "$@" 2>&1 | tee -a "$INSTALL_LOG"
  ) &
  local wrapper_pid=$!
  local elapsed=0

  log "Started PID $wrapper_pid. Heartbeat interval: ${HEARTBEAT_SECONDS}s"

  while kill -0 "$wrapper_pid" 2>/dev/null; do
    sleep 1
    elapsed=$((elapsed + 1))
    if (( elapsed % HEARTBEAT_SECONDS == 0 )) && kill -0 "$wrapper_pid" 2>/dev/null; then
      print_process_tree "$wrapper_pid" "$elapsed"
    fi
  done

  local status=0
  wait "$wrapper_pid" || status=$?
  if (( status != 0 )); then
    log "${RED}Command failed with exit status $status.${RESET}"
    log "Installation log: $INSTALL_LOG"
    return "$status"
  fi
}

on_error() {
  local status=$?
  log ""
  log "${RED}${BOLD}Installation failed during step $CURRENT_STEP/$TOTAL_STEPS: $CURRENT_STEP_NAME${RESET}"
  log "Exit status: $status"
  log "Log file: $INSTALL_LOG"
  log "Last 40 log lines:"
  tail -n 40 "$INSTALL_LOG" 2>/dev/null || true
  exit "$status"
}

on_interrupt() {
  log ""
  log "${RED}Installation interrupted during: $CURRENT_STEP_NAME${RESET}"
  log "Log file: $INSTALL_LOG"
  exit 130
}

trap on_error ERR
trap on_interrupt INT TERM

find_manager() {
  local candidate
  for candidate in mamba micromamba conda; do
    if command -v "$candidate" >/dev/null 2>&1; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

resolve_packed_environment_prefix() {
  local root_prefix="" candidate="" manager_path=""

  # Prefer explicit prefixes. Validate every candidate before accepting it:
  # some Mamba versions print an error message to stdout for unsupported
  # info options, which must not be mistaken for a filesystem path.
  if [[ -n "${MAMBA_ROOT_PREFIX:-}" && -d "$MAMBA_ROOT_PREFIX" ]]; then
    root_prefix="$MAMBA_ROOT_PREFIX"
  fi
  if [[ -z "$root_prefix" && -n "${CONDA_EXE:-}" ]]; then
    candidate="$(cd "$(dirname "$CONDA_EXE")/.." 2>/dev/null && pwd || true)"
    if [[ -n "$candidate" && -d "$candidate" ]]; then
      root_prefix="$candidate"
    fi
  fi
  if [[ -z "$root_prefix" ]]; then
    manager_path="$(command -v "$MANAGER")"
    candidate="$(cd "$(dirname "$manager_path")/.." 2>/dev/null && pwd || true)"
    if [[ -n "$candidate" && -d "$candidate" ]]; then
      root_prefix="$candidate"
    fi
  fi
  if [[ -z "$root_prefix" ]]; then
    candidate="$("$MANAGER" info --base 2>/dev/null | tail -n 1 || true)"
    if [[ -n "$candidate" && -d "$candidate" ]]; then
      root_prefix="$candidate"
    fi
  fi
  if [[ -z "$root_prefix" ]]; then
    candidate="$("$MANAGER" info --root-prefix 2>/dev/null | tail -n 1 || true)"
    if [[ -n "$candidate" && -d "$candidate" ]]; then
      root_prefix="$candidate"
    fi
  fi

  if [[ -z "$root_prefix" || ! -d "$root_prefix" ]]; then
    log "${RED}Error: failed to resolve the Conda/Mamba root prefix for the packed-environment fallback.${RESET}"
    return 4
  fi

  printf '%s/envs/%s\n' "${root_prefix%/}" "$ENV_NAME"
}

restore_after_packed_failure() {
  local target="$1"
  local backup="$2"
  local failed="${target}.packed-failed-$(date +%Y%m%d-%H%M%S)"

  if [[ -e "$target" ]]; then
    mv -- "$target" "$failed"
    log "Failed packed environment preserved at: $failed"
  fi
  if [[ -n "$backup" && -e "$backup" ]]; then
    mv -- "$backup" "$target"
    log "Previous environment restored at: $target"
  fi
}

install_packed_environment() {
  local target parent temporary backup=""
  local timestamp

  if [[ ! -r "$PACKED_ENV_ARCHIVE" ]]; then
    log "${RED}Error: packed fallback archive is missing or unreadable: $PACKED_ENV_ARCHIVE${RESET}"
    return 4
  fi
  if ! command -v tar >/dev/null 2>&1; then
    log "${RED}Error: tar is required to install the packed fallback environment.${RESET}"
    return 4
  fi

  target="$(resolve_packed_environment_prefix)"
  parent="$(dirname "$target")"
  temporary="${target}.packed-tmp.$$"
  timestamp="$(date +%Y%m%d-%H%M%S)"

  mkdir -p "$parent"
  if [[ -e "$temporary" ]]; then
    log "${RED}Error: temporary fallback path already exists: $temporary${RESET}"
    return 4
  fi
  mkdir -p "$temporary"

  if [[ "$INSTALL_MODE" == "online" ]]; then
    log "Normal Conda/Mamba installation failed; using the bundled conda-pack environment."
  else
    log "Installing directly from the bundled conda-pack environment."
  fi
  log "Packed environment archive: $PACKED_ENV_ARCHIVE"
  log "Packed environment target: $target"
  if ! run_visible tar -xzf "$PACKED_ENV_ARCHIVE" -C "$temporary"; then
    log "${RED}Error: failed to extract the packed environment.${RESET}"
    log "Incomplete extraction preserved at: $temporary"
    return 4
  fi

  if [[ ! -x "$temporary/bin/python" || ! -d "$temporary/conda-meta" ]]; then
    log "${RED}Error: the packed archive does not look like a complete Conda environment.${RESET}"
    log "Extracted files preserved at: $temporary"
    return 4
  fi

  if [[ -e "$target" ]]; then
    backup="${target}.pre-packed-backup-${timestamp}"
    mv -- "$target" "$backup"
    log "Existing or partially created environment backed up to: $backup"
  fi
  mv -- "$temporary" "$target"

  if [[ ! -x "$target/bin/conda-unpack" ]]; then
    log "${RED}Error: conda-unpack is absent from the packed environment.${RESET}"
    restore_after_packed_failure "$target" "$backup"
    return 4
  fi
  if ! run_visible "$target/bin/conda-unpack"; then
    log "${RED}Error: conda-unpack failed for the packed environment.${RESET}"
    restore_after_packed_failure "$target" "$backup"
    return 4
  fi

  PACKED_FALLBACK_USED=true
  PACKED_PREFIX="$target"
  log "Bundled conda-pack environment installed successfully."
  if [[ -n "$backup" ]]; then
    log "Previous environment backup retained at: $backup"
  fi
}

begin_step "checking prerequisites"
MANAGER="$(find_manager || true)"
if [[ -z "$MANAGER" ]]; then
  log "${RED}Error: mamba, micromamba, or conda is required.${RESET}"
  exit 2
fi
log "Environment manager: $MANAGER ($(command -v "$MANAGER"))"
if [[ "$MANAGER" == "mamba" || "$MANAGER" == "micromamba" ]]; then
  log "Mamba download settings: threads=$MAMBA_DOWNLOAD_THREADS, retries=$MAMBA_REMOTE_MAX_RETRIES, backoff=$MAMBA_REMOTE_BACKOFF_FACTOR"
fi
log "Target environment: $ENV_NAME"
log "Installation mode: $INSTALL_MODE"
log "Normal installation retries before packed fallback: $ENV_INSTALL_RETRIES"
log "Installer directory: $ROOT_DIR"
log "Log file: $INSTALL_LOG"
finish_step

begin_step "checking the target Conda environment"
ENV_EXISTS=false
if "$MANAGER" env list 2>/dev/null | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  ENV_EXISTS=true
  if [[ "$INSTALL_MODE" == "local" ]]; then
    log "Environment '$ENV_NAME' already exists and will be replaced from the bundled archive."
  else
    log "Environment '$ENV_NAME' already exists and will be updated."
  fi
else
  if [[ "$INSTALL_MODE" == "local" ]]; then
    log "Environment '$ENV_NAME' does not exist and will be installed from the bundled archive."
  else
    log "Environment '$ENV_NAME' does not exist and will be created."
  fi
fi
finish_step

begin_step "creating or updating the Conda environment"
if [[ "$INSTALL_MODE" == "local" ]]; then
  install_packed_environment
else
  if $ENV_EXISTS; then
    if ! run_visible env \
      MAMBA_REMOTE_MAX_RETRIES="$ENV_INSTALL_RETRIES" \
      CONDA_REMOTE_MAX_RETRIES="$ENV_INSTALL_RETRIES" \
      "$MANAGER" env update -y -n "$ENV_NAME" -f "$ROOT_DIR/environment.yml" --prune; then
      install_packed_environment
    fi
  else
    if ! run_visible env \
      MAMBA_REMOTE_MAX_RETRIES="$ENV_INSTALL_RETRIES" \
      CONDA_REMOTE_MAX_RETRIES="$ENV_INSTALL_RETRIES" \
      "$MANAGER" env create -y -n "$ENV_NAME" -f "$ROOT_DIR/environment.yml"; then
      install_packed_environment
    fi
  fi
fi
finish_step

begin_step "resolving the environment installation prefix"
if $PACKED_FALLBACK_USED; then
  PREFIX="$PACKED_PREFIX"
  log "Using packed environment prefix: $PREFIX"
else
  print_command "$MANAGER" run -n "$ENV_NAME" python -c 'import sys; print(sys.prefix)'
  PREFIX="$($MANAGER run -n "$ENV_NAME" python -c 'import sys; print(sys.prefix)')"
fi
if [[ -z "$PREFIX" || ! -d "$PREFIX" ]]; then
  log "${RED}Error: failed to resolve a valid environment prefix.${RESET}"
  exit 3
fi
log "Environment prefix: $PREFIX"
mkdir -p "$PREFIX/bin" "$PREFIX/share/doc/cnv-vaf-validator" \
  "$PREFIX/share/cnv-vaf-validator/models" "$PREFIX/share/cnv-vaf-validator/data" \
  "$PREFIX/share/cnv-vaf-validator/example"
finish_step

begin_step "installing the validator and documentation"
install -m 0755 "$ROOT_DIR/cnv_vaf_validator.py" "$PREFIX/bin/cnv-vaf-validator"
install -m 0755 "$ROOT_DIR/cnv_vaf_validator.py" "$PREFIX/bin/cnv-vaf-validator.py"
install -m 0644 "$ROOT_DIR/readme.md" "$PREFIX/share/doc/cnv-vaf-validator/readme.md"
install -m 0644 "$ROOT_DIR/help.md" "$PREFIX/share/doc/cnv-vaf-validator/help.md"
install -m 0644 "$ROOT_DIR/hg38_centromeres_telomeres.bed" \
  "$PREFIX/share/cnv-vaf-validator/hg38_centromeres_telomeres.bed"
install -m 0644 "$ROOT_DIR/models/CNN_coverage.pt" \
  "$PREFIX/share/cnv-vaf-validator/models/CNN_coverage.pt"
install -m 0644 "$ROOT_DIR/models/svc_vaf.pkl" \
  "$PREFIX/share/cnv-vaf-validator/models/svc_vaf.pkl"
install -m 0644 "$ROOT_DIR/data/VAF_train.tsv" \
  "$PREFIX/share/cnv-vaf-validator/data/VAF_train.tsv"
install -m 0644 "$ROOT_DIR/example/"* \
  "$PREFIX/share/cnv-vaf-validator/example/"
log "Installed command: $PREFIX/bin/cnv-vaf-validator"
log "Installed alias: $PREFIX/bin/cnv-vaf-validator.py"
log "Installed documentation: $PREFIX/share/doc/cnv-vaf-validator"
log "Installed annotation BED: $PREFIX/share/cnv-vaf-validator/hg38_centromeres_telomeres.bed"
log "Installed CNN checkpoint: $PREFIX/share/cnv-vaf-validator/models/CNN_coverage.pt"
log "Installed VAF SVC model: $PREFIX/share/cnv-vaf-validator/models/svc_vaf.pkl"
log "Installed VAF training table: $PREFIX/share/cnv-vaf-validator/data/VAF_train.tsv"
log "Installed examples: $PREFIX/share/cnv-vaf-validator/example"
finish_step

begin_step "running installation self-tests"
run_visible "$MANAGER" run -n "$ENV_NAME" cnv-vaf-validator --version
run_visible "$MANAGER" run -n "$ENV_NAME" cnv-vaf-validator.py --version
run_visible "$MANAGER" run -n "$ENV_NAME" cnv-vaf-validator --help
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'import joblib, matplotlib, numpy, pandas, PIL, sklearn, torch; assert hasattr(numpy, "trapezoid"); print("Python dependencies: OK")'
run_visible "$MANAGER" run -n "$ENV_NAME" bcftools --version
run_visible "$MANAGER" run -n "$ENV_NAME" samtools --version
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'from pathlib import Path; import sys; path = Path(sys.argv[1]); rows = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]; assert rows and all(len(line.split("\t")) >= 4 for line in rows); print(f"Annotation BED: OK ({len(rows)} intervals)")' "$PREFIX/share/cnv-vaf-validator/hg38_centromeres_telomeres.bed"
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'from pathlib import Path; import sys; path = Path(sys.argv[1]); assert path.stat().st_size > 0; print("CNN checkpoint: OK")' "$PREFIX/share/cnv-vaf-validator/models/CNN_coverage.pt"
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'import joblib, sys; from sklearn.svm import SVC; model = joblib.load(sys.argv[1]); params = model.get_params(deep=False); assert isinstance(model, SVC); assert params["C"] == 1.0 and params["class_weight"] is None and params["coef0"] == 0.0 and params["gamma"] == "scale" and params["kernel"] == "sigmoid"; assert model.probability is True; assert model.feature_names_in_.tolist() == ["X0.5", "X1", "X0.33_X0.66"]; assert model.classes_.tolist() == [0, 1, 2]; print("VAF SVC model: OK")' "$PREFIX/share/cnv-vaf-validator/models/svc_vaf.pkl"
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'import csv, sys; rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t")); assert rows; assert set(rows[0]) == {"SVTYPE", "X0.33_X0.66", "X0.5", "X1"}; assert {row["SVTYPE"] for row in rows} <= {"DEL", "DIP", "DUP"}; print(f"VAF training table: OK ({len(rows)} rows)")' "$PREFIX/share/cnv-vaf-validator/data/VAF_train.tsv"
run_visible "$MANAGER" run -n "$ENV_NAME" python -c 'from pathlib import Path; import sys; root = Path(sys.argv[1]); required = ["ClinCNV_output_chr22.tsv", "HG00577_example.cram", "HG00577_example.vcf.gz"]; assert all((root / name).is_file() for name in required); print("Bundled examples: OK")' "$PREFIX/share/cnv-vaf-validator/example"
finish_step

log ""
log "${GREEN}${BOLD}Installation completed.${RESET}"
log ""
log "Activate the environment:"
log "  $MANAGER activate $ENV_NAME"
log ""
log "Commands:"
log "  cnv-vaf-validator --help"
log "  cnv-vaf-validator.py --help"
log ""
log "Documentation:"
log "  $PREFIX/share/doc/cnv-vaf-validator/readme.md"
log "  $PREFIX/share/doc/cnv-vaf-validator/help.md"
log ""
log "Full installation log:"
log "  $INSTALL_LOG"
