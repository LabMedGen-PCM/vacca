#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v conda-build >/dev/null 2>&1; then
  echo "conda-build is required. Install it with:" >&2
  echo "  mamba install -n base -c conda-forge conda-build conda-index" >&2
  exit 2
fi
conda build "$ROOT_DIR/conda-recipe"
