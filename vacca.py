#!/usr/bin/env python3
"""Validate CNV calls from VAF distributions and plot CNV coverage profiles.

The program:
  1. auto-detects a standard CNV table, ClinCNV output, or Manta VCF;
  2. optionally validates CNVs from VAF distributions with the bundled SVC;
  3. optionally validates the same CNVs with a coverage-image CNN;
  4. creates full coverage profiles by default;
  5. writes separate per-sample VAF, CNN, and combined TSV files.

Modes:
  --mode validate
      Run VAF validation, CNN validation, or both, and optionally create full
      coverage profiles.

  --mode plot-cov
      Generate coverage profiles from BAM/CRAM files with samtools depth.

VCF naming examples:
  --vcf-dir /path/to/vcf_files
  --vcf-mask '{sample}/{sample}.output.vcf.gz'
      SRR123 -> /path/to/vcf_files/SRR123/SRR123.output.vcf.gz

Alignment naming example for plot-cov mode:
  --alignment-dir /path/to/bam_or_cram
  --alignment-mask '{sample}.markdup.sorted.cram'
      SRR123 -> /path/to/bam_or_cram/SRR123.markdup.sorted.cram

Absolute templates are also accepted.

Exit codes:
  0: completed (individual rows may still have a non-ok status)
  2: invalid command-line arguments or input table
  3: no rows could be predicted
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    import numpy as np
    import pandas as pd
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing Python dependencies. Install them with:\n"
        "  python -m pip install numpy pandas\n"
        f"Original error: {exc}"
    ) from exc


VERSION = "1.0.0"
FEATURE_COLUMNS = ["X0.33", "X0.66", "X0.5", "X1", "X0.33_X0.66"]
MODEL_FEATURE_COLUMNS = ["X0.5", "X1", "X0.33_X0.66"]
VAF_MODEL_FILENAME = "svc_vaf.pkl"
VAF_TRAIN_FILENAME = "VAF_train.tsv"
NN_MODEL_FILENAME = "CNN_coverage.pt"
VAF_MODEL_CLASSES = ("DEL", "DIP", "DUP")
VAF_MODEL_PARAMETERS = {
    "C": 1.0,
    "class_weight": None,
    "coef0": 0.0,
    "gamma": "scale",
    "kernel": "sigmoid",
}
NN_CLASS_NAMES = ("DEL", "DIP", "DUP")
VALIDATION_SVTYPES = frozenset({"DEL", "DIP", "DUP"})
PLOT_COV_SVTYPES = VALIDATION_SVTYPES
CELL_BLOCK_LIMIT_MESSAGE = "Exceeded cell block limit"
FIRST_COVERAGE_WINDOW_BP = 100
VERBOSE_ENABLED = True

DEFAULT_INTERVALS = {
    "X0.33": (0.28, 0.38),
    "X0.66": (0.61, 0.71),
    "X0.5": (0.45, 0.55),
    "X1": (0.95, 1.00),
}


class ToolError(RuntimeError):
    """Expected user-facing error."""


@dataclass(frozen=True)
class TableInput:
    frame: pd.DataFrame
    single_sample: bool
    base_columns: tuple[str, ...]


@dataclass(frozen=True)
class Interval:
    row_index: int
    start: int
    stop: int


@dataclass(frozen=True)
class NNPlotTask:
    """All picklable inputs needed to create one CNN coverage image."""

    row_index: int
    run: str
    cnv_id: str
    alignment_path: Path
    chrom: str
    view_start: int
    view_stop: int
    cnv_start: int
    cnv_stop: int
    output_path: Path
    coverage_cap_multiplier: float
    dpi: int
    samtools_threads: int


@dataclass(frozen=True)
class NNPlotResult:
    """Result returned by an NN image worker process."""

    row_index: int
    output_path: Path
    error: str | None = None


class SweepAssigner:
    """Assign sorted variant positions to every overlapping interval."""

    def __init__(self, intervals: Sequence[Interval], values: dict[int, list[float]]):
        self.intervals = sorted(intervals, key=lambda item: (item.start, item.stop, item.row_index))
        self.values = values
        self.next_interval = 0
        self.active: set[int] = set()
        self.end_heap: list[tuple[int, int]] = []
        self.last_position = -1

    def add(self, position: int, vaf: float) -> None:
        if position < self.last_position:
            raise ToolError("VCF/query records are not sorted by position")
        self.last_position = position

        while (
            self.next_interval < len(self.intervals)
            and self.intervals[self.next_interval].start <= position
        ):
            interval = self.intervals[self.next_interval]
            self.active.add(interval.row_index)
            heapq.heappush(self.end_heap, (interval.stop, interval.row_index))
            self.next_interval += 1

        while self.end_heap and self.end_heap[0][0] < position:
            _, row_index = heapq.heappop(self.end_heap)
            self.active.discard(row_index)

        for row_index in self.active:
            self.values[row_index].append(vaf)


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr, flush=True)


def vprint(*args: object) -> None:
    """Print an immediately visible progress message when verbose mode is enabled."""
    if VERBOSE_ENABLED:
        eprint("[INFO]", *args)


def canonical_chrom(value: object) -> str:
    chrom = str(value).strip()
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    chrom = chrom.upper()
    if chrom == "M":
        chrom = "MT"
    return chrom


def make_masked_path(base_dir: Path, mask: str, sample: str) -> Path:
    """Resolve a sample-aware relative or absolute path template."""
    try:
        if "{sample}" in mask or "{Run}" in mask:
            filename = mask.format(sample=sample, Run=sample)
        elif "{}" in mask:
            filename = mask.format(sample)
        else:
            filename = f"{sample}{mask}"
    except (KeyError, IndexError, ValueError) as exc:
        raise ToolError(f"Invalid path mask {mask!r}: {exc}") from exc
    return base_dir / filename


def make_vcf_path(vcf_dir: Path, mask: str, sample: str) -> Path:
    return make_masked_path(vcf_dir, mask, sample)


def _column_is_numeric(frame: pd.DataFrame, position: int) -> bool:
    if position >= frame.shape[1]:
        return False
    values = pd.to_numeric(frame.iloc[:, position], errors="coerce")
    return bool(values.notna().all())


def _column_looks_like_chrom(frame: pd.DataFrame, position: int) -> bool:
    if position >= frame.shape[1]:
        return False
    values = frame.iloc[:, position].astype(str).str.strip()
    return bool(
        values.str.fullmatch(
            r"(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)",
            case=False,
        ).all()
    )

def _column_is_svtype(
    frame: pd.DataFrame,
    position: int,
    allowed_svtypes: frozenset[str],
) -> bool:
    if position >= frame.shape[1]:
        return False
    values = frame.iloc[:, position].astype(str).str.upper().str.strip()
    return bool(values.isin(allowed_svtypes).all())


def _output_base_columns(single_sample: bool) -> tuple[str, ...]:
    if single_sample:
        return ("CHR", "START", "END", "INPUT_SVTYPE")
    return ("SAMPLE_ID", "CHR", "START", "END", "INPUT_SVTYPE")


def _finalize_cnv_frame(
    frame: pd.DataFrame,
    *,
    single_sample: bool,
    allowed_svtypes: frozenset[str],
) -> TableInput:
    """Validate an internally normalized CNV frame and add derived columns."""
    frame = frame.copy().reset_index(drop=True)
    if frame.empty:
        raise ToolError("CNV input contains no supported records")

    if single_sample:
        frame["Run"] = "single_sample"
    else:
        if "Run" not in frame:
            raise ToolError("The CNV input does not contain sample identifiers")
        frame["Run"] = frame["Run"].astype("string").str.strip()
        if frame["Run"].isna().any() or (frame["Run"] == "").any():
            raise ToolError("The sample_id data contain an empty value")

    frame["chrom"] = frame["chrom"].astype("string").str.strip()
    if frame["chrom"].isna().any() or (frame["chrom"] == "").any():
        raise ToolError("The chr data contain an empty value")

    frame["SVTYPE"] = frame["SVTYPE"].astype("string").str.upper().str.strip()
    unsupported = ~frame["SVTYPE"].isin(allowed_svtypes)
    if unsupported.any():
        values = ", ".join(sorted(frame.loc[unsupported, "SVTYPE"].astype(str).unique()))
        raise ToolError(f"Unsupported SVTYPE value(s): {values}")

    for column in ("start", "stop"):
        converted = pd.to_numeric(frame[column], errors="coerce")
        invalid = converted.isna() | (converted % 1 != 0)
        if invalid.any():
            raise ToolError(f"The {column} data contain invalid coordinates")
        frame[column] = converted.astype(np.int64)

    invalid_coordinates = (frame["start"] < 1) | (frame["stop"] < frame["start"])
    if invalid_coordinates.any():
        examples = frame.loc[
            invalid_coordinates, ["Run", "chrom", "start", "stop"]
        ].head(5)
        raise ToolError("Invalid CNV coordinates:\n" + examples.to_string(index=False))

    if "length" in frame:
        lengths = pd.to_numeric(frame["length"], errors="coerce")
        invalid_lengths = lengths.isna() | (lengths % 1 != 0) | (lengths < 1)
        if invalid_lengths.any():
            raise ToolError("The CNV length data contain invalid values")
        frame["length"] = lengths.astype(np.int64)
    else:
        frame["length"] = np.maximum(1, frame["stop"] - frame["start"])

    frame["ID"] = [f"cnv_{number}" for number in range(1, len(frame) + 1)]
    return TableInput(
        frame=frame,
        single_sample=single_sample,
        base_columns=_output_base_columns(single_sample),
    )


def detect_cnv_input_type(path: Path) -> str:
    """Detect a manifest, standard table, native ClinCNV table, or Manta VCF."""
    try:
        with open_text_auto(path) as handle:
            lines = []
            for _, line in zip(range(512), handle):
                lines.append(line.rstrip("\n"))
    except (OSError, UnicodeError) as exc:
        raise ToolError(f"Cannot inspect CNV input {path}: {exc}") from exc

    nonempty = [line for line in lines if line.strip()]
    if not nonempty:
        raise ToolError("CNV input is empty")
    first_line = nonempty[0].lower()
    if (
        "sample_id" in first_line
        and "cnvs" in first_line
        and "vcf" in first_line
        and any(name in first_line for name in ("bam/cram", "bam_cram", "alignment"))
    ):
        return "manifest"
    if nonempty[0].startswith("##fileformat=VCF"):
        metadata = "\n".join(line for line in nonempty if line.startswith("##"))
        manta_record = any(
            len(fields := line.split("\t")) >= 3 and fields[2].startswith("Manta")
            for line in nonempty
            if not line.startswith("#")
        )
        metadata_lower = metadata.lower()
        if (
            "manta" in metadata_lower
            or "generatesvcandidates" in metadata_lower
            or manta_record
        ):
            return "manta"
        raise ToolError(
            "The CNV input is a VCF, but its header does not identify it as Manta output"
        )
    lowered = "\n".join(nonempty).lower()
    if "clincnv" in lowered or any(
        line.lstrip("#").lower().startswith("chr\tstart\tend\tcn_change")
        for line in nonempty
    ):
        return "clincnv"
    return "table"


def _parse_info_field(raw_info: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for token in raw_info.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            info[key] = value
        else:
            info[token] = ""
    return info


def _called_manta_samples(fields: list[str], sample_names: list[str]) -> list[str]:
    if not sample_names or len(fields) < 10:
        return sample_names
    format_keys = fields[8].split(":")
    if "GT" not in format_keys:
        return sample_names
    gt_offset = format_keys.index("GT")
    called: list[str] = []
    has_defined_genotype = False
    for sample_name, sample_field in zip(sample_names, fields[9:]):
        values = sample_field.split(":")
        genotype = values[gt_offset] if gt_offset < len(values) else "."
        alleles = re.split(r"[/|]", genotype)
        if any(allele not in {".", ""} for allele in alleles):
            has_defined_genotype = True
        if any(allele not in {"0", ".", ""} for allele in alleles):
            called.append(sample_name)
    # Candidate and somatic Manta outputs can carry missing GT values. In that
    # case the record itself is still a call and must not disappear.
    return called if has_defined_genotype else sample_names


def read_manta_vcf(
    path: Path,
    *,
    single_vcf: bool,
    allowed_svtypes: frozenset[str],
) -> TableInput:
    vprint(f"Reading Manta CNV calls from {path}")
    records: list[dict[str, object]] = []
    sample_names: list[str] = []
    fallback_sample = re.sub(r"\.vcf(?:\.bgz|\.gz)?$", "", path.name, flags=re.I)
    try:
        with open_text_auto(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if line_number % 100_000 == 0:
                    vprint(
                        f"Reading Manta VCF {path}: scanned "
                        f"{line_number:,} line(s)"
                    )
                if not line.strip() or line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                    sample_names = header[9:]
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ToolError(f"Manta VCF line {line_number} has fewer than 8 columns")
                try:
                    start = int(fields[1])
                except ValueError as exc:
                    raise ToolError(
                        f"Manta VCF line {line_number} has an invalid POS"
                    ) from exc
                info = _parse_info_field(fields[7])
                svtype = info.get("SVTYPE", "").upper().split(":", 1)[0]
                if svtype not in allowed_svtypes:
                    # Manta represents inversions as BND by default; neither BND nor
                    # INV is a class supported by the validator models.
                    continue
                raw_end = info.get("END")
                try:
                    end = int(raw_end) if raw_end not in {None, "", "."} else start
                except ValueError as exc:
                    raise ToolError(
                        f"Manta VCF line {line_number} has an invalid END"
                    ) from exc
                if svtype == "DEL" and raw_end in {None, "", "."}:
                    end = start + max(0, len(fields[3]) - 1)

                called_samples = _called_manta_samples(fields, sample_names)
                if single_vcf:
                    if sample_names and not called_samples:
                        continue
                    records.append(
                        {"chrom": fields[0], "start": start, "stop": end, "SVTYPE": svtype}
                    )
                else:
                    targets = called_samples or ([] if sample_names else [fallback_sample])
                    for sample_name in targets:
                        records.append(
                            {
                                "Run": sample_name,
                                "chrom": fields[0],
                                "start": start,
                                "stop": end,
                                "SVTYPE": svtype,
                            }
                        )
    except OSError as exc:
        raise ToolError(f"Cannot read Manta VCF {path}: {exc}") from exc

    table = _finalize_cnv_frame(
        pd.DataFrame.from_records(records),
        single_sample=single_vcf,
        allowed_svtypes=allowed_svtypes,
    )
    vprint(f"Loaded {len(table.frame)} supported CNV call(s) from Manta VCF {path}")
    return table


def read_clincnv_table(
    path: Path,
    *,
    single_vcf: bool,
    allowed_svtypes: frozenset[str],
) -> TableInput:
    vprint(f"Reading ClinCNV calls from {path}")
    header: str | None = None
    data_lines: list[str] = []
    metadata_sample: str | None = None
    sample_gender: str | None = None
    try:
        with open_text_auto(path) as handle:
            for line in handle:
                stripped = line.rstrip("\n")
                if stripped.startswith("##"):
                    match = re.match(
                        r"##(?:sample|sample_id)\s*[=:]\s*(.+)$",
                        stripped,
                        flags=re.IGNORECASE,
                    )
                    if match:
                        metadata_sample = match.group(1).strip()
                    gender_match = re.match(
                        r"##gender\s+of\s+sample\s*:\s*([FM])\s*$",
                        stripped,
                        flags=re.IGNORECASE,
                    )
                    if gender_match:
                        sample_gender = gender_match.group(1).upper()
                    continue
                if stripped.startswith("#") and "CN_change" in stripped:
                    header = stripped.lstrip("#")
                    continue
                if stripped and not stripped.startswith("#"):
                    data_lines.append(stripped)
    except OSError as exc:
        raise ToolError(f"Cannot read ClinCNV table {path}: {exc}") from exc

    if header is None:
        raise ToolError("ClinCNV input does not contain its #chr/.../CN_change header")
    try:
        raw = pd.read_csv(
            StringIO("\n".join([header, *data_lines])),
            sep="\t",
            dtype=str,
        )
    except Exception as exc:
        raise ToolError(f"Cannot parse ClinCNV table {path}: {exc}") from exc

    columns = {str(column).strip().lower(): str(column) for column in raw.columns}
    required = {"chr", "start", "end", "cn_change"}
    missing = sorted(required - set(columns))
    if missing:
        raise ToolError("ClinCNV input is missing column(s): " + ", ".join(missing))

    copy_number = pd.to_numeric(raw[columns["cn_change"]], errors="coerce")
    if copy_number.isna().any():
        raise ToolError("ClinCNV CN_change contains a non-numeric value")
    canonical_chromosomes = raw[columns["chr"]].map(canonical_chrom)
    expected_copy_number = pd.Series(2.0, index=raw.index)
    if sample_gender == "M":
        expected_copy_number.loc[canonical_chromosomes.isin({"X", "Y"})] = 1.0
    elif sample_gender == "F":
        # Female X follows the autosomal baseline; absence of Y is normal.
        expected_copy_number.loc[canonical_chromosomes == "Y"] = 0.0
    elif canonical_chromosomes.isin({"X", "Y"}).any():
        eprint(
            f"ClinCNV input {path} has sex-chromosome calls but no recognized "
            "'##gender of sample: F|M' metadata; using diploid baseline 2"
        )
    svtype = np.where(
        copy_number < expected_copy_number,
        "DEL",
        np.where(copy_number > expected_copy_number, "DUP", "DIP"),
    )
    frame = pd.DataFrame(
        {
            "chrom": raw[columns["chr"]],
            "start": raw[columns["start"]],
            "stop": raw[columns["end"]],
            "SVTYPE": svtype,
        }
    )
    if not single_vcf:
        sample_column = next(
            (
                columns[name]
                for name in ("sample_id", "sample", "sampleid")
                if name in columns
            ),
            None,
        )
        if sample_column is not None:
            frame["Run"] = raw[sample_column]
        elif metadata_sample:
            frame["Run"] = metadata_sample
        else:
            raise ToolError(
                "ClinCNV output has no sample identifier. Use --vcf-path for "
                "single-sample mode, or use the sample_id/cnvs/bam-or-cram/vcf "
                "manifest for multi-sample mode."
            )

    table = _finalize_cnv_frame(
        frame,
        single_sample=single_vcf,
        allowed_svtypes=allowed_svtypes,
    )
    vprint(f"Loaded {len(table.frame)} CNV call(s) from ClinCNV table {path}")
    return table


def _manifest_path(value: object, base_dir: Path) -> Path:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "none"}:
        raise ToolError("manifest contains an empty path")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def read_multi_sample_manifest(
    path: Path,
    sep: str | None,
    *,
    allowed_svtypes: frozenset[str],
) -> TableInput:
    """Expand sample_id/cnvs/bam-or-cram/vcf rows into normalized CNV rows."""
    vprint(f"Reading multi-sample manifest {path}")
    try:
        header_line = next(
            line for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()
        )
        if sep is None and "|" in header_line:
            raw = pd.read_csv(
                path,
                sep=r"\s*\|\s*",
                engine="python",
                dtype=str,
                encoding="utf-8-sig",
            )
        elif sep is None:
            raw = pd.read_csv(
                path,
                sep=None,
                engine="python",
                dtype=str,
                encoding="utf-8-sig",
            )
        else:
            raw = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8-sig")
    except Exception as exc:
        raise ToolError(f"Cannot read multi-sample manifest {path}: {exc}") from exc

    raw = raw.dropna(axis=1, how="all")
    raw.columns = [str(column).strip() for column in raw.columns]
    columns = {column.lower(): column for column in raw.columns}
    aliases = {
        "sample_id": ("sample_id",),
        "cnvs": ("cnvs",),
        "alignment": ("bam/cram", "bam_cram", "alignment", "alignment_path"),
        "vcf": ("vcf", "vcf_path"),
    }
    selected: dict[str, str] = {}
    for field, names in aliases.items():
        selected_column = next((columns[name] for name in names if name in columns), None)
        if selected_column is None:
            raise ToolError(
                f"Multi-sample manifest is missing {field}; expected columns: "
                "sample_id, cnvs, bam/cram, vcf"
            )
        selected[field] = selected_column

    if raw.empty:
        raise ToolError("Multi-sample manifest is empty")
    sample_ids = raw[selected["sample_id"]].astype("string").str.strip()
    if sample_ids.isna().any() or (sample_ids == "").any():
        raise ToolError("Multi-sample manifest contains an empty sample_id")
    if sample_ids.duplicated().any():
        duplicates = ", ".join(sample_ids[sample_ids.duplicated()].astype(str).unique())
        raise ToolError(f"Multi-sample manifest contains duplicate sample_id values: {duplicates}")

    vprint(f"Manifest contains {len(raw)} sample row(s)")

    frames: list[pd.DataFrame] = []
    base_dir = path.parent
    for row_number, (_, row) in enumerate(raw.iterrows(), start=2):
        sample_id = str(row[selected["sample_id"]]).strip()
        vprint(
            f"[{sample_id}] processing manifest row {row_number} "
            f"({row_number - 1}/{len(raw)})"
        )
        try:
            cnv_path = _manifest_path(row[selected["cnvs"]], base_dir)
            alignment_path = _manifest_path(row[selected["alignment"]], base_dir)
            vcf_path = _manifest_path(row[selected["vcf"]], base_dir)
            if not cnv_path.is_file():
                raise ToolError(f"CNV caller output does not exist: {cnv_path}")
            vprint(f"[{sample_id}] inspecting CNV caller output {cnv_path}")
            caller_type = detect_cnv_input_type(cnv_path)
            vprint(f"[{sample_id}] detected caller format: {caller_type}")
            if caller_type == "clincnv":
                sample_table = read_clincnv_table(
                    cnv_path,
                    single_vcf=True,
                    allowed_svtypes=allowed_svtypes,
                )
            elif caller_type == "manta":
                sample_table = read_manta_vcf(
                    cnv_path,
                    single_vcf=True,
                    allowed_svtypes=allowed_svtypes,
                )
            else:
                raise ToolError(
                    f"cnvs must point to ClinCNV or Manta output, detected {caller_type}"
                )
        except ToolError as exc:
            eprint(f"[{sample_id}] skipping manifest row {row_number}: {exc}")
            continue

        sample_frame = sample_table.frame.copy()
        sample_frame["Run"] = sample_id
        sample_frame["_vcf_path"] = str(vcf_path)
        sample_frame["_alignment_path"] = str(alignment_path)
        frames.append(sample_frame)
        vprint(
            f"[{sample_id}] accepted {len(sample_frame)} CNV call(s); "
            f"VCF={vcf_path}; alignment={alignment_path}"
        )

    if not frames:
        raise ToolError("No usable samples remain in the multi-sample manifest")
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop(columns=["ID"], errors="ignore")
    table = _finalize_cnv_frame(
        combined,
        single_sample=False,
        allowed_svtypes=allowed_svtypes,
    )
    vprint(
        f"Finished manifest expansion: {table.frame['Run'].nunique()} sample(s), "
        f"{len(table.frame)} CNV call(s)"
    )
    return table


def read_table(
    path: Path,
    sep: str | None,
    *,
    single_vcf: bool = False,
    allowed_svtypes: frozenset[str] = VALIDATION_SVTYPES,
) -> TableInput:
    """Auto-detect and read manifest, standard, ClinCNV, or Manta CNV input."""
    if not path.exists():
        raise ToolError(f"CNV input does not exist: {path}")

    vprint(f"Inspecting CNV input {path}")
    input_type = detect_cnv_input_type(path)
    vprint(f"Detected CNV input format: {input_type}")
    if input_type == "manifest":
        if single_vcf:
            raise ToolError(
                "A multi-sample manifest already defines one VCF per sample; "
                "do not combine it with --vcf-path"
            )
        return read_multi_sample_manifest(
            path,
            sep,
            allowed_svtypes=allowed_svtypes,
        )
    if input_type == "manta":
        return read_manta_vcf(
            path,
            single_vcf=True,
            allowed_svtypes=allowed_svtypes,
        )
    if input_type == "clincnv":
        return read_clincnv_table(
            path,
            single_vcf=True,
            allowed_svtypes=allowed_svtypes,
        )

    try:
        if sep is None:
            raw = pd.read_csv(path, sep=None, engine="python", dtype=str)
        else:
            raw = pd.read_csv(path, sep=sep, dtype=str)
    except Exception as exc:
        raise ToolError(f"Cannot read CNV table {path}: {exc}") from exc

    if raw.empty:
        raise ToolError("CNV table is empty")

    # New layouts are positional and header-name independent:
    #   multi:  sample_id, chr, start, end, SVTYPE
    #   single: chr, start, end, SVTYPE
    # The former legacy layouts with a length column before SVTYPE are accepted.
    layouts = {
        "multi": (
            raw.shape[1] >= 5
            and _column_is_numeric(raw, 2)
            and _column_is_numeric(raw, 3)
            and _column_is_svtype(raw, 4, allowed_svtypes)
        ),
        "single": (
            raw.shape[1] >= 4
            and _column_is_numeric(raw, 1)
            and _column_is_numeric(raw, 2)
            and _column_is_svtype(raw, 3, allowed_svtypes)
        ),
        "legacy_multi": (
            raw.shape[1] >= 6
            and _column_is_numeric(raw, 2)
            and _column_is_numeric(raw, 3)
            and _column_is_numeric(raw, 4)
            and _column_is_svtype(raw, 5, allowed_svtypes)
        ),
        "legacy_single": (
            raw.shape[1] >= 5
            and _column_is_numeric(raw, 1)
            and _column_is_numeric(raw, 2)
            and _column_is_numeric(raw, 3)
            and _column_is_svtype(raw, 4, allowed_svtypes)
        ),
    }

    legacy_single_likely = (
        layouts["legacy_single"] and _column_looks_like_chrom(raw, 0)
    )
    if single_vcf and layouts["single"]:
        layout = "single"
    elif single_vcf and legacy_single_likely:
        layout = "legacy_single"
    elif layouts["multi"]:
        layout = "multi"
    elif layouts["legacy_multi"]:
        layout = "legacy_multi"
    elif layouts["single"]:
        layout = "single"
    elif layouts["legacy_single"]:
        layout = "legacy_single"
    else:
        raise ToolError(
            "Cannot identify required positional CNV data. Expected either "
            "sample_id, chr, start, end, SVTYPE or chr, start, end, SVTYPE."
        )

    input_has_sample = layout in {"multi", "legacy_multi"}
    if input_has_sample:
        sample_position, chrom_position, start_position, stop_position = 0, 1, 2, 3
        svtype_position = 4 if layout == "multi" else 5
        length_position = None if layout == "multi" else 4
    else:
        sample_position = None
        chrom_position, start_position, stop_position = 0, 1, 2
        svtype_position = 3 if layout == "single" else 4
        length_position = None if layout == "single" else 3

    frame = pd.DataFrame(index=raw.index)
    if not input_has_sample:
        frame["Run"] = "single_sample"
    else:
        assert sample_position is not None
        sample_values = raw.iloc[:, sample_position]
        if sample_values.isna().any():
            raise ToolError("The leading sample_id column contains an empty value")
        frame["Run"] = sample_values.astype(str).str.strip()
        if (frame["Run"] == "").any():
            raise ToolError("The leading sample_id column contains an empty value")

    chrom_values = raw.iloc[:, chrom_position]
    if chrom_values.isna().any():
        raise ToolError("The positional chr data contain an empty value")
    frame["chrom"] = chrom_values.astype(str).str.strip()
    if (frame["chrom"] == "").any():
        raise ToolError("The positional chr data contain an empty value")
    frame["SVTYPE"] = raw.iloc[:, svtype_position]
    frame["start"] = raw.iloc[:, start_position]
    frame["stop"] = raw.iloc[:, stop_position]
    if length_position is not None:
        frame["length"] = raw.iloc[:, length_position]

    table = _finalize_cnv_frame(
        frame,
        single_sample=not input_has_sample,
        allowed_svtypes=allowed_svtypes,
    )
    vprint(
        f"Loaded standard CNV table using {layout} layout: "
        f"{len(table.frame)} row(s)"
    )
    return table

def open_text_auto(path: Path):
    if path.suffix in {".gz", ".bgz"} or path.name.endswith((".vcf.gz", ".vcf.bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_numeric_values(
    raw: str | None,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> list[float]:
    if raw is None:
        return []
    raw = raw.strip()
    if not raw or raw == ".":
        return []
    values: list[float] = []
    for token in raw.split(","):
        try:
            value = float(token)
        except ValueError:
            continue
        if not math.isfinite(value):
            continue
        if minimum is not None and value < minimum:
            continue
        if maximum is not None and value > maximum:
            continue
        values.append(value)
    return values


def parse_vaf_values(raw: str | None) -> list[float]:
    return parse_numeric_values(raw, minimum=0.0, maximum=1.0)


def vaf_from_sample(format_keys: list[str], sample_value: str) -> float | None:
    fields = sample_value.split(":")
    data = dict(zip(format_keys, fields))

    for tag in ("VAF", "AF"):
        values = parse_vaf_values(data.get(tag))
        if values:
            return values[0]

    ad = parse_numeric_values(data.get("AD"), minimum=0.0)
    if len(ad) >= 2:
        total = sum(ad)
        if total > 0:
            return ad[1] / total
    return None


def choose_vcf_sample(sample_names: list[str], run: str) -> int:
    if not sample_names:
        raise ToolError("VCF has no sample columns")
    if run in sample_names:
        return sample_names.index(run)
    if len(sample_names) == 1:
        return 0
    raise ToolError(
        f"VCF contains multiple samples ({', '.join(sample_names[:8])}) and none is named {run}"
    )


def build_assigners(
    rows: pd.DataFrame,
    values: dict[int, list[float]],
) -> dict[str, SweepAssigner]:
    by_chrom: dict[str, list[Interval]] = {}
    for row_index, row in rows.iterrows():
        chrom = canonical_chrom(row["chrom"])
        by_chrom.setdefault(chrom, []).append(
            Interval(row_index=int(row_index), start=int(row["start"]), stop=int(row["stop"]))
        )
    return {chrom: SweepAssigner(intervals, values) for chrom, intervals in by_chrom.items()}


def parse_vcf_python(
    vcf_path: Path,
    run: str,
    rows: pd.DataFrame,
    values: dict[int, list[float]],
    qual_cutoff: float | None,
    filter_name: str | None,
) -> int:
    assigners = build_assigners(rows, values)
    sample_names: list[str] | None = None
    sample_column = -1
    accepted = 0
    scanned_records = 0

    try:
        handle = open_text_auto(vcf_path)
    except OSError as exc:
        raise ToolError(f"Cannot open VCF {vcf_path}: {exc}") from exc

    with handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                sample_column = 9 + choose_vcf_sample(sample_names, run)
                continue
            if line.startswith("#") or not line.strip():
                continue
            if sample_names is None:
                raise ToolError(f"VCF header line #CHROM was not found in {vcf_path}")

            scanned_records += 1
            if scanned_records % 100_000 == 0:
                vprint(
                    f"[{run}] Python VCF scan: {scanned_records:,} record(s), "
                    f"{accepted:,} accepted"
                )

            columns = line.rstrip("\n").split("\t")
            if len(columns) <= sample_column or len(columns) < 10:
                continue

            chrom = canonical_chrom(columns[0])
            assigner = assigners.get(chrom)
            if assigner is None:
                continue

            try:
                position = int(columns[1])
            except ValueError:
                continue

            if qual_cutoff is not None:
                try:
                    quality = float(columns[5])
                except ValueError:
                    continue
                if not quality > qual_cutoff:
                    continue

            if filter_name is not None:
                filters = columns[6].split(";")
                if filter_name not in filters:
                    continue

            format_keys = columns[8].split(":")
            vaf = vaf_from_sample(format_keys, columns[sample_column])
            if vaf is None:
                continue

            assigner.add(position, vaf)
            accepted += 1

    return accepted


def run_command(command: Sequence[str], *, text: bool = True) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=text,
        )
    except OSError as exc:
        raise ToolError(f"Cannot run {' '.join(command)}: {exc}") from exc


def bcftools_metadata(
    vcf_path: Path,
    threads: int,
) -> tuple[list[str], list[str], set[str]]:
    samples_result = run_command(["bcftools", "query", "-l", str(vcf_path)])
    if samples_result.returncode != 0:
        raise ToolError(samples_result.stderr.strip() or "bcftools query -l failed")
    samples = [line.strip() for line in samples_result.stdout.splitlines() if line.strip()]

    header_result = run_command(
        ["bcftools", "view", "--threads", str(threads), "-h", str(vcf_path)]
    )
    if header_result.returncode != 0:
        raise ToolError(header_result.stderr.strip() or "bcftools view -h failed")

    contigs: list[str] = []
    format_tags: set[str] = set()
    for line in header_result.stdout.splitlines():
        contig_match = re.match(r"##contig=<ID=([^,>]+)", line)
        if contig_match:
            contigs.append(contig_match.group(1))
        format_match = re.match(r"##FORMAT=<ID=([^,>]+)", line)
        if format_match:
            format_tags.add(format_match.group(1))
    return samples, contigs, format_tags


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, stop in ordered[1:]:
        current = merged[-1]
        if start <= current[1] + 1:
            current[1] = max(current[1], stop)
        else:
            merged.append([start, stop])
    return [(start, stop) for start, stop in merged]


def parse_bcftools_line(line: str, tags: list[str]) -> tuple[str, int, float] | None:
    columns = line.rstrip("\n").split("\t")
    if len(columns) < 2 + len(tags):
        return None
    try:
        position = int(columns[1])
    except ValueError:
        return None

    data = dict(zip(tags, columns[2:]))
    for tag in ("VAF", "AF"):
        values = parse_vaf_values(data.get(tag))
        if values:
            return columns[0], position, values[0]

    ad = parse_numeric_values(data.get("AD"), minimum=0.0)
    if len(ad) >= 2:
        total = sum(ad)
        if total > 0:
            return columns[0], position, ad[1] / total
    return None


def parse_vcf_bcftools(
    vcf_path: Path,
    run: str,
    rows: pd.DataFrame,
    values: dict[int, list[float]],
    qual_cutoff: float | None,
    filter_name: str | None,
    threads: int,
) -> int:
    samples, contigs, format_tags = bcftools_metadata(vcf_path, threads)
    selected_sample = samples[choose_vcf_sample(samples, run)]

    tags = [tag for tag in ("VAF", "AF", "AD") if tag in format_tags]
    if not tags:
        raise ToolError("VCF contains none of the supported FORMAT fields: VAF, AF, AD")

    canonical_to_contig: dict[str, str] = {}
    for contig in contigs:
        canonical_to_contig.setdefault(canonical_chrom(contig), contig)

    assigners = build_assigners(rows, values)
    bed_intervals: dict[str, list[tuple[int, int]]] = {}
    for _, row in rows.iterrows():
        canonical = canonical_chrom(row["chrom"])
        actual_contig = canonical_to_contig.get(canonical)
        if actual_contig is None:
            continue
        bed_intervals.setdefault(actual_contig, []).append((int(row["start"]), int(row["stop"])))

    if not bed_intervals:
        return 0

    with tempfile.NamedTemporaryFile("wt", suffix=".bed", delete=False) as bed_handle:
        bed_path = Path(bed_handle.name)
        for contig in contigs:
            for start, stop in merge_intervals(bed_intervals.get(contig, [])):
                # BED is 0-based, half-open; CNV table is interpreted as 1-based, inclusive.
                bed_handle.write(f"{contig}\t{start - 1}\t{stop}\n")

    format_string = "%CHROM\t%POS[" + "".join(f"\t%{tag}" for tag in tags) + "]\n"
    command = [
        "bcftools", "query",
        "-R", str(bed_path), "-s", selected_sample,
    ]

    expressions: list[str] = []
    if qual_cutoff is not None:
        expressions.append(f"QUAL>{qual_cutoff:g}")
    if filter_name is not None:
        escaped = filter_name.replace('"', '\\"')
        expressions.append(f'FILTER="{escaped}"')
    if expressions:
        command.extend(["-i", " && ".join(expressions)])
    command.extend(["-f", format_string, str(vcf_path)])

    accepted = 0
    scanned_records = 0
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        bed_path.unlink(missing_ok=True)
        raise ToolError(f"Cannot run bcftools: {exc}") from exc

    assert process.stdout is not None
    try:
        for line in process.stdout:
            scanned_records += 1
            if scanned_records % 100_000 == 0:
                vprint(
                    f"[{run}] bcftools VCF scan: {scanned_records:,} record(s), "
                    f"{accepted:,} accepted"
                )
            parsed = parse_bcftools_line(line, tags)
            if parsed is None:
                continue
            chrom, position, vaf = parsed
            assigner = assigners.get(canonical_chrom(chrom))
            if assigner is None:
                continue
            assigner.add(position, vaf)
            accepted += 1
    finally:
        process.stdout.close()

    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    bed_path.unlink(missing_ok=True)
    if return_code != 0:
        raise ToolError(stderr.strip() or f"bcftools query failed with exit code {return_code}")
    return accepted


def bandwidth_nrd0(values: np.ndarray) -> float:
    """R-compatible bw.nrd0 bandwidth used by stats::density."""
    n = values.size
    if n < 2:
        raise ToolError("At least two VAF values are required for KDE")

    standard_deviation = float(np.std(values, ddof=1))
    q25, q75 = np.percentile(values, [25, 75], method="linear")
    iqr = float(q75 - q25)
    low = min(standard_deviation, iqr / 1.34)
    numerical_zero = 1e-12
    if not math.isfinite(low) or low <= numerical_zero:
        if math.isfinite(standard_deviation) and standard_deviation > numerical_zero:
            low = standard_deviation
        elif abs(float(values[0])) > numerical_zero:
            low = abs(float(values[0]))
        else:
            low = 1.0
    return max(0.9 * low * n ** (-0.2), 1e-8)


def gaussian_kde_grid(values: Sequence[float], grid_size: int = 512) -> tuple[np.ndarray, np.ndarray]:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data) & (data >= 0.0) & (data <= 1.0)]
    if data.size < 2:
        raise ToolError("At least two valid VAF values are required for KDE")

    bandwidth = bandwidth_nrd0(data)
    grid = np.linspace(0.0, 1.0, grid_size)
    density = np.zeros(grid_size, dtype=float)

    normalizer = 1.0 / (data.size * bandwidth * math.sqrt(2.0 * math.pi))
    chunk_size = 20_000
    for offset in range(0, data.size, chunk_size):
        chunk = data[offset : offset + chunk_size]
        z = (grid[:, None] - chunk[None, :]) / bandwidth
        density += np.exp(-0.5 * z * z).sum(axis=1)
    density *= normalizer

    total_area = float(np.trapezoid(density, grid))
    if not math.isfinite(total_area) or total_area <= 0:
        raise ToolError("KDE has zero or invalid area on [0, 1]")
    density /= total_area
    return grid, density


def integrate_linear_density(
    grid: np.ndarray,
    density: np.ndarray,
    lower: float,
    upper: float,
) -> float:
    lower = max(0.0, lower)
    upper = min(1.0, upper)
    if upper <= lower:
        return 0.0

    internal = (grid > lower) & (grid < upper)
    x = np.concatenate(([lower], grid[internal], [upper]))
    y = np.concatenate(
        ([np.interp(lower, grid, density)], density[internal], [np.interp(upper, grid, density)])
    )
    return float(np.trapezoid(y, x))


def calculate_features(values: Sequence[float]) -> dict[str, float]:
    grid, density = gaussian_kde_grid(values, grid_size=512)
    features = {
        name: integrate_linear_density(grid, density, lower, upper)
        for name, (lower, upper) in DEFAULT_INTERVALS.items()
    }
    features["X0.33_X0.66"] = features["X0.33"] + features["X0.66"]
    return features


def bundled_asset_candidates(*parts: str) -> tuple[Path, ...]:
    """Return development-tree and installed-prefix candidates for an asset."""
    script_path = Path(__file__).resolve()
    relative = Path(*parts)
    return (
        script_path.parent / relative,
        script_path.parent.parent / "share" / "vacca" / relative,
        Path(sys.prefix) / "share" / "vacca" / relative,
    )


def resolve_vaf_model() -> Path:
    """Find the bundled fitted VAF SVC."""
    for path in bundled_asset_candidates("models", VAF_MODEL_FILENAME):
        if path.is_file():
            return path
    raise ToolError(
        f"Bundled VAF model models/{VAF_MODEL_FILENAME} was not found in the "
        "project or installation share directory."
    )


def load_vaf_model(model_path: Path) -> Any:
    """Load and validate the fitted VAF SVC."""
    try:
        import joblib
        from sklearn.svm import SVC
    except ImportError as exc:
        raise ToolError(
            "VAF validation requires joblib and scikit-learn. Install the "
            "dependencies from environment.yml."
        ) from exc

    try:
        model = joblib.load(model_path)
    except Exception as exc:
        raise ToolError(f"Cannot load VAF model {model_path}: {exc}") from exc

    if not isinstance(model, SVC):
        raise ToolError(
            f"VAF model {model_path} contains {type(model).__name__}, expected SVC"
        )

    parameters = model.get_params(deep=False)
    mismatches = [
        f"{name}={parameters.get(name)!r} (expected {expected!r})"
        for name, expected in VAF_MODEL_PARAMETERS.items()
        if parameters.get(name) != expected
    ]
    if mismatches:
        raise ToolError(
            f"VAF model {model_path} has unexpected parameters: " + ", ".join(mismatches)
        )

    feature_count = int(getattr(model, "n_features_in_", -1))
    if feature_count != len(MODEL_FEATURE_COLUMNS):
        raise ToolError(
            f"VAF model expects {feature_count} feature(s), but the validator provides "
            f"{len(MODEL_FEATURE_COLUMNS)}"
        )
    feature_names = getattr(model, "feature_names_in_", None)
    if feature_names is not None and tuple(map(str, feature_names)) != tuple(MODEL_FEATURE_COLUMNS):
        raise ToolError(
            f"VAF model feature order is {tuple(map(str, feature_names))!r}, expected "
            f"{tuple(MODEL_FEATURE_COLUMNS)!r}"
        )

    try:
        encoded_classes = tuple(map(int, np.asarray(model.classes_).tolist()))
    except Exception as exc:
        raise ToolError(f"Cannot read VAF model classes: {exc}") from exc
    expected_encoded_classes = tuple(range(len(VAF_MODEL_CLASSES)))
    if encoded_classes != expected_encoded_classes:
        raise ToolError(
            f"VAF model classes are {encoded_classes!r}, expected "
            f"{expected_encoded_classes!r} for DEL, DIP, DUP"
        )
    if not bool(getattr(model, "probability", False)):
        raise ToolError(
            "VAF SVC was fitted without probability=True; MODEL_SCORE columns "
            "require predict_proba()"
        )
    return model


def model_predict(
    model: Any,
    feature_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Predict SVC classes and probabilities for ordered KDE features."""
    feature_frame = pd.DataFrame(feature_rows, columns=MODEL_FEATURE_COLUMNS)
    try:
        encoded_predictions = np.asarray(model.predict(feature_frame), dtype=int)
        probabilities = np.asarray(model.predict_proba(feature_frame), dtype=float)
        predictions = np.asarray(
            [VAF_MODEL_CLASSES[int(value)] for value in encoded_predictions],
            dtype=object,
        )
    except Exception as exc:
        raise ToolError(f"VAF model inference failed: {exc}") from exc

    try:
        class_names = tuple(
            VAF_MODEL_CLASSES[int(value)]
            for value in np.asarray(model.classes_, dtype=int)
        )
    except (IndexError, TypeError, ValueError) as exc:
        raise ToolError(f"Cannot map VAF model classes: {exc}") from exc
    expected_shape = (len(feature_frame), len(class_names))
    if probabilities.shape != expected_shape:
        raise ToolError(
            f"VAF model returned probability shape {probabilities.shape}, expected "
            f"{expected_shape}"
        )
    if not np.isfinite(probabilities).all():
        raise ToolError("VAF model returned non-finite probabilities")
    return predictions, probabilities, class_names


def resolve_vaf_train() -> Path:
    """Find the bundled VAF training feature table."""
    for path in bundled_asset_candidates("data", VAF_TRAIN_FILENAME):
        if path.is_file():
            return path
    raise ToolError(
        f"Bundled VAF training table data/{VAF_TRAIN_FILENAME} was not found in "
        "the project or installation share directory."
    )


def load_vaf_train(path: Path) -> pd.DataFrame:
    """Load and validate the feature table used for train-space plots."""
    try:
        train = pd.read_csv(path, sep="\t", dtype={"SVTYPE": "string"})
    except Exception as exc:
        raise ToolError(f"Cannot read VAF training table {path}: {exc}") from exc

    required_columns = {"SVTYPE", "X0.5", "X1", "X0.33_X0.66"}
    missing_columns = sorted(required_columns.difference(train.columns))
    if missing_columns:
        raise ToolError(
            f"VAF training table {path} is missing column(s): "
            + ", ".join(missing_columns)
        )
    if train.empty:
        raise ToolError(f"VAF training table {path} is empty")

    train = train.loc[:, ["SVTYPE", "X0.5", "X1", "X0.33_X0.66"]].copy()
    train["SVTYPE"] = train["SVTYPE"].astype("string").str.upper().str.strip()
    invalid_classes = sorted(
        set(train["SVTYPE"].dropna().astype(str)).difference(VAF_MODEL_CLASSES)
    )
    if invalid_classes or train["SVTYPE"].isna().any():
        details = ", ".join(invalid_classes) if invalid_classes else "missing SVTYPE"
        raise ToolError(f"VAF training table {path} has invalid classes: {details}")

    for column in MODEL_FEATURE_COLUMNS:
        train[column] = pd.to_numeric(train[column], errors="coerce")
        values = train[column].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ToolError(
                f"VAF training table {path} contains non-numeric or non-finite "
                f"values in {column}"
            )
    return train


def render_vaf_train_dimension_plot(
    train: pd.DataFrame,
    *,
    target_x: float,
    target_y: float,
    target_label: str,
    output_path: Path,
    dpi: int,
) -> None:
    """Plot one CNV in the X0.33+X0.66 versus X0.5 train projection."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ToolError(
            "VAF train-space plotting requires matplotlib. Install the "
            "dependencies from environment.yml."
        ) from exc

    colors = {"DEL": "red", "DIP": "green", "DUP": "blue"}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8, 6), constrained_layout=True)
    try:
        for class_name in VAF_MODEL_CLASSES:
            class_rows = train.loc[train["SVTYPE"] == class_name]
            if class_rows.empty:
                continue
            axis.scatter(
                class_rows["X0.33_X0.66"],
                class_rows["X0.5"],
                c=colors[class_name],
                s=24,
                alpha=0.85,
                label=f"train_{class_name}",
                edgecolors="none",
            )
        axis.scatter(
            [target_x],
            [target_y],
            c="black",
            s=62,
            label="CNV",
            zorder=5,
        )
        axis.annotate(
            target_label,
            (target_x, target_y),
            xytext=(7, 7),
            textcoords="offset points",
            fontsize=8,
            color="black",
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "alpha": 0.8, "ec": "none"},
            zorder=6,
        )
        axis.set_xlabel("X0.33_X0.66 (X0.33 + X0.66)")
        axis.set_ylabel("X0.5")
        axis.grid(color="#d9d9d9", linewidth=0.6, alpha=0.8)
        axis.legend(title="SVTYPE", frameon=True)
        figure.savefig(
            output_path,
            format="jpeg",
            dpi=dpi,
            facecolor="white",
            edgecolor="white",
            transparent=False,
        )
    finally:
        plt.close(figure)


def run_vaf_train_dimension_plots(
    frame: pd.DataFrame,
    vaf_result: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.Series:
    """Create per-CNV VAF train-space JPEGs and return their absolute paths."""
    plot_paths = pd.Series(pd.NA, index=frame.index, dtype="string")
    train_path = resolve_vaf_train()
    vprint(f"Loading VAF training features from {train_path}")
    train = load_vaf_train(train_path)

    eligible_indices: list[int] = []
    for row_index in frame.index:
        try:
            target_x = float(vaf_result.at[row_index, "X0.33_X0.66"])
            target_y = float(vaf_result.at[row_index, "X0.5"])
        except (TypeError, ValueError):
            continue
        if math.isfinite(target_x) and math.isfinite(target_y):
            eligible_indices.append(int(row_index))
    vprint(
        f"VAF train-space plot preparation: {len(eligible_indices)}/{len(frame)} "
        "CNV row(s) have plottable features"
    )

    for completed, row_index in enumerate(eligible_indices, start=1):
        row = frame.loc[row_index]
        run = str(row["Run"])
        chrom = str(row["chrom"])
        cnv_start = int(row["start"])
        cnv_stop = int(row["stop"])
        output_dir = (
            resolve_sample_output_dir(args.output, run)
            / "vaf_train_dimension"
            / plot_svtype_folder_name(row["SVTYPE"])
        )
        target_label = f"{row['SVTYPE']}_{chrom}_{cnv_start}_{cnv_stop}"
        output_path = output_dir / (
            f"{safe_path_token(run)}_{safe_path_token(chrom)}_"
            f"{cnv_start}_{cnv_stop}.jpeg"
        )
        try:
            render_vaf_train_dimension_plot(
                train,
                target_x=float(vaf_result.at[row_index, "X0.33_X0.66"]),
                target_y=float(vaf_result.at[row_index, "X0.5"]),
                target_label=target_label,
                output_path=output_path,
                dpi=args.plot_dpi,
            )
        except Exception as exc:
            eprint(f"[{run}:{row['ID']}] VAF train-space plot error: {exc}")
            continue
        plot_paths.at[row_index] = str(output_path.resolve())
        vprint(
            f"[{run}:{row['ID']}] wrote VAF train-space plot "
            f"{output_path} ({completed}/{len(eligible_indices)})"
        )
    return plot_paths


def output_label(model_label: str) -> str:
    return model_label


def validation_code(original_type: str, predicted_model_label: str) -> str:
    """Compare submitted and predicted classes using public output semantics."""
    original_model_label = output_label(original_type)
    predicted_model_label = output_label(predicted_model_label)
    if predicted_model_label == original_model_label:
        return "CONSISTENT"
    if predicted_model_label == "DIP":
        return "NOT_CNV"
    return "INCONSISTENT"


def process_sample(
    run: str,
    rows: pd.DataFrame,
    values: dict[int, list[float]],
    args: argparse.Namespace,
) -> tuple[str, str]:
    if "_vcf_path" in rows:
        manifest_paths = rows["_vcf_path"].drop_duplicates().astype(str).tolist()
        if len(manifest_paths) != 1:
            return "vcf_error", f"sample {run} has multiple VCF paths in the manifest"
        vcf_path = Path(manifest_paths[0])
    elif args.vcf_path is not None:
        vcf_path = args.vcf_path.expanduser().resolve()
    else:
        vcf_path = make_vcf_path(args.vcf_dir, args.vcf_mask, run)
    if not vcf_path.exists():
        return "missing_vcf", f"VCF does not exist: {vcf_path}"

    backend = args.backend
    if backend == "auto":
        backend = "bcftools" if shutil.which("bcftools") else "python"
    vprint(
        f"[{run}] extracting VAF values from {vcf_path} with {backend} backend "
        f"for {len(rows)} CNV interval(s)"
    )

    try:
        if backend == "bcftools":
            if not shutil.which("bcftools"):
                raise ToolError("bcftools is not available in PATH")
            try:
                accepted = parse_vcf_bcftools(
                    vcf_path, run, rows, values, args.qual_cutoff, args.filter_name,
                    args.threads,
                )
            except ToolError as exc:
                if args.backend == "auto":
                    eprint(f"[{run}] bcftools failed; using Python fallback: {exc}")
                    accepted = parse_vcf_python(
                        vcf_path, run, rows, values, args.qual_cutoff, args.filter_name
                    )
                else:
                    raise
        else:
            accepted = parse_vcf_python(
                vcf_path, run, rows, values, args.qual_cutoff, args.filter_name
            )
    except (OSError, ToolError) as exc:
        return "vcf_error", str(exc)

    return "parsed", f"accepted {accepted} VCF records"


def safe_path_token(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    return text.strip("._-") or "unknown"


def plot_svtype_folder_name(svtype: object) -> str:
    svtype_text = str(svtype).upper().strip()
    if svtype_text in {"DEL", "DIP", "DUP"}:
        return svtype_text
    raise ToolError(f"Unsupported SVTYPE for coverage plot: {svtype!r}")


def plot_output_dir(plot_dir: Path, sample: str, svtype: object) -> Path:
    sample_dir = plot_dir / safe_path_token(sample)
    type_dir = sample_dir / plot_svtype_folder_name(svtype)
    type_dir.mkdir(parents=True, exist_ok=True)
    return type_dir


def coverage_plot_window(
    cnv_start: int,
    cnv_stop: int,
    cnv_length: int,
) -> tuple[int, int]:
    """Return the requested window using the CNV table's length value."""
    flank_multiplier = 3 if cnv_length >= 100_000 else 4
    flank = flank_multiplier * cnv_length
    return max(1, cnv_start - flank), cnv_stop + flank


def validate_alignment_path(path: Path) -> None:
    """Validate a coordinate-sorted BAM/CRAM input and its expected index."""
    suffix = path.suffix.lower()
    if suffix not in {".bam", ".cram"}:
        raise ToolError(f"Alignment file must have a .bam or .cram extension: {path}")
    if not path.is_file():
        raise ToolError(f"BAM/CRAM file does not exist: {path}")

    index_suffixes = (".bai", ".csi") if suffix == ".bam" else (".crai",)
    index_candidates = tuple(
        candidate
        for index_suffix in index_suffixes
        for candidate in (Path(f"{path}{index_suffix}"), path.with_suffix(index_suffix))
    )
    if not any(index_path.is_file() for index_path in index_candidates):
        expected = " or ".join(str(index_path) for index_path in index_candidates)
        raise ToolError(f"Alignment index was not found; expected {expected}")


def samtools_depth(
    alignment_path: Path,
    chrom: str,
    start: int,
    stop: int,
    threads: int,
) -> tuple[list[int], list[float]]:
    region = f"{chrom}:{start}-{stop}"
    command = [
        "samtools", "depth", "-@", str(threads), "-aa", "-r", region,
        str(alignment_path),
    ]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        raise ToolError(f"Cannot run samtools depth: {exc}") from exc

    positions: list[int] = []
    coverages: list[float] = []
    assert process.stdout is not None
    try:
        for line_number, line in enumerate(process.stdout, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                position = int(fields[1])
                coverage = float(fields[2])
            except ValueError as exc:
                process.kill()
                raise ToolError(
                    f"Malformed samtools depth output at line {line_number} for {region}"
                ) from exc
            if not math.isfinite(coverage) or coverage < 0:
                process.kill()
                raise ToolError(f"Invalid coverage {coverage!r} returned for {region}")
            positions.append(position)
            coverages.append(coverage)
    finally:
        process.stdout.close()

    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise ToolError(stderr.strip() or f"samtools depth failed for {region}")
    if len(positions) == 0:
        raise ToolError(f"samtools depth returned no positions for {region}")
    return positions, coverages


def render_coverage_plot(
    positions: Sequence[int],
    coverages: Sequence[float],
    *,
    chrom: str,
    view_start: int,
    view_stop: int,
    cnv_start: int,
    cnv_stop: int,
    output_path: Path,
    coverage_cap_multiplier: float,
    dpi: int,
    legend: bool = True,
    cnv_borders: bool = True,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, MaxNLocator
    except ImportError as exc:
        raise ToolError(
            "Plotting requires matplotlib. Install the dependencies from environment.yml."
        ) from exc

    coverage_array = np.asarray(coverages, dtype=float)
    median_coverage = float(np.median(coverage_array))
    maximum_coverage = median_coverage * coverage_cap_multiplier
    if not math.isfinite(maximum_coverage) or maximum_coverage <= 0:
        maximum_coverage = float(np.max(coverage_array))
    if not math.isfinite(maximum_coverage) or maximum_coverage <= 0:
        maximum_coverage = 1.0
    plotted_coverage = np.minimum(coverage_array, maximum_coverage)

    figure, axis = plt.subplots(figsize=(12, 4), constrained_layout=True)
    axis.fill_between(
        positions,
        plotted_coverage,
        0,
        color="black",
        linewidth=0,
        zorder=1,
    )
    if cnv_borders:
        axis.axvline(cnv_start, color="red", linewidth=2.0, zorder=3)
        axis.axvline(cnv_stop, color="red", linewidth=2.0, zorder=3)
    axis.set_xlim(view_start, view_stop)
    axis.set_ylim(0, maximum_coverage)
    axis.margins(x=0, y=0)

    if legend:
        axis.set_xlabel(f"{chrom} coordinate (kb)")
        axis.set_ylabel("Coverage")
        axis.xaxis.set_major_locator(MaxNLocator(nbins=7, min_n_ticks=4))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1000:,.0f}"))
        axis.grid(axis="x", color="#d9d9d9", linewidth=0.6, zorder=0)
    else:
        axis.set_axis_off()

    figure.savefig(
        output_path,
        format=output_path.suffix.lstrip("."),
        dpi=dpi,
        facecolor="white",
        edgecolor="white",
        transparent=False,
    )
    plt.close(figure)


def median_coverage_windows(
    positions: Sequence[int],
    coverages: Sequence[float],
    window_size: int,
) -> tuple[Sequence[int], Sequence[float]]:
    """Collapse per-base coverage to median coverage in genomic windows."""
    if window_size < 1:
        raise ValueError("window_size must be at least 1")
    if len(positions) != len(coverages):
        raise ValueError("positions and coverages must have the same length")
    if len(positions) == 0:
        raise ValueError("positions and coverages cannot be empty")

    coverage_array = np.asarray(coverages, dtype=float)
    position_count = len(positions)
    first_position = int(positions[0])
    last_position = int(positions[-1])

    # `samtools depth -aa` emits one row per reference base, so the usual path
    # can use a vectorized reshape without constructing another per-base
    # position array. Keep a genomic-bin fallback for defensive robustness.
    if last_position - first_position + 1 == position_count:
        full_window_count = position_count // window_size
        full_value_count = full_window_count * window_size
        median_parts: list[np.ndarray] = []
        position_parts: list[np.ndarray] = []

        if full_window_count:
            full_medians = np.median(
                coverage_array[:full_value_count].reshape(full_window_count, window_size),
                axis=1,
            )
            full_positions = (
                first_position
                + np.arange(full_window_count, dtype=np.int64) * window_size
                + (window_size - 1) // 2
            )
            median_parts.append(full_medians)
            position_parts.append(full_positions)

        if full_value_count < position_count:
            tail_median = np.array(
                [float(np.median(coverage_array[full_value_count:]))], dtype=float
            )
            tail_start = int(positions[full_value_count])
            tail_position = np.array(
                [tail_start + (last_position - tail_start) // 2], dtype=np.int64
            )
            median_parts.append(tail_median)
            position_parts.append(tail_position)

        return np.concatenate(position_parts), np.concatenate(median_parts)

    position_array = np.asarray(positions, dtype=np.int64)
    bin_ids = (position_array - first_position) // window_size
    boundaries = np.flatnonzero(np.diff(bin_ids)) + 1
    starts = np.concatenate((np.array([0]), boundaries))
    stops = np.concatenate((boundaries, np.array([position_count])))
    window_positions = np.empty(len(starts), dtype=np.int64)
    window_coverages = np.empty(len(starts), dtype=float)
    for output_index, (start_index, stop_index) in enumerate(zip(starts, stops)):
        window_positions[output_index] = (
            int(position_array[start_index]) + int(position_array[stop_index - 1])
        ) // 2
        window_coverages[output_index] = float(
            np.median(coverage_array[start_index:stop_index])
        )
    return window_positions, window_coverages


def is_cell_block_limit_error(error: BaseException) -> bool:
    """Return whether an exception chain contains Matplotlib's cell-block overflow."""
    current: BaseException | None = error
    visited: set[int] = set()
    while current is not None and id(current) not in visited:
        visited.add(id(current))
        if isinstance(current, OverflowError) and CELL_BLOCK_LIMIT_MESSAGE in str(current):
            return True
        current = current.__cause__ or current.__context__
    return False


def close_failed_plot_figures() -> None:
    """Release figures left open when Matplotlib raises during savefig()."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plt.close("all")


def render_coverage_plot_adaptive(
    positions: Sequence[int],
    coverages: Sequence[float],
    **plot_options: Any,
) -> int:
    """Render coverage, retrying cell-block overflows with wider median windows.

    Returns 1 for an unchanged per-base plot, or the successful window size in
    base pairs when downsampling was necessary.
    """
    try:
        render_coverage_plot(positions, coverages, **plot_options)
        return 1
    except Exception as exc:
        if not is_cell_block_limit_error(exc):
            raise
        close_failed_plot_figures()

    window_size = FIRST_COVERAGE_WINDOW_BP
    while True:
        window_positions, window_coverages = median_coverage_windows(
            positions,
            coverages,
            window_size,
        )
        eprint(
            f"Coverage plot exceeded Matplotlib's cell block limit; retrying with "
            f"median coverage in {window_size}-bp windows "
            f"({len(coverages):,} -> {len(window_coverages):,} points)"
        )
        try:
            render_coverage_plot(
                window_positions,
                window_coverages,
                **plot_options,
            )
            return window_size
        except Exception as exc:
            if not is_cell_block_limit_error(exc):
                raise
            close_failed_plot_figures()
            if len(window_coverages) <= 1:
                raise
            window_size *= 10


def resolve_nn_model(requested_path: Path | None) -> Path:
    """Find the bundled CNN checkpoint or validate an explicit path."""
    if requested_path is not None:
        path = requested_path.expanduser().resolve()
        if not path.is_file():
            raise ToolError(f"CNN model checkpoint does not exist: {path}")
        return path

    for path in bundled_asset_candidates("models", NN_MODEL_FILENAME):
        if path.is_file():
            return path
    raise ToolError(
        f"Bundled models/{NN_MODEL_FILENAME} was not found. Supply it explicitly "
        "with --nn-model."
    )


def load_nn_model(model_path: Path, threads: int) -> tuple[Any, Any, Any, Any]:
    """Load the embedded 600x200 three-class CNN architecture and checkpoint."""
    try:
        import torch
        import torch.nn as nn
        from PIL import Image
    except ImportError as exc:
        raise ToolError(
            "NN validation requires PyTorch and Pillow. Install the dependencies "
            "from environment.yml."
        ) from exc

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type == "cpu":
        # PyTorch convolutions use intra-op threads. Configure this before any
        # model work; the surrounding image pipeline is finished before CNN
        # inference starts, so it does not compete with the worker processes.
        torch.set_num_threads(threads)
        try:
            torch.set_num_interop_threads(max(1, min(threads, 4)))
        except RuntimeError:
            # PyTorch permits setting inter-op threads only once per process.
            # This can happen when main() is invoked repeatedly by a caller.
            pass

    class CNVDeepCNN(nn.Module):
        def __init__(self, num_classes: int = 3):
            super().__init__()
            self.body = nn.Sequential(
                nn.Conv2d(3, 32, kernel_size=3, padding=1),
                nn.ReLU(),
                nn.Conv2d(32, 32, kernel_size=3, padding=1),
                nn.ReLU(),
                nn.MaxPool2d(2, 2),
                nn.Conv2d(32, 64, kernel_size=3, padding=1),
                nn.ReLU(),
                nn.Conv2d(64, 64, kernel_size=3, padding=1),
                nn.ReLU(),
                nn.MaxPool2d(2, 2),
                nn.Conv2d(64, 128, kernel_size=3, padding=1),
                nn.ReLU(),
            )
            self.adaptive_pool = nn.AdaptiveAvgPool2d((4, 12))
            self.head = nn.Sequential(
                nn.Flatten(),
                nn.Linear(128 * 4 * 12, 64),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(64, num_classes),
            )

        def forward(self, tensor):
            tensor = self.body(tensor)
            tensor = self.adaptive_pool(tensor)
            return self.head(tensor)

    try:
        try:
            checkpoint = torch.load(model_path, map_location=device, weights_only=True)
        except TypeError:  # PyTorch versions predating weights_only.
            checkpoint = torch.load(model_path, map_location=device)
        state_dict = checkpoint.get("state_dict", checkpoint)
        model = CNVDeepCNN(num_classes=len(NN_CLASS_NAMES))
        model.load_state_dict(state_dict)
        model.to(device)
        model.eval()
    except Exception as exc:
        raise ToolError(f"Cannot load CNN checkpoint {model_path}: {exc}") from exc
    return model, torch, Image, device


def load_nn_image_array(
    image_path: Path,
    Image: Any,
) -> tuple[np.ndarray | None, str | None]:
    """Read and normalize one saved coverage image for batched inference."""
    try:
        with Image.open(image_path) as image:
            rgb_image = image.convert("RGB").resize(
                (600, 200),
                resample=Image.Resampling.BILINEAR,
            )
            image_array = np.asarray(rgb_image, dtype=np.float32).copy() / 255.0
        return image_array.transpose(2, 0, 1), None
    except Exception as exc:
        return None, f"Cannot read CNN image {image_path}: {exc}"


def predict_nn_images(
    model_bundle: tuple[Any, Any, Any, Any],
    image_paths: Sequence[Path],
    *,
    batch_size: int,
    loader_threads: int,
) -> tuple[dict[Path, tuple[str, np.ndarray]], dict[Path, str]]:
    """Preprocess images concurrently and classify them in PyTorch batches."""
    model, torch, Image, device = model_bundle
    predictions: dict[Path, tuple[str, np.ndarray]] = {}
    errors: dict[Path, str] = {}
    mean = torch.tensor([0.485, 0.456, 0.406], dtype=torch.float32).view(1, 3, 1, 1)
    std = torch.tensor([0.229, 0.224, 0.225], dtype=torch.float32).view(1, 3, 1, 1)

    loader: ThreadPoolExecutor | None = None
    if loader_threads > 1:
        loader = ThreadPoolExecutor(max_workers=loader_threads)
    try:
        for offset in range(0, len(image_paths), batch_size):
            batch_paths = list(image_paths[offset:offset + batch_size])
            if loader is None:
                loaded = [load_nn_image_array(path, Image) for path in batch_paths]
            else:
                futures = [loader.submit(load_nn_image_array, path, Image) for path in batch_paths]
                loaded = [future.result() for future in futures]

            valid_paths: list[Path] = []
            valid_arrays: list[np.ndarray] = []
            for image_path, (image_array, error) in zip(batch_paths, loaded):
                if error is not None or image_array is None:
                    errors[image_path] = error or f"Cannot read CNN image {image_path}"
                    continue
                valid_paths.append(image_path)
                valid_arrays.append(image_array)
            if not valid_arrays:
                continue

            tensor = torch.from_numpy(np.stack(valid_arrays, axis=0))
            tensor = ((tensor - mean) / std).to(device)
            try:
                with torch.inference_mode():
                    batch_probabilities = torch.softmax(model(tensor), dim=1)
                probability_rows = batch_probabilities.cpu().numpy()
            except Exception as exc:
                for image_path in valid_paths:
                    errors[image_path] = f"CNN inference failed: {exc}"
                continue

            for image_path, probability_array in zip(valid_paths, probability_rows):
                predicted_class = NN_CLASS_NAMES[int(probability_array.argmax())]
                predictions[image_path] = (predicted_class, probability_array)
    finally:
        if loader is not None:
            loader.shutdown(wait=True)

    return predictions, errors


def create_nn_plot(task: NNPlotTask) -> NNPlotResult:
    """Create one NN input image in an isolated worker process."""
    try:
        positions, coverages = samtools_depth(
            task.alignment_path,
            task.chrom,
            task.view_start,
            task.view_stop,
            task.samtools_threads,
        )
        render_coverage_plot_adaptive(
            positions,
            coverages,
            chrom=task.chrom,
            view_start=task.view_start,
            view_stop=task.view_stop,
            cnv_start=task.cnv_start,
            cnv_stop=task.cnv_stop,
            output_path=task.output_path,
            coverage_cap_multiplier=task.coverage_cap_multiplier,
            dpi=task.dpi,
            legend=False,
            cnv_borders=True,
        )
        return NNPlotResult(task.row_index, task.output_path)
    except Exception as exc:
        return NNPlotResult(task.row_index, task.output_path, str(exc))


def resolve_nn_alignment_paths(
    frame: pd.DataFrame,
    args: argparse.Namespace,
    single_sample: bool,
) -> dict[str, Path]:
    if not shutil.which("samtools"):
        raise ToolError("samtools is not available in PATH")

    if single_sample:
        path = args.alignment_path.expanduser().resolve()
        validate_alignment_path(path)
        return {"single_sample": path}

    alignment_dir = (
        args.alignment_dir.expanduser().resolve()
        if args.alignment_dir is not None
        else None
    )
    paths: dict[str, Path] = {}
    for run in frame["Run"].drop_duplicates().astype(str):
        rows = frame.loc[frame["Run"].astype(str) == run]
        if "_alignment_path" in rows:
            manifest_paths = rows["_alignment_path"].drop_duplicates().astype(str).tolist()
            if len(manifest_paths) != 1:
                eprint(f"[{run}] skipping alignment-dependent processing: multiple paths")
                continue
            path = Path(manifest_paths[0])
        else:
            assert alignment_dir is not None
            path = make_masked_path(alignment_dir, args.alignment_mask, run)
        try:
            validate_alignment_path(path)
        except ToolError as exc:
            eprint(f"[{run}] skipping alignment-dependent processing: {exc}")
            continue
        paths[run] = path
    return paths


def run_nn_validation(
    frame: pd.DataFrame,
    args: argparse.Namespace,
    *,
    single_sample: bool,
    inference_root: Path,
) -> pd.DataFrame:
    """Create coverage images in parallel and classify them in CNN batches."""
    vprint("Resolving and validating BAM/CRAM files for CNN validation")
    alignment_paths = resolve_nn_alignment_paths(frame, args, single_sample)
    output_root = inference_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    nn_result = pd.DataFrame(index=frame.index)
    nn_result["NN_predicted_SVTYPE"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    nn_result["NN_is_val"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    nn_result["NN_status"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    for class_name in NN_CLASS_NAMES:
        nn_result[f"NN_probability_{class_name}"] = np.nan

    tasks: list[NNPlotTask] = []
    task_by_row: dict[int, NNPlotTask] = {}
    eligible_count = int(frame["Run"].astype(str).isin(alignment_paths).sum())
    worker_count = min(args.threads, max(1, eligible_count))
    # samtools -@ counts additional decompression threads.
    samtools_threads = args.threads if worker_count == 1 else 0

    for run, rows in frame.groupby("Run", sort=False):
        run = str(run)
        if run not in alignment_paths:
            nn_result.loc[rows.index, "NN_status"] = "missing_alignment"
            continue
        alignment_path = alignment_paths[run]
        for row_number, (row_index, row) in enumerate(rows.iterrows(), start=1):
            chrom = str(row["chrom"])
            cnv_start = int(row["start"])
            cnv_stop = int(row["stop"])
            view_start, view_stop = coverage_plot_window(
                cnv_start,
                cnv_stop,
                int(row["length"]),
            )
            output_dir = plot_output_dir(output_root, run, row["SVTYPE"])
            output_path = output_dir / (
                f"{safe_path_token(run)}_{safe_path_token(row['ID'])}_{row_number}_"
                f"{safe_path_token(chrom)}_{cnv_start}_{cnv_stop}.jpeg"
            )
            task = NNPlotTask(
                row_index=int(row_index),
                run=run,
                cnv_id=str(row["ID"]),
                alignment_path=alignment_path,
                chrom=chrom,
                view_start=view_start,
                view_stop=view_stop,
                cnv_start=cnv_start,
                cnv_stop=cnv_stop,
                output_path=output_path,
                coverage_cap_multiplier=args.coverage_cap_multiplier,
                dpi=args.plot_dpi,
                samtools_threads=samtools_threads,
            )
            tasks.append(task)
            task_by_row[task.row_index] = task

    if not tasks:
        eprint("NN validation skipped: no valid alignment files were found")
        return nn_result

    vprint(
        f"NN image preparation: {len(tasks)} CNV(s), {worker_count} worker(s)"
    )
    plot_results: list[NNPlotResult] = []
    if worker_count == 1:
        for completed_count, task in enumerate(tasks, start=1):
            vprint(
                f"[{task.run}:{task.cnv_id}] preparing CNN coverage image "
                f"({completed_count}/{len(tasks)})"
            )
            plot_results.append(create_nn_plot(task))
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            future_to_task = {
                executor.submit(create_nn_plot, task): task for task in tasks
            }
            for completed_count, future in enumerate(
                as_completed(future_to_task),
                start=1,
            ):
                task = future_to_task[future]
                try:
                    plot_results.append(future.result())
                except Exception as exc:
                    plot_results.append(
                        NNPlotResult(task.row_index, task.output_path, str(exc))
                    )
                vprint(
                    f"[{task.run}:{task.cnv_id}] CNN image task finished "
                    f"({completed_count}/{len(tasks)})"
                )

    successful_paths: list[Path] = []
    row_by_path: dict[Path, int] = {}
    for plot_result in plot_results:
        task = task_by_row[plot_result.row_index]
        if plot_result.error is not None:
            nn_result.at[task.row_index, "NN_status"] = "nn_error"
            eprint(
                f"[{task.run}:{task.cnv_id}] NN image preparation error: "
                f"{plot_result.error}"
            )
            continue
        successful_paths.append(plot_result.output_path)
        row_by_path[plot_result.output_path] = task.row_index

    if not successful_paths:
        eprint("NN inference skipped: no coverage image was created successfully")
        return nn_result

    vprint(
        f"CNN image preparation complete: {len(successful_paths)}/{len(tasks)} "
        "image(s) ready"
    )

    # Import and initialize PyTorch only after process-based image preparation.
    # Forking a process after a threaded PyTorch runtime has been initialized can
    # deadlock on some platforms.
    model_path = resolve_nn_model(args.nn_model)
    vprint(f"Loading CNN model from {model_path}")
    model_bundle = load_nn_model(model_path, args.threads)
    device = model_bundle[3]
    thread_details = (
        f", PyTorch intra-op threads={args.threads}" if device.type == "cpu" else ""
    )
    vprint(f"NN inference: device={device}{thread_details}, batch size={args.nn_batch_size}")
    predictions, inference_errors = predict_nn_images(
        model_bundle,
        successful_paths,
        batch_size=args.nn_batch_size,
        loader_threads=min(args.threads, args.nn_batch_size),
    )
    for output_path, error in inference_errors.items():
        row_index = row_by_path[output_path]
        task = task_by_row[row_index]
        nn_result.at[row_index, "NN_status"] = "nn_error"
        eprint(f"[{task.run}:{task.cnv_id}] NN validation error: {error}")

    for output_path, (prediction, probabilities) in predictions.items():
        row_index = row_by_path[output_path]
        row = frame.loc[row_index]
        nn_result.at[row_index, "NN_predicted_SVTYPE"] = prediction
        nn_result.at[row_index, "NN_is_val"] = validation_code(
            str(row["SVTYPE"]),
            prediction,
        )
        nn_result.at[row_index, "NN_status"] = "ok"
        for class_name, probability in zip(NN_CLASS_NAMES, probabilities):
            nn_result.at[row_index, f"NN_probability_{class_name}"] = float(probability)

    vprint(f"CNN inference complete: {len(predictions)} prediction(s)")

    return nn_result


def resolve_sample_output_dir(output_template: Path, sample_id: str) -> Path:
    """Resolve a per-sample output directory, honoring {sample_id} templates."""
    raw_template = str(output_template.expanduser())
    safe_sample = safe_path_token(sample_id)
    if "{sample_id}" in raw_template:
        try:
            rendered = raw_template.format(sample_id=safe_sample)
        except (KeyError, IndexError, ValueError) as exc:
            raise ToolError(f"Invalid --output mask {raw_template!r}: {exc}") from exc
    else:
        if "{" in raw_template or "}" in raw_template:
            raise ToolError(
                "Unsupported placeholder in --output; only {sample_id} is accepted"
            )
        rendered = str(Path(raw_template) / safe_sample)
    return Path(rendered).resolve()


def run_validation_coverage_profiles(
    frame: pd.DataFrame,
    args: argparse.Namespace,
    *,
    single_sample: bool,
) -> pd.Series:
    """Create user-facing full coverage profiles and return their absolute paths."""
    profile_paths = pd.Series(pd.NA, index=frame.index, dtype="string")
    vprint("Resolving and validating BAM/CRAM files for coverage profiles")
    alignment_paths = resolve_nn_alignment_paths(frame, args, single_sample)
    total_profiles = int(frame["Run"].astype(str).isin(alignment_paths).sum())
    completed_profiles = 0
    vprint(
        f"Coverage profile preparation: {total_profiles} CNV(s), "
        f"{len(alignment_paths)} alignment file(s)"
    )

    for run, rows in frame.groupby("Run", sort=False):
        run = str(run)
        if run not in alignment_paths:
            continue
        alignment_path = alignment_paths[run]
        sample_output_dir = resolve_sample_output_dir(args.output, run)
        vprint(f"[{run}] creating {len(rows)} coverage profile(s) from {alignment_path}")
        for row_number, (row_index, row) in enumerate(rows.iterrows(), start=1):
            chrom = str(row["chrom"])
            cnv_start = int(row["start"])
            cnv_stop = int(row["stop"])
            view_start, view_stop = coverage_plot_window(
                cnv_start,
                cnv_stop,
                int(row["length"]),
            )
            output_dir = (
                sample_output_dir
                / "coverage_profiles"
                / plot_svtype_folder_name(row["SVTYPE"])
            )
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / (
                f"{safe_path_token(run)}_{safe_path_token(row['ID'])}_{row_number}_"
                f"{safe_path_token(chrom)}_{cnv_start}_{cnv_stop}.{args.plot_format}"
            )
            try:
                vprint(
                    f"[{run}:{row['ID']}] reading coverage for "
                    f"{chrom}:{view_start}-{view_stop} "
                    f"({completed_profiles + 1}/{total_profiles})"
                )
                positions, coverages = samtools_depth(
                    alignment_path,
                    chrom,
                    view_start,
                    view_stop,
                    args.threads,
                )
                render_coverage_plot_adaptive(
                    positions,
                    coverages,
                    chrom=chrom,
                    view_start=view_start,
                    view_stop=view_stop,
                    cnv_start=cnv_start,
                    cnv_stop=cnv_stop,
                    output_path=output_path,
                    coverage_cap_multiplier=args.coverage_cap_multiplier,
                    dpi=args.plot_dpi,
                    legend=args.legend,
                    cnv_borders=args.cnv_borders,
                )
            except Exception as exc:
                eprint(f"[{run}:{row['ID']}] coverage profile error: {exc}")
                completed_profiles += 1
                continue
            profile_paths.at[row_index] = str(output_path.resolve())
            completed_profiles += 1
            vprint(f"[{run}:{row['ID']}] wrote coverage profile {output_path}")

    vprint(
        f"Coverage profile stage complete: {profile_paths.notna().sum()}/"
        f"{total_profiles} profile(s) created"
    )

    return profile_paths


def empty_nn_result(frame: pd.DataFrame, status: str = "nn_error") -> pd.DataFrame:
    """Return a complete row-aligned CNN result after a global CNN failure."""
    result = pd.DataFrame(index=frame.index)
    result["NN_predicted_SVTYPE"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    result["NN_is_val"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    result["NN_status"] = status
    for class_name in NN_CLASS_NAMES:
        result[f"NN_probability_{class_name}"] = np.nan
    return result


def empty_vaf_result(frame: pd.DataFrame, status: str = "model_error") -> pd.DataFrame:
    """Return a complete row-aligned VAF result after a global VAF failure."""
    result = frame.copy()
    for column in FEATURE_COLUMNS:
        result[column] = np.nan
    result["N_VAF"] = 0
    result["predicted_SVTYPE"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    result["is_val"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    result["status"] = status
    for class_name in VAF_MODEL_CLASSES:
        result[f"model_score_{class_name}"] = np.nan
    return result


def run_vaf_validation(frame: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Extract VAF features and perform SVC inference without aborting other modes."""
    try:
        model_path = resolve_vaf_model()
        vprint(f"Loading VAF SVC model from {model_path}")
        model = load_vaf_model(model_path)
    except ToolError as exc:
        eprint(f"VAF model error: {exc}")
        return empty_vaf_result(frame)

    values: dict[int, list[float]] = {int(index): [] for index in frame.index}
    sample_state: dict[str, tuple[str, str]] = {}
    for run, rows in frame.groupby("Run", sort=False):
        vprint(f"[{run}] processing {len(rows)} CNV interval(s) for VAF")
        state, message = process_sample(str(run), rows, values, args)
        sample_state[str(run)] = (state, message)
        if state == "parsed":
            vprint(f"[{run}] {state}: {message}")
        else:
            eprint(f"[{run}] {state}: {message}")

    result = empty_vaf_result(frame, status="unprocessed")
    predictable_indices: list[int] = []
    predictable_features: list[list[float]] = []
    for row_index, row in result.iterrows():
        run = str(row["Run"])
        sample_status, _ = sample_state[run]
        region_values = values[int(row_index)]
        result.at[row_index, "N_VAF"] = len(region_values)

        if sample_status != "parsed":
            result.at[row_index, "status"] = sample_status
            continue
        if len(region_values) < args.min_variants:
            result.at[row_index, "status"] = "insufficient_vaf"
            continue
        try:
            features = calculate_features(region_values)
        except ToolError as exc:
            result.at[row_index, "status"] = "kde_error"
            eprint(f"[{run}:{row['ID']}] KDE error: {exc}")
            continue

        for column, value in features.items():
            result.at[row_index, column] = value
        result.at[row_index, "status"] = (
            "ok" if len(region_values) >= 50 else "ok_low_vaf_count"
        )
        predictable_indices.append(int(row_index))
        predictable_features.append([features[column] for column in MODEL_FEATURE_COLUMNS])

    if not predictable_indices:
        vprint("VAF feature extraction complete: no CNV rows are eligible for inference")
        return result

    matrix = np.asarray(predictable_features, dtype=float)
    vprint(
        f"VAF feature extraction complete: {len(predictable_indices)}/{len(frame)} "
        "CNV row(s) eligible; running SVC inference"
    )
    try:
        predictions, probabilities, probability_classes = model_predict(
            model,
            matrix,
        )
    except ToolError as exc:
        eprint(f"VAF model inference error: {exc}")
        result.loc[predictable_indices, "status"] = "model_error"
        return result

    probability_offsets = {
        class_name: offset
        for offset, class_name in enumerate(probability_classes)
    }
    for offset, row_index in enumerate(predictable_indices):
        prediction = str(predictions[offset])
        result.at[row_index, "predicted_SVTYPE"] = output_label(prediction)
        result.at[row_index, "is_val"] = validation_code(
            str(result.at[row_index, "SVTYPE"]),
            prediction,
        )
        for class_name in VAF_MODEL_CLASSES:
            result.at[row_index, f"model_score_{class_name}"] = probabilities[
                offset,
                probability_offsets[class_name],
            ]
    vprint(f"VAF inference complete: {len(predictable_indices)} prediction(s)")
    return result


def output_base_frame(frame: pd.DataFrame, single_sample: bool) -> pd.DataFrame:
    """Create the common public columns with their final uppercase names."""
    result = pd.DataFrame(index=frame.index)
    if not single_sample:
        result["SAMPLE_ID"] = frame["Run"].astype("string")
    result["CHR"] = frame["chrom"].astype("string")
    result["START"] = frame["start"].astype(np.int64)
    result["END"] = frame["stop"].astype(np.int64)
    result["INPUT_SVTYPE"] = frame["SVTYPE"].astype("string")
    return result


def build_vaf_output(
    frame: pd.DataFrame,
    vaf_result: pd.DataFrame,
    annotations: pd.Series,
    coverage_paths: pd.Series,
    train_dimension_paths: pd.Series,
    *,
    single_sample: bool,
) -> pd.DataFrame:
    result = output_base_frame(frame, single_sample)
    statuses = vaf_result["status"].astype("string")
    result["PREDICTED_SVTYPE"] = vaf_result["predicted_SVTYPE"]
    result["INPUT_OUTPUT_CONSISTENCY"] = vaf_result["is_val"]
    result["MODEL_ELIGIBILITY"] = statuses.map(
        lambda value: (
            "ELIGIBLE"
            if value in {"ok", "ok_low_vaf_count"}
            else "INSUFFICIENT_NUM_VAR"
            if value == "insufficient_vaf"
            else "ERROR"
        )
    )
    result["LOW_MAPPABILITY_REGION"] = annotations
    result["INPUT_REGION_COVERAGE"] = coverage_paths
    result["NUM_OF_VARIANTS"] = vaf_result["N_VAF"].astype(np.int64)
    for source, target in (
        ("X0.33", "MARKER_0.33"),
        ("X0.66", "MARKER_0.66"),
        ("X0.33_X0.66", "MARKER_0.33+0.66"),
        ("X0.5", "MARKER_0.5"),
        ("X1", "MARKER_1"),
        ("model_score_DEL", "MODEL_SCORE_DEL"),
        ("model_score_DIP", "MODEL_SCORE_DIP"),
        ("model_score_DUP", "MODEL_SCORE_DUP"),
    ):
        result[target] = vaf_result[source]
    result["CNV_POSITION_IN_TRAIN_DIMENSION"] = train_dimension_paths
    return result


def build_cnn_output(
    frame: pd.DataFrame,
    nn_result: pd.DataFrame,
    annotations: pd.Series,
    coverage_paths: pd.Series,
    *,
    single_sample: bool,
) -> pd.DataFrame:
    result = output_base_frame(frame, single_sample)
    result["PREDICTED_SVTYPE"] = nn_result["NN_predicted_SVTYPE"]
    result["INPUT_OUTPUT_CONSISTENCY"] = nn_result["NN_is_val"]
    result["MODEL_ELIGIBILITY"] = nn_result["NN_status"].map(
        lambda value: "ELIGIBLE" if value == "ok" else "ERROR"
    )
    result["LOW_MAPPABILITY_REGION"] = annotations
    result["INPUT_REGION_COVERAGE"] = coverage_paths
    for class_name in NN_CLASS_NAMES:
        result[f"NN_PROBABILITY_{class_name}"] = nn_result[
            f"NN_probability_{class_name}"
        ]
    return result


def combined_final_decision(
    input_svtype: str,
    vaf_prediction: object,
    nn_prediction: object,
) -> object:
    predictions = [
        str(value)
        for value in (vaf_prediction, nn_prediction)
        if not pd.isna(value)
    ]
    if not predictions:
        return pd.NA
    if all(prediction == input_svtype for prediction in predictions):
        return input_svtype
    return "CNV_not_input type"


def build_combined_output(
    frame: pd.DataFrame,
    vaf_result: pd.DataFrame,
    nn_result: pd.DataFrame,
    *,
    single_sample: bool,
) -> pd.DataFrame:
    result = output_base_frame(frame, single_sample)
    result["VAF_PREDICTED_SVTYPE"] = vaf_result["predicted_SVTYPE"]
    result["NN_PREDICTED_SVTYPE"] = nn_result["NN_predicted_SVTYPE"]
    result["FINAL_DECISION"] = [
        combined_final_decision(input_svtype, vaf_prediction, nn_prediction)
        for input_svtype, vaf_prediction, nn_prediction in zip(
            result["INPUT_SVTYPE"],
            result["VAF_PREDICTED_SVTYPE"],
            result["NN_PREDICTED_SVTYPE"],
        )
    ]
    return result


def write_sample_outputs(
    frame: pd.DataFrame,
    output_template: Path,
    tables: dict[str, pd.DataFrame],
    *,
    precision: int,
) -> None:
    float_format = f"%.{precision}f"
    for run, rows in frame.groupby("Run", sort=False):
        sample_id = str(run)
        sample_output_dir = resolve_sample_output_dir(output_template, sample_id)
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        for filename, table in tables.items():
            output_path = sample_output_dir / filename
            table.loc[rows.index].to_csv(
                output_path,
                sep="\t",
                index=False,
                na_rep="NA",
                float_format=float_format,
                quoting=csv.QUOTE_MINIMAL,
            )
            vprint(f"[{sample_id}] wrote {len(rows)} row(s) to {output_path}")


def run_plot_cov_mode(frame: pd.DataFrame, args: argparse.Namespace) -> int:
    if not shutil.which("samtools"):
        raise ToolError("samtools is not available in PATH")

    alignment_dir = (
        args.alignment_dir.expanduser().resolve()
        if args.alignment_dir is not None
        else None
    )
    output_root = args.plot_dir.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    vprint(f"Coverage-only output directory: {output_root}")
    vprint("Resolving and validating BAM/CRAM files for coverage-only mode")

    sample_paths: dict[str, Path] = {}
    for run in frame["Run"].drop_duplicates().astype(str):
        rows = frame.loc[frame["Run"].astype(str) == run]
        if "_alignment_path" in rows:
            manifest_paths = rows["_alignment_path"].drop_duplicates().astype(str).tolist()
            if len(manifest_paths) != 1:
                eprint(f"[{run}] skipping coverage plots: multiple alignment paths")
                continue
            path = Path(manifest_paths[0])
        else:
            assert alignment_dir is not None
            path = make_masked_path(alignment_dir, args.alignment_mask, run)
        try:
            validate_alignment_path(path)
        except ToolError as exc:
            eprint(f"[{run}] skipping coverage plots: {exc}")
            continue
        sample_paths[run] = path
        vprint(f"[{run}] using alignment file {path}")

    plot_count = 0
    for run, rows in frame.groupby("Run", sort=False):
        run = str(run)
        if run not in sample_paths:
            continue
        alignment_path = sample_paths[run]
        for row_number, (_, row) in enumerate(rows.iterrows(), start=1):
            chrom = str(row["chrom"])
            cnv_start = int(row["start"])
            cnv_stop = int(row["stop"])
            cnv_length = int(row["length"])
            view_start, view_stop = coverage_plot_window(
                cnv_start, cnv_stop, cnv_length
            )
            output_dir = plot_output_dir(output_root, run, row["SVTYPE"])
            filename = (
                f"{safe_path_token(run)}_{safe_path_token(row['ID'])}_{row_number}_"
                f"{safe_path_token(chrom)}_{cnv_start}_{cnv_stop}.{args.plot_format}"
            )
            output_path = output_dir / filename
            vprint(
                f"[{run}:{row['ID']}] plotting {chrom}:{view_start}-{view_stop} "
                f"from {alignment_path}"
            )
            positions, coverages = samtools_depth(
                alignment_path,
                chrom,
                view_start,
                view_stop,
                args.threads,
            )
            render_coverage_plot_adaptive(
                positions,
                coverages,
                chrom=chrom,
                view_start=view_start,
                view_stop=view_stop,
                cnv_start=cnv_start,
                cnv_stop=cnv_stop,
                output_path=output_path,
                coverage_cap_multiplier=args.coverage_cap_multiplier,
                dpi=args.plot_dpi,
                legend=args.legend,
                cnv_borders=args.cnv_borders,
            )
            plot_count += 1

    vprint(f"Created {plot_count} coverage plot(s) in {output_root}")
    return 0


def run_single_plot_cov_mode(args: argparse.Namespace) -> int:
    """Render exactly one coverage image for an explicit genomic interval."""
    if not shutil.which("samtools"):
        raise ToolError("samtools is not available in PATH")

    alignment_path = args.alignment_path.expanduser().resolve()
    validate_alignment_path(alignment_path)

    chrom = str(args.single_chrom).strip()
    cnv_start = int(args.single_start)
    cnv_stop = int(args.single_end)
    if not chrom:
        raise ToolError("--chr cannot be empty")
    if cnv_start < 1 or cnv_stop < cnv_start:
        raise ToolError("Single-interval coordinates must satisfy 1 <= start <= end")

    cnv_length = max(1, cnv_stop - cnv_start)
    view_start, view_stop = coverage_plot_window(
        cnv_start,
        cnv_stop,
        cnv_length,
    )
    output_root = args.plot_dir.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    filename = (
        f"{safe_path_token(alignment_path.stem)}_{safe_path_token(chrom)}_"
        f"{cnv_start}_{cnv_stop}.{args.plot_format}"
    )
    output_path = output_root / filename

    vprint(
        f"Plotting {chrom}:{view_start}-{view_stop} from {alignment_path}"
    )
    positions, coverages = samtools_depth(
        alignment_path,
        chrom,
        view_start,
        view_stop,
        args.threads,
    )
    render_coverage_plot_adaptive(
        positions,
        coverages,
        chrom=chrom,
        view_start=view_start,
        view_stop=view_stop,
        cnv_start=cnv_start,
        cnv_stop=cnv_stop,
        output_path=output_path,
        coverage_cap_multiplier=args.coverage_cap_multiplier,
        dpi=args.plot_dpi,
        legend=args.legend,
        cnv_borders=args.cnv_borders,
    )
    vprint(f"Created coverage plot: {output_path}")
    return 0

def resolve_annotation_bed(requested_path: Path | None) -> Path:
    if requested_path is not None:
        path = requested_path.expanduser().resolve()
        if not path.is_file():
            raise ToolError(f"Centromere/telomere BED file does not exist: {path}")
        return path

    filename = "hg38_centromeres_telomeres.bed"
    for path in bundled_asset_candidates(filename):
        if path.exists():
            return path
    raise ToolError(
        f"Bundled {filename} was not found. Supply it explicitly with "
        "--centromere-telomere-bed."
    )


def read_annotation_bed(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    features: dict[str, list[tuple[int, int, str]]] = {}
    try:
        handle = path.open("rt", encoding="utf-8")
    except OSError as exc:
        raise ToolError(f"Cannot open annotation BED {path}: {exc}") from exc

    with handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ToolError(f"Annotation BED line {line_number} has fewer than 4 columns")
            try:
                start = int(fields[1])
                stop = int(fields[2])
            except ValueError as exc:
                raise ToolError(
                    f"Annotation BED line {line_number} has invalid coordinates"
                ) from exc
            if start < 0 or stop <= start:
                raise ToolError(f"Annotation BED line {line_number} has an invalid interval")
            feature_name = fields[3].strip().lower()
            if feature_name == "centromere":
                output_name = "CENTROMERIC"
            elif feature_name == "telomere":
                output_name = "TELOMERIC"
            else:
                raise ToolError(
                    f"Unsupported annotation {fields[3]!r} on BED line {line_number}; "
                    "expected centromere or telomere"
                )
            features.setdefault(canonical_chrom(fields[0]), []).append(
                (start, stop, output_name)
            )
    if not features:
        raise ToolError(f"Annotation BED contains no centromere/telomere intervals: {path}")
    return features


def annotate_special_regions(
    frame: pd.DataFrame,
    features: dict[str, list[tuple[int, int, str]]],
) -> pd.Series:
    annotations = pd.Series(pd.NA, index=frame.index, dtype="string")
    label_order = ("CENTROMERIC", "TELOMERIC")
    for row_index, row in frame.iterrows():
        # CNV coordinates are 1-based inclusive; BED intervals are 0-based half-open.
        cnv_start = int(row["start"]) - 1
        cnv_stop = int(row["stop"])
        overlaps = {
            label
            for feature_start, feature_stop, label in features.get(
                canonical_chrom(row["chrom"]), []
            )
            if cnv_start < feature_stop and feature_start < cnv_stop
        }
        if overlaps:
            annotations.at[row_index] = ";".join(
                label for label in label_order if label in overlaps
            )
    return annotations

def make_single_cnv_input(
    chrom: str,
    start: int,
    stop: int,
    svtype: str,
) -> TableInput:
    normalized_chrom = chrom.strip()
    if not normalized_chrom:
        raise ToolError("--chr cannot be empty")
    normalized_type = svtype.upper().strip()
    if normalized_type not in VALIDATION_SVTYPES:
        raise ToolError("--svtype must be DEL, DIP, or DUP")
    if start < 1 or stop < start:
        raise ToolError("Single-CNV coordinates must satisfy 1 <= start <= end")
    frame = pd.DataFrame(
        [
            {
                "Run": "single_sample",
                "chrom": normalized_chrom,
                "start": start,
                "stop": stop,
                "length": max(1, stop - start),
                "SVTYPE": normalized_type,
                "ID": "cnv_1",
            }
        ]
    )
    return TableInput(
        frame=frame,
        single_sample=True,
        base_columns=_output_base_columns(True),
    )


def resolve_example_dir() -> Path:
    for path in bundled_asset_candidates("example"):
        if path.is_dir():
            return path
    raise ToolError(
        "Bundled example directory was not found in the project or installation "
        "share directory"
    )


def ensure_test_alignment_index(alignment_path: Path, threads: int) -> None:
    """Create the bundled example CRAM index on first test-run when needed."""
    try:
        validate_alignment_path(alignment_path)
        return
    except ToolError as exc:
        if "Alignment index was not found" not in str(exc):
            raise
    if not shutil.which("samtools"):
        raise ToolError("samtools is required to index the example CRAM")
    command = ["samtools", "index", "-@", str(threads), str(alignment_path)]
    vprint(f"Creating missing test CRAM index for {alignment_path}")
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise ToolError(f"Cannot run samtools index for test-run: {exc}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise ToolError(f"Cannot index example CRAM: {detail}")
    validate_alignment_path(alignment_path)
    vprint(f"Created test CRAM index for {alignment_path}")


def run_test_mode(threads: int, verbose: bool) -> int:
    """Run the complete single-sample workflow on bundled example data."""
    example_dir = resolve_example_dir()
    cnvs = example_dir / "ClinCNV_output_chr22.tsv"
    vcf = example_dir / "HG00577_example.vcf.gz"
    alignment = example_dir / "HG00577_example.cram"
    missing = [path for path in (cnvs, vcf, alignment) if not path.is_file()]
    if missing:
        raise ToolError(
            "Test-run file(s) are missing: " + ", ".join(str(path) for path in missing)
        )
    ensure_test_alignment_index(alignment, threads)

    output_dir = Path.cwd() / "test_run"
    vprint(f"Running bundled example from {example_dir}")
    vprint(f"Test output directory: {output_dir}")
    return main(
        [
            "--cnvs", str(cnvs),
            "--vcf-path", str(vcf),
            "--output", str(output_dir),
            "--backend", "python",
            "--cnn-validation", "true",
            "--alignment-path", str(alignment),
            "--threads", str(threads),
            "--verbose", str(verbose).lower(),
        ]
    )


def parse_bool(value: str) -> bool:
    """Parse an explicit true/false command-line value."""
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "on"}:
        return True
    if normalized in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def parse_sep(value: str | None) -> str | None:
    if value is None or value.lower() == "auto":
        return None
    aliases = {"tab": "\t", "tsv": "\t", "comma": ",", "csv": ","}
    return aliases.get(value.lower(), value.encode("utf-8").decode("unicode_escape"))


class HelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate CNVs from VAF distributions or plot CNV coverage profiles.",
        formatter_class=HelpFormatter,
        epilog="""Standard input table layouts (the first row is a header; names are ignored):
  multi-sample:  sample_id  chr  start  end  SVTYPE
  single-sample: chr        start end    SVTYPE

Native ClinCNV tables and Manta VCF output are detected automatically.
Coordinates are 1-based inclusive. SVTYPE is DEL, DIP, or DUP.
See readme.md and help.md for the complete input/output contract.""",
    )
    parser.add_argument(
        "--test-run",
        action="store_true",
        help=(
            "Run the bundled ClinCNV, VCF, and CRAM example and write "
            "per-model TSV files plus full coverage profiles under test_run"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("validate", "plot-cov"),
        default="validate",
        help="Operation mode",
    )
    parser.add_argument(
        "--verbose",
        type=parse_bool,
        default=True,
        metavar="{true,false}",
        help=(
            "Print progress details for input parsing, coverage generation, "
            "model inference, and output writing"
        ),
    )
    parser.add_argument(
        "-c", "--cnvs", type=Path,
        help=(
            "Standard CNV TSV/CSV, multi-sample manifest, native ClinCNV output, "
            "or Manta VCF; omit for explicit single-CNV validation or explicit "
            "single-interval plot-cov"
        ),
    )
    parser.add_argument(
        "--input-sep",
        default="auto",
        help="Input separator: auto, tab, comma, or a literal separator",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help=(
            "Parallelism limit for NN image workers and PyTorch CPU inference; "
            "also passed to supported samtools/bcftools operations"
        ),
    )

    interval = parser.add_argument_group("Single genomic interval")
    interval.add_argument("--chr", dest="single_chrom", help="Chromosome")
    interval.add_argument("--start", dest="single_start", type=int, help="Interval start")
    interval.add_argument("--end", dest="single_end", type=int, help="Interval end")

    validation = parser.add_argument_group("Validation")
    validation.add_argument(
        "--vaf-validation",
        type=parse_bool,
        default=True,
        metavar="{true,false}",
        help="Run validation with the bundled VAF SVC model",
    )
    validation.add_argument(
        "--plot-vaf-train-dimension",
        type=parse_bool,
        default=True,
        metavar="{true,false}",
        help=(
            "Create one JPEG per plottable CNV in the bundled VAF training "
            "projection; applies only when VAF validation is enabled"
        ),
    )
    validation.add_argument(
        "-v", "--vcf-dir", type=Path,
        help="Directory containing VCF files for multi-sample VAF validation",
    )
    validation.add_argument(
        "-m",
        "--vcf-mask",
        default="{sample}/{sample}.output.vcf.gz",
        help=(
            "VCF path template relative to --vcf-dir, or an absolute template. "
            "The {sample} placeholder may occur more than once."
        ),
    )
    validation.add_argument(
        "-o", "--output", type=Path,
        help=(
            "Output directory or sample-aware mask such as "
            "'/path/to/{sample_id}/output_dir'; required in validation mode"
        ),
    )
    validation.add_argument(
        "--vcf-path", type=Path,
        help="One small-variant VCF used for single-sample VAF validation",
    )
    validation.add_argument(
        "--svtype", type=lambda value: value.upper(), choices=("DEL", "DIP", "DUP"),
        help="Original type of one explicit CNV",
    )
    validation.add_argument(
        "--backend",
        choices=("auto", "bcftools", "python"),
        default="auto",
        help="VCF extraction backend; auto prefers bcftools and falls back to Python",
    )
    validation.add_argument(
        "--qual-cutoff",
        type=float,
        default=0.29,
        help="Keep records with QUAL strictly greater than this value; use a negative value to disable",
    )
    validation.add_argument(
        "--filter",
        dest="filter_name",
        default="PASS",
        help="Required FILTER token; use ANY to disable FILTER filtering",
    )
    validation.add_argument(
        "--min-variants",
        type=int,
        default=50,
        help="Minimum VAF observations required to calculate KDE features and predict",
    )
    validation.add_argument(
        "--precision",
        type=int,
        default=6,
        help="Number of decimal places used for feature and score columns",
    )
    validation.add_argument(
        "--centromere-telomere-bed",
        type=Path,
        help=(
            "BED file with centromere/telomere intervals. By default the bundled "
            "hg38_centromeres_telomeres.bed is used"
        ),
    )
    validation.add_argument(
        "--cnn-validation",
        "--NN-validation",
        "--nn-validation",
        dest="nn_validation",
        type=parse_bool,
        default=False,
        metavar="{true,false}",
        help="Run validation with the bundled coverage-image CNN",
    )
    validation.add_argument(
        "--nn-model",
        type=Path,
        help=f"CNN checkpoint; defaults to the bundled models/{NN_MODEL_FILENAME}",
    )
    validation.add_argument(
        "--nn-batch-size",
        type=int,
        default=4,
        help="Number of coverage images classified by the CNN per inference batch",
    )

    plotting = parser.add_argument_group("Coverage plot mode")
    plotting.add_argument(
        "--plot-cov",
        type=parse_bool,
        default=True,
        metavar="{true,false}",
        help="Create full user-facing coverage profiles during validation",
    )
    plotting.add_argument(
        "--alignment-path",
        "--bam-path",
        "--cram-path",
        dest="alignment_path",
        type=Path,
        help=(
            "One BAM/CRAM file for explicit single-interval plot-cov or "
            "single-sample validation with CNN/profiles"
        ),
    )
    plotting.add_argument(
        "--alignment-dir",
        "--bam-dir",
        "--cram-dir",
        dest="alignment_dir",
        type=Path,
        help=(
            "Base directory containing BAM/CRAM files; required for table plot-cov "
            "and multi-sample validation with CNN/profiles"
        ),
    )
    plotting.add_argument(
        "--alignment-mask",
        "--bam-mask",
        "--cram-mask",
        dest="alignment_mask",
        default="{sample}.markdup.sorted.cram",
        help=(
            "BAM/CRAM path template relative to --alignment-dir, or an absolute template. "
            "Supports {sample} and {Run}."
        ),
    )
    plotting.add_argument(
        "--plot-dir",
        type=Path,
        default=Path("coverage_plots"),
        help="Output directory used by standalone --mode plot-cov",
    )
    plotting.add_argument(
        "--plot-format",
        "--image-format",
        choices=("jpeg", "svg"),
        default="jpeg",
        help="Coverage image format",
    )
    plotting.add_argument(
        "--coverage-cap-multiplier",
        type=float,
        default=5.0,
        help="Clip the y axis at median window coverage multiplied by this value",
    )
    plotting.add_argument(
        "--plot-dpi",
        type=int,
        default=200,
        help="Resolution for JPEG plots; SVG output is vector-based",
    )
    plotting.add_argument(
        "--legend", type=parse_bool, default=True, metavar="{true,false}",
        help="Show coordinate/coverage axes, tick labels, and x-grid",
    )
    plotting.add_argument(
        "--cnv-borders", type=parse_bool, default=True, metavar="{true,false}",
        help="Draw red vertical lines at the CNV boundaries",
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return parser

def main(argv: Sequence[str] | None = None) -> int:
    global VERBOSE_ENABLED

    parser = build_parser()
    args = parser.parse_args(argv)
    VERBOSE_ENABLED = args.verbose
    vprint(
        f"Starting vacca {VERSION}: mode={args.mode}, "
        f"threads={args.threads}"
    )

    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.nn_batch_size < 1:
        parser.error("--nn-batch-size must be at least 1")
    if args.test_run:
        try:
            return run_test_mode(args.threads, args.verbose)
        except ToolError as exc:
            parser.error(str(exc))

    interval_arguments = (
        args.single_chrom,
        args.single_start,
        args.single_end,
    )
    validation_coordinate_arguments = interval_arguments + (args.svtype,)
    if args.mode == "plot-cov":
        if not math.isfinite(args.coverage_cap_multiplier) or args.coverage_cap_multiplier <= 0:
            parser.error("--coverage-cap-multiplier must be a positive finite number")
        if args.plot_dpi < 1:
            parser.error("--plot-dpi must be at least 1")
        invalid_validation_options = {
            "--vcf-path": args.vcf_path,
            "--vcf-dir": args.vcf_dir,
            "--output": args.output,
            "--svtype": args.svtype,
            "--centromere-telomere-bed": args.centromere_telomere_bed,
            "--cnn-validation": True if args.nn_validation else None,
            "--nn-model": args.nn_model,
        }
        used_invalid_options = [
            name for name, value in invalid_validation_options.items() if value is not None
        ]
        if used_invalid_options:
            parser.error(
                "Validation-only option(s) cannot be used in plot-cov mode: "
                + ", ".join(used_invalid_options)
            )

        if args.cnvs is None:
            required_plot = {
                "--chr": args.single_chrom,
                "--start": args.single_start,
                "--end": args.single_end,
                "--alignment-path": args.alignment_path,
            }
            missing = [name for name, value in required_plot.items() if value is None]
            if missing:
                parser.error(
                    "Explicit single-interval plot-cov requires: " + ", ".join(missing)
                )
            if args.alignment_dir is not None:
                parser.error("Use --alignment-path, not --alignment-dir, for one interval")
            try:
                vprint("Starting explicit single-interval coverage mode")
                return run_single_plot_cov_mode(args)
            except ToolError as exc:
                parser.error(str(exc))

        if any(value is not None for value in interval_arguments):
            parser.error("--chr/--start/--end cannot be combined with --cnvs in plot-cov mode")
        if args.alignment_path is not None:
            parser.error("--alignment-path cannot be combined with --cnvs")
        if args.alignment_dir is None:
            parser.error("--alignment-dir is required for table plot-cov")
        try:
            vprint(f"Loading CNV intervals for coverage-only mode from {args.cnvs}")
            table_input = read_table(
                args.cnvs.resolve(),
                parse_sep(args.input_sep),
                single_vcf=False,
                allowed_svtypes=PLOT_COV_SVTYPES,
            )
            vprint(
                f"Starting coverage-only mode for {len(table_input.frame)} CNV row(s)"
            )
            return run_plot_cov_mode(table_input.frame, args)
        except ToolError as exc:
            parser.error(str(exc))

    if not args.vaf_validation and not args.nn_validation:
        parser.error(
            "Validation mode requires --vaf-validation true, "
            "--cnn-validation true, or both"
        )
    if args.output is None:
        parser.error("--output is required in validation mode")
    if args.vcf_path is not None and args.vcf_dir is not None:
        parser.error("--vcf-path and --vcf-dir are mutually exclusive")
    if args.cnvs is None:
        required_single = {
            "--chr": args.single_chrom,
            "--start": args.single_start,
            "--end": args.single_end,
            "--svtype": args.svtype,
        }
        missing = [name for name, value in required_single.items() if value is None]
        if missing:
            parser.error(
                "Explicit single-CNV validation requires: " + ", ".join(missing)
            )
        try:
            vprint("Building input for one explicit CNV")
            table_input = make_single_cnv_input(
                args.single_chrom,
                args.single_start,
                args.single_end,
                args.svtype,
            )
        except ToolError as exc:
            parser.error(str(exc))
    else:
        if any(value is not None for value in validation_coordinate_arguments):
            parser.error("--chr/--start/--end/--svtype cannot be combined with --cnvs")
        try:
            vprint(f"Loading validation input from {args.cnvs}")
            table_input = read_table(
                args.cnvs.resolve(),
                parse_sep(args.input_sep),
                single_vcf=args.vcf_path is not None,
            )
        except ToolError as exc:
            parser.error(str(exc))

    frame = table_input.frame
    vprint(
        f"Input ready: {len(frame)} CNV row(s), "
        f"{frame['Run'].nunique()} sample(s), "
        f"layout={'single-sample' if table_input.single_sample else 'multi-sample'}"
    )
    manifest_vcf = "_vcf_path" in frame
    manifest_alignment = "_alignment_path" in frame
    if args.vaf_validation:
        if manifest_vcf and (args.vcf_path is not None or args.vcf_dir is not None):
            parser.error(
                "The multi-sample manifest already defines VCF paths; do not use "
                "--vcf-path or --vcf-dir"
            )
        if table_input.single_sample:
            if args.vcf_path is None:
                parser.error("--vcf-path is required for single-sample VAF validation")
            if args.vcf_dir is not None:
                parser.error("Use --vcf-path, not --vcf-dir, for single-sample VAF validation")
        else:
            if args.vcf_path is not None:
                parser.error("Use --vcf-dir, not --vcf-path, for multi-sample VAF validation")
            if args.vcf_dir is None and not manifest_vcf:
                parser.error("--vcf-dir is required for multi-sample VAF validation")
    elif args.vcf_path is not None or args.vcf_dir is not None:
        parser.error("--vcf-path/--vcf-dir require --vaf-validation true")

    if args.plot_dpi < 1:
        parser.error("--plot-dpi must be at least 1")

    needs_alignment = args.plot_cov or args.nn_validation
    if needs_alignment:
        if not math.isfinite(args.coverage_cap_multiplier) or args.coverage_cap_multiplier <= 0:
            parser.error("--coverage-cap-multiplier must be a positive finite number")
        if table_input.single_sample:
            if args.alignment_path is None:
                parser.error(
                    "--alignment-path is required for single-sample CNN validation "
                    "or coverage profiles"
                )
            if args.alignment_dir is not None:
                parser.error("Use --alignment-path, not --alignment-dir, for single-sample validation")
        else:
            if manifest_alignment:
                if args.alignment_dir is not None or args.alignment_path is not None:
                    parser.error(
                        "The multi-sample manifest already defines BAM/CRAM paths; "
                        "do not use --alignment-dir or --alignment-path"
                    )
            elif args.alignment_dir is None:
                parser.error(
                    "--alignment-dir is required for multi-sample CNN validation "
                    "or coverage profiles"
                )
            elif args.alignment_path is not None:
                parser.error("Use --alignment-dir, not --alignment-path, for multi-sample validation")
    else:
        if args.alignment_path is not None or args.alignment_dir is not None:
            parser.error(
                "--alignment-path/--alignment-dir require --plot-cov true or "
                "--cnn-validation true"
            )
    if args.nn_model is not None and not args.nn_validation:
        parser.error("--nn-model requires --cnn-validation true")
    if args.min_variants < 2:
        parser.error("--min-variants must be at least 2")
    if args.precision < 1 or args.precision > 15:
        parser.error("--precision must be between 1 and 15")

    if args.vcf_dir is not None:
        args.vcf_dir = args.vcf_dir.resolve()
    if args.vcf_path is not None:
        args.vcf_path = args.vcf_path.expanduser().resolve()
    args.qual_cutoff = None if args.qual_cutoff < 0 else args.qual_cutoff
    if str(args.filter_name).upper() in {"ANY", "NONE", "*"}:
        args.filter_name = None

    try:
        vprint("Resolving per-sample output directories")
        sample_output_dirs = {
            str(run): resolve_sample_output_dir(args.output, str(run))
            for run in frame["Run"].drop_duplicates().astype(str)
        }
        if len(set(sample_output_dirs.values())) != len(sample_output_dirs):
            raise ToolError(
                "Two sample IDs resolve to the same output directory; adjust --output"
            )
        annotation_path = resolve_annotation_bed(args.centromere_telomere_bed)
        vprint(f"Loading centromere/telomere annotations from {annotation_path}")
        annotation_features = read_annotation_bed(annotation_path)
    except ToolError as exc:
        parser.error(str(exc))

    vprint("Annotating CNV intervals with centromere/telomere overlaps")
    annotations = annotate_special_regions(frame, annotation_features)
    coverage_paths = pd.Series(pd.NA, index=frame.index, dtype="string")
    if args.plot_cov:
        vprint("Starting full coverage profile stage")
        try:
            coverage_paths = run_validation_coverage_profiles(
                frame,
                args,
                single_sample=table_input.single_sample,
            )
        except ToolError as exc:
            eprint(f"Coverage profile error: {exc}")
    else:
        vprint("Full coverage profile stage is disabled")

    vaf_result: pd.DataFrame | None = None
    train_dimension_paths = pd.Series(pd.NA, index=frame.index, dtype="string")
    if args.vaf_validation:
        vprint("Starting VAF validation stage")
        vaf_result = run_vaf_validation(frame, args)
        if args.plot_vaf_train_dimension:
            vprint("Starting VAF train-space plot stage")
            try:
                train_dimension_paths = run_vaf_train_dimension_plots(
                    frame,
                    vaf_result,
                    args,
                )
            except ToolError as exc:
                eprint(f"VAF train-space plot error: {exc}")
        else:
            vprint("VAF train-space plot stage is disabled")
    else:
        vprint("VAF validation stage is disabled")

    nn_result: pd.DataFrame | None = None
    if args.nn_validation:
        vprint("Starting CNN validation stage")
        with tempfile.TemporaryDirectory(prefix="vacca-cnn-") as temporary_dir:
            try:
                nn_result = run_nn_validation(
                    frame,
                    args,
                    single_sample=table_input.single_sample,
                    inference_root=Path(temporary_dir),
                )
            except ToolError as exc:
                eprint(f"CNN validation error: {exc}")
                nn_result = empty_nn_result(frame)
    else:
        vprint("CNN validation stage is disabled")

    vprint("Building output table(s)")
    tables: dict[str, pd.DataFrame] = {}
    if vaf_result is not None:
        tables["OUTPUT_VAF.tsv"] = build_vaf_output(
            frame,
            vaf_result,
            annotations,
            coverage_paths,
            train_dimension_paths,
            single_sample=table_input.single_sample,
        )
    if nn_result is not None:
        tables["OUTPUT_CNN.tsv"] = build_cnn_output(
            frame,
            nn_result,
            annotations,
            coverage_paths,
            single_sample=table_input.single_sample,
        )
    if vaf_result is not None and nn_result is not None:
        tables["OUTPUT_COMBINED.tsv"] = build_combined_output(
            frame,
            vaf_result,
            nn_result,
            single_sample=table_input.single_sample,
        )

    vprint(
        "Writing output table(s): " + ", ".join(tables)
    )
    write_sample_outputs(
        frame,
        args.output,
        tables,
        precision=args.precision,
    )
    vaf_predicted_count = (
        int(vaf_result["is_val"].notna().sum()) if vaf_result is not None else 0
    )
    nn_predicted_count = (
        int(nn_result["NN_is_val"].notna().sum()) if nn_result is not None else 0
    )
    vprint(
        f"Completed {len(frame)} CNV row(s); VAF predictions: {vaf_predicted_count}; "
        f"CNN predictions: {nn_predicted_count}"
    )
    return 0 if vaf_predicted_count > 0 or nn_predicted_count > 0 else 3

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        eprint("Interrupted")
        raise SystemExit(130)
