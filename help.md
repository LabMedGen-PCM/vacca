# CNV VAF Validator 1.0.0 — command reference

This document describes `cnv_vaf_validator.py` version **1.0.0**.

## Installation

```bash
bash install.sh
mamba activate cnv-vaf-validator
```

The installer creates or updates the `cnv-vaf-validator` Conda environment and installs the equivalent commands `cnv-vaf-validator` and `cnv-vaf-validator.py`. See `readme.md` for installation prerequisites, manual environment creation, and troubleshooting.

## Built-in test run

```bash
cnv-vaf-validator.py --test-run
```

The test runs VAF and CNN validation on the bundled ClinCNV, VCF, and CRAM data and writes:

```text
test_run/
└── single_sample/
    ├── OUTPUT_VAF.tsv
    ├── OUTPUT_CNN.tsv
    ├── OUTPUT_COMBINED.tsv
    ├── coverage_profiles/
    └── vaf_train_dimension/{DEL,DIP,DUP}/
```

## Validation switches

```text
--verbose true|false          default: true
--vaf-validation true|false   default: true
--cnn-validation true|false   default: false
--plot-cov true|false         default: true
--plot-vaf-train-dimension true|false  default: true
```

Verbose mode reports the current input, manifest row, sample, coverage task,
model stage, and output file to standard error. Set `--verbose false` to hide
progress messages. Warnings and errors are always printed.

At least one validation model must be enabled in `--mode validate`. The legacy spellings `--NN-validation` and `--nn-validation` remain aliases for `--cnn-validation`.

`--plot-cov` controls full user-facing profiles during validation. The separate `--mode plot-cov` still creates profiles without running either model.

`--plot-vaf-train-dimension` controls per-CNV JPEGs in the bundled VAF training projection. It applies only when VAF validation is enabled, is independent of `--plot-cov`, and requires no BAM/CRAM.

## Minimal commands

### Multi-sample VAF with default coverage profiles

```bash
cnv-vaf-validator.py \
  --cnvs CALLS.tsv \
  --vcf-dir VCF_DIR \
  --alignment-dir ALIGNMENT_DIR \
  --output OUTPUT_DIR
```

### Multi-sample VAF without profiles

```bash
cnv-vaf-validator.py \
  --cnvs CALLS.tsv \
  --vcf-dir VCF_DIR \
  --output OUTPUT_DIR \
  --plot-cov false
```

### CNN only

```bash
cnv-vaf-validator.py \
  --cnvs CALLS.tsv \
  --vaf-validation false \
  --cnn-validation true \
  --alignment-dir ALIGNMENT_DIR \
  --output OUTPUT_DIR
```

### VAF and CNN

```bash
cnv-vaf-validator.py \
  --cnvs CALLS.tsv \
  --vcf-dir VCF_DIR \
  --alignment-dir ALIGNMENT_DIR \
  --cnn-validation true \
  --output OUTPUT_DIR
```

### Multi-sample ClinCNV/Manta manifest

```bash
cnv-vaf-validator.py --cnvs samples.tsv --cnn-validation true --output OUTPUT_DIR
```

The manifest supplies CNV, VCF, and BAM/CRAM paths.

### Single-sample validation

```bash
cnv-vaf-validator.py \
  --cnvs CALLS_OR_CALLER_OUTPUT \
  --vcf-path SAMPLE.small_variants.vcf.gz \
  --alignment-path SAMPLE.cram \
  --output OUTPUT_DIR
```

For CNN-only, omit `--vcf-path` and set `--vaf-validation false --cnn-validation true`. For VAF-only without profiles, omit `--alignment-path` and set `--plot-cov false`.

### One explicit CNV

```bash
cnv-vaf-validator.py \
  --chr chr1 --start 100000 --end 200000 --svtype DEL \
  --vcf-path SAMPLE.small_variants.vcf.gz \
  --alignment-path SAMPLE.cram \
  --output OUTPUT_DIR
```

### Coverage only

```bash
cnv-vaf-validator.py \
  --mode plot-cov \
  --cnvs CALLS.tsv \
  --alignment-dir ALIGNMENT_DIR \
  --plot-dir coverage_plots
```

## CNV input auto-detection

`--cnvs` is detected by content as a standard table, a native ClinCNV table, a Manta VCF, or a multi-sample native-output manifest.

Standard multi-sample table:

```text
sample_id  chr  start  end  SVTYPE
```

Standard single-sample table:

```text
chr  start  end  SVTYPE
```

Fields are positional; header names are ignored. Coordinates are 1-based inclusive. Legacy tables with `length` immediately before `SVTYPE` remain accepted.

Only `DEL`, `DIP`, and `DUP` are supported. `INS` is no longer accepted. Manta `INS`, `BND`, and `INV` records are skipped because the models do not contain those classes.

### Multi-sample ClinCNV/Manta manifest

```text
sample_id | cnvs | bam/cram | vcf
SAMPLE_1 | calls/SAMPLE_1.clincnv.tsv | bam/SAMPLE_1.cram | vcf/SAMPLE_1.vcf.gz
SAMPLE_2 | calls/SAMPLE_2.manta.vcf.gz | bam/SAMPLE_2.bam | vcf/SAMPLE_2.vcf.gz
```

Relative paths are resolved relative to the manifest. A missing caller output skips that manifest row. Missing VCF or alignment files affect only the corresponding validation stage and sample.

### ClinCNV sex handling

| Sample/chromosome | Expected CN | Classification |
|---|---:|---|
| autosomes | 2 | `<2 DEL`, `=2 DIP`, `>2 DUP` |
| female X | 2 | `<2 DEL`, `=2 DIP`, `>2 DUP` |
| female Y | 0 | `=0 DIP`, `>0 DUP` |
| male X or Y | 1 | `<1 DEL`, `=1 DIP`, `>1 DUP` |

If `##gender of sample: F|M` is missing or unrecognized, sex chromosomes use baseline 2 with a warning.

## Output directory resolution

`--output` is required for validation.

Without a placeholder:

```bash
--output /path/to/output_dir
```

Result:

```text
/path/to/output_dir/SAMPLE_ID/
├── OUTPUT_VAF.tsv
├── OUTPUT_CNN.tsv
├── OUTPUT_COMBINED.tsv
├── coverage_profiles/{DEL,DIP,DUP}/...
└── vaf_train_dimension/{DEL,DIP,DUP}/SAMPLE_ID_chr_start_end.jpeg
```

With the literal `{sample_id}` placeholder:

```bash
--output "/path/to/{sample_id}/output_dir"
```

For `SAMPLE_1`, result files are written directly to `/path/to/SAMPLE_1/output_dir/`; no additional sample directory is added.

In single-sample mode the internal directory token is `single_sample`, and `SAMPLE_ID` is omitted from the tables.

## Output files by mode

| VAF | CNN | Profiles | Files |
|---:|---:|---:|---|
| true | false | true | `OUTPUT_VAF.tsv` plus profiles |
| false | true | true | `OUTPUT_CNN.tsv` plus profiles |
| true | true | true | VAF, CNN, combined, plus profiles |
| true | false | false | `OUTPUT_VAF.tsv` |
| false | true | false | `OUTPUT_CNN.tsv` |
| true | true | false | VAF, CNN, and combined |
| false | false | any | Argument error |
| standalone `plot-cov` | — | true | Profiles only |

When VAF validation is enabled, `--plot-vaf-train-dimension true` additionally creates train-space JPEGs regardless of the Profiles column above.

## OUTPUT_VAF.tsv

Multi-sample columns, in order:

```text
SAMPLE_ID
CHR
START
END
INPUT_SVTYPE
PREDICTED_SVTYPE
INPUT_OUTPUT_CONSISTENCY
MODEL_ELIGIBILITY
LOW_MAPPABILITY_REGION
INPUT_REGION_COVERAGE
NUM_OF_VARIANTS
MARKER_0.33
MARKER_0.66
MARKER_0.33+0.66
MARKER_0.5
MARKER_1
MODEL_SCORE_DEL
MODEL_SCORE_DIP
MODEL_SCORE_DUP
CNV_POSITION_IN_TRAIN_DIMENSION
```

- `MODEL_ELIGIBILITY`: `ELIGIBLE`, `INSUFFICIENT_NUM_VAR`, `ERROR`.
- `PREDICTED_SVTYPE`: `DEL`, `DIP`, `DUP`, `NA`.
- `INPUT_OUTPUT_CONSISTENCY`: `CONSISTENT`, `INCONSISTENT`, `NOT_CNV`, `NA`.
- `LOW_MAPPABILITY_REGION`: `CENTROMERIC`, `TELOMERIC`, `CENTROMERIC;TELOMERIC`, `NA`.
- `MARKER_0.33+0.66`: sum of `MARKER_0.33` and `MARKER_0.66`.
- `CNV_POSITION_IN_TRAIN_DIMENSION`: absolute path to the per-CNV train-space JPEG, or `NA` when disabled or unavailable.
- SVC scores are probabilities in `[0, 1]` or `NA`.

## OUTPUT_CNN.tsv

```text
SAMPLE_ID
CHR
START
END
INPUT_SVTYPE
PREDICTED_SVTYPE
INPUT_OUTPUT_CONSISTENCY
MODEL_ELIGIBILITY
LOW_MAPPABILITY_REGION
INPUT_REGION_COVERAGE
NN_PROBABILITY_DEL
NN_PROBABILITY_DIP
NN_PROBABILITY_DUP
```

CNN `MODEL_ELIGIBILITY` is `ELIGIBLE` or `ERROR`. Probabilities are in `[0, 1]` or `NA`.

## OUTPUT_COMBINED.tsv

Created only when both models are enabled:

```text
SAMPLE_ID
CHR
START
END
INPUT_SVTYPE
VAF_PREDICTED_SVTYPE
NN_PREDICTED_SVTYPE
FINAL_DECISION
```

`FINAL_DECISION` is `DEL`, `DIP`, or `DUP` when every available prediction agrees with the input; `CNV_not_input type` when an available prediction conflicts with the input; and `NA` when neither model produced a prediction.

`INPUT_REGION_COVERAGE` is an absolute path to the full public profile. It is `NA` when profiles are disabled or the profile failed. Legend-free CNN inputs are temporary internal files and never appear in the output tables.

## Models and VAF extraction

The bundled `models/svc_vaf.pkl` contains:

```python
SVC(
    C=1.0,
    class_weight=None,
    coef0=0.0,
    gamma="scale",
    kernel="sigmoid",
    probability=True,
)
```

SVC feature order is `X0.5`, `X1`, `X0.33_X0.66`, where `X0.33_X0.66 = X0.33 + X0.66`; encoded classes are `0 = DEL`, `1 = DIP`, `2 = DUP`.

The bundled `data/VAF_train.tsv` supplies the `SVTYPE`, `X0.33_X0.66`, `X0.5`, and `X1` values used for train-space plots. Each plot shows `X0.33_X0.66` on the X axis and `X0.5` on the Y axis. Training classes are labelled `train_DEL`, `train_DIP`, and `train_DUP`; the black CNV point is annotated as `INPUT_SVTYPE_chr_start_end`.

Small-variant VCF FORMAT values are tried in this order: `VAF`, `AF`, then `AD`. For `AD`, VAF is `AD[1] / sum(AD)`.

## Coverage behavior

BAM/CRAM files must be coordinate-sorted and indexed. CRAM decoding must have access to the matching reference. Full profiles include axes, grid, legend, and CNV boundary lines by default.

For very large regions, per-base rendering automatically falls back to median coverage in 100-bp windows, then 1,000 bp, 10,000 bp, and subsequent powers of ten.

## Argument reference

### General

| Argument | Meaning |
|---|---|
| `--test-run` | Run the bundled complete example. |
| `--mode {validate,plot-cov}` | Operation mode; default `validate`. |
| `-c, --cnvs FILE` | Standard table, manifest, ClinCNV output, or Manta VCF. |
| `--input-sep VALUE` | `auto`, `tab`, `tsv`, `comma`, `csv`, or a literal delimiter. |
| `-t, --threads N` | Parallelism limit; integer `>=1`. |
| `--version` | Print version `1.0.0`. |

### Validation

| Argument | Meaning |
|---|---|
| `--vaf-validation true|false` | Enable VAF SVC; default `true`. |
| `--plot-vaf-train-dimension true|false` | Create per-CNV VAF train-space JPEGs; default `true`. |
| `--cnn-validation true|false` | Enable CNN; default `false`. Aliases: `--NN-validation`, `--nn-validation`. |
| `--plot-cov true|false` | Create full public profiles; default `true`. |
| `-o, --output PATH` | Output directory or mask containing `{sample_id}`. |
| `-v, --vcf-dir DIR` | Multi-sample small-variant VCF root. |
| `-m, --vcf-mask MASK` | Default `{sample}/{sample}.output.vcf.gz`. |
| `--vcf-path FILE` | Single-sample small-variant VCF. |
| `--backend {auto,bcftools,python}` | VAF extraction backend. |
| `--qual-cutoff FLOAT` | Strict `QUAL > cutoff`; negative disables. |
| `--filter TOKEN` | Required FILTER; `ANY`, `NONE`, or `*` disables. |
| `--min-variants N` | Minimum VAF count; default 50, minimum 2. |
| `--precision N` | Output decimals, 1–15. |
| `--centromere-telomere-bed FILE` | Override bundled hg38 annotation. |
| `--nn-model FILE` | Override `models/CNN_coverage.pt`. |
| `--nn-batch-size N` | CNN inference batch size; default 4. |

### Explicit interval

| Argument | Meaning |
|---|---|
| `--chr CHROM` | Chromosome/contig. |
| `--start START` | 1-based inclusive start. |
| `--end END` | 1-based inclusive end. |
| `--svtype DEL|DIP|DUP` | Input CNV type. |

### Alignment and plotting

| Argument | Meaning |
|---|---|
| `--alignment-path FILE` | Single BAM/CRAM; aliases `--bam-path`, `--cram-path`. |
| `--alignment-dir DIR` | Multi-sample alignment root; aliases `--bam-dir`, `--cram-dir`. |
| `--alignment-mask MASK` | Default `{sample}.markdup.sorted.cram`. |
| `--plot-dir DIR` | Destination for standalone `--mode plot-cov`. |
| `--plot-format jpeg|svg` | Full-profile image format. |
| `--coverage-cap-multiplier FLOAT` | Positive Y-axis cap multiplier; default 5.0. |
| `--plot-dpi N` | JPEG resolution; default 200. |
| `--legend true|false` | Show axes, ticks, labels, and grid; default `true`. |
| `--cnv-borders true|false` | Draw CNV boundaries; default `true`. |

Run `cnv-vaf-validator.py --help` for argparse-generated usage.
