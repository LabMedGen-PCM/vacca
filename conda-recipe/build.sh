#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$PREFIX/bin" "$PREFIX/share/doc/cnv-vaf-validator" \
  "$PREFIX/share/cnv-vaf-validator/models" "$PREFIX/share/cnv-vaf-validator/data" \
  "$PREFIX/share/cnv-vaf-validator/example"
install -m 0755 "$SRC_DIR/cnv_vaf_validator.py" "$PREFIX/bin/cnv-vaf-validator"
install -m 0755 "$SRC_DIR/cnv_vaf_validator.py" "$PREFIX/bin/cnv-vaf-validator.py"
install -m 0644 "$SRC_DIR/readme.md" "$PREFIX/share/doc/cnv-vaf-validator/readme.md"
install -m 0644 "$SRC_DIR/help.md" "$PREFIX/share/doc/cnv-vaf-validator/help.md"
install -m 0644 "$SRC_DIR/hg38_centromeres_telomeres.bed" \
  "$PREFIX/share/cnv-vaf-validator/hg38_centromeres_telomeres.bed"
install -m 0644 "$SRC_DIR/models/CNN_coverage.pt" \
  "$PREFIX/share/cnv-vaf-validator/models/CNN_coverage.pt"
install -m 0644 "$SRC_DIR/models/svc_vaf.pkl" \
  "$PREFIX/share/cnv-vaf-validator/models/svc_vaf.pkl"
install -m 0644 "$SRC_DIR/data/VAF_train.tsv" \
  "$PREFIX/share/cnv-vaf-validator/data/VAF_train.tsv"
install -m 0644 "$SRC_DIR/example/"* \
  "$PREFIX/share/cnv-vaf-validator/example/"
