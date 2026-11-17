#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$PREFIX/bin" "$PREFIX/share/doc/vacca" \
  "$PREFIX/share/vacca/models" "$PREFIX/share/vacca/data" \
  "$PREFIX/share/vacca/example"
install -m 0755 "$SRC_DIR/vacca.py" "$PREFIX/bin/vacca"
install -m 0644 "$SRC_DIR/readme.md" "$PREFIX/share/doc/vacca/readme.md"
install -m 0644 "$SRC_DIR/help.md" "$PREFIX/share/doc/vacca/help.md"
install -m 0644 "$SRC_DIR/hg38_centromeres_telomeres.bed" \
  "$PREFIX/share/vacca/hg38_centromeres_telomeres.bed"
install -m 0644 "$SRC_DIR/models/CNN_coverage.pt" \
  "$PREFIX/share/vacca/models/CNN_coverage.pt"
install -m 0644 "$SRC_DIR/models/svc_vaf.pkl" \
  "$PREFIX/share/vacca/models/svc_vaf.pkl"
install -m 0644 "$SRC_DIR/data/VAF_train.tsv" \
  "$PREFIX/share/vacca/data/VAF_train.tsv"
install -m 0644 "$SRC_DIR/example/"* \
  "$PREFIX/share/vacca/example/"
