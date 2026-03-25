#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$SCRIPT_DIR/out"
OUT_PDF="$OUT_DIR/BaneMECP_Project_Manual_zh.pdf"

mkdir -p "$OUT_DIR"

pandoc "$SCRIPT_DIR/BaneMECP_Project_Manual_zh.md" \
  --defaults "$SCRIPT_DIR/pandoc.yaml" \
  --output "$OUT_PDF"

echo "[OK] wrote $OUT_PDF"
