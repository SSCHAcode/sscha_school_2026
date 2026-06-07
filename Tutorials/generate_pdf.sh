#!/usr/bin/env bash
set -euo pipefail

# ---------- dependency checks ----------
missing_deps=()

if ! command -v pandoc &>/dev/null; then
  missing_deps+=("pandoc")
fi

if ! command -v xelatex &>/dev/null; then
  missing_deps+=("xelatex (part of texlive-xetex)")
fi

if ! command -v pandoc-crossref &>/dev/null; then
  missing_deps+=("pandoc-crossref")
fi

if [ ${#missing_deps[@]} -gt 0 ]; then
  echo "ERROR: The following dependencies are missing:" >&2
  for dep in "${missing_deps[@]}"; do
    echo "  - $dep" >&2
  done
  echo "" >&2
  echo "Please install them before running this script." >&2
  echo "" >&2
  echo "Install missing dependencies with:" >&2
  echo "  conda install -c conda-forge pandoc texlive-core pandoc-crossref" >&2
  exit 1
fi
# ----------------------------------------

OUTPUT="SSCHA_Tutorials.pdf"
TMP_MD=$(mktemp /tmp/sscha_tutorials_XXXXXX.md)

echo "Building combined markdown with adjusted paths..."

> "$TMP_MD"

setups_md="00-Setup/main.md"
if [ -f "$setups_md" ]; then
  echo "" >> "$TMP_MD"
  echo '\newpage' >> "$TMP_MD"
  echo "" >> "$TMP_MD"
  sed -E 's#(!\[.*\]\()([^)]+)\)#\100-Setup/\2)#g' "$setups_md" >> "$TMP_MD"
fi

for dir in 01-* 02-* 03-* 04-* 05-*; do
  md="$dir/main.md"
  if [ ! -f "$md" ]; then
    echo "Warning: $md not found, skipping."
    continue
  fi

  # Append a page break before each chapter (except the first)
  echo "" >> "$TMP_MD"
  echo '\newpage' >> "$TMP_MD"
  echo "" >> "$TMP_MD"

  # Rewrite image paths: prepend the subdirectory so pandoc finds them
  # Prepends the directory to any relative image path in ![alt](path)
  sed -E "s#(!\[.*\]\()([^)]+)\)#\1${dir}/\2)#g" "$md" >> "$TMP_MD"
done

echo "Generating PDF with pandoc..."
echo "  Output: $OUTPUT"

pandoc "$TMP_MD" \
  --pdf-engine=xelatex \
  --include-in-header=header.tex \
  --include-before-body=title.tex \
  --filter pandoc-crossref \
  --lua-filter=boxed-blockquotes.lua \
  -o "$OUTPUT" \
  --toc \
  --toc-depth=3 \
  -V documentclass=report \
  -V numbersections \
  -V secnumdepth=2 \
  --syntax-highlighting=breezeDark \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -V geometry:margin=2.5cm

rm -f "$TMP_MD"

echo "Done: $OUTPUT"
