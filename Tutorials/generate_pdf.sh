#!/usr/bin/env bash
set -euo pipefail

OUTPUT="SSCHA_Tutorials.pdf"
TMP_MD=$(mktemp /tmp/sscha_tutorials_XXXXXX.md)

echo "Building combined markdown with adjusted paths..."

> "$TMP_MD"

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
  --lua-filter=boxed-blockquotes.lua \
  -o "$OUTPUT" \
  --toc \
  --toc-depth=3 \
  -V documentclass=report \
  -V numbersections \
  -V secnumdepth=2 \
  --syntax-highlighting=idiomatic \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -V geometry:margin=2.5cm \
  --metadata title="SSCHA School Tutorials" \
  --metadata subtitle="Hands-on sessions on the Stochastic Self-Consistent Harmonic Approximation"

rm -f "$TMP_MD"

echo "Done: $OUTPUT"
