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
  # Examples:
  #   ![alt](figures_01/foo.png)  ->  ![alt](01-Free_Energy_Structural_Relaxations/figures_01/foo.png)
  #   ![alt](05_raman_ir/foo.png) ->  ![alt](03-Spectral_Functions_Raman_IR_Spectra/05_raman_ir/foo.png)
  # We handle both patterns: figures_NN/ and 05_raman_ir/
  sed -E "s#(!\[.*\]\()(figures_|05_raman_ir/)#\1${dir}/\2#g" "$md" >> "$TMP_MD"
done

echo "Generating PDF with pandoc..."
echo "  Output: $OUTPUT"

pandoc "$TMP_MD" \
  --pdf-engine=xelatex \
  --include-in-header=header.tex \
  -o "$OUTPUT" \
  --toc \
  --toc-depth=2 \
  --highlight-style=tango \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -V geometry:margin=2.5cm \
  -V mainfont="DejaVu Serif" \
  -V monofont="DejaVu Sans Mono" \
  --metadata title="SSCHA School Tutorials" \
  --metadata subtitle="Hands-on sessions on the Stochastic Self-Consistent Harmonic Approximation"

rm -f "$TMP_MD"

echo "Done: $OUTPUT"
