#!/bin/bash
# eps2pdf.sh <file.eps> [...] -- convert to a tightly-cropped PDF, removing the EPS.
#
# Why not have ROOT write the PDF directly: ROOT prints onto a full A4 page and tags
# landscape canvases with /Rotate 90, so \includegraphics scales the surrounding whitespace
# and the plot lands far smaller than the requested width. Cropping that PDF afterwards is
# unreliable because Ghostscript's PageOffset applies in unrotated space while the bounding
# box is reported in rotated space, which displaces the content.
#
# ROOT's EPS output, by contrast, carries a correct tight %%BoundingBox, so -dEPSCrop gives
# a page equal to the plot. -dAutoRotatePages=/None keeps Ghostscript from reintroducing
# the rotation.
set -euo pipefail

for f in "$@"; do
  [ -f "$f" ] || { echo "  [skip] no such file: $f" >&2; continue; }
  out="${f%.eps}.pdf"
  gs -q -o "$out" -sDEVICE=pdfwrite -dEPSCrop -dAutoRotatePages=/None -f "$f"
  rm -f "$f"
  echo "  $out"
done
