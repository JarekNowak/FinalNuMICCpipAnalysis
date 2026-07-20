#!/usr/bin/env bash
# run_unfolder.sh
# Stage 3 of the xsec_analyzer chain: extract the unfolded CCpi+ dsigma/dpmu
# cross section from the universe file, with full systematics.
#
# Two binaries share the same CrossSectionExtractor (identical systematics,
# unfolding, and NuMI flux conversion) but differ in output:
#   Unfolder      -> structured ROOT file (EventCountUnits/XsecUnits + covariances
#                    + per-slice histograms).  This is the measurement product.
#   UnfolderNuMI  -> diagnostic PDFs into unfold_output; writes NO ROOT file:
#                    the unfolded slices, the regularization matrix, and the
#                    analysis-step plots (plot_step1..4: reco spectrum,
#                    background subtraction, smearceptance, efficiency) drawn
#                    in the same style for side-by-side comparison.
# This driver runs the plain Unfolder for the ROOT file, then UnfolderNuMI for
# the plots (set PLOTS=0 to skip the plots step).
#
# CLI (both):  <bin> XSEC_CONFIG SLICE_CONFIG OUTPUT_FILE
#
# Usage:
#   ./run_unfolder.sh
#   PLOTS=0 ./run_unfolder.sh
#   XSEC=configs/foo.txt SLICE=configs/bar.txt OUT=/path/out.root ./run_unfolder.sh

set -euo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"

XSEC="${XSEC:-configs/ccpi_xsec_config_numi.txt}"
SLICE="${SLICE:-configs/ccpi_pmu_slice_config.txt}"
OUT="${OUT:-unfold_output/ccpi_Run1_pmu_xsec.root}"

set +u; source ./setup_xsec_analyzer.sh; set -u
command -v root &>/dev/null || { echo "ERROR: ROOT not in PATH" >&2; exit 1; }

PLOTS="${PLOTS:-1}"
need=(bin/Unfolder)
[[ "$PLOTS" == "1" ]] && need+=(bin/UnfolderNuMI)
for b in "${need[@]}"; do
  if [[ ! -x "$b" ]]; then
    echo "[build] $b not found — running make ..."
    make -j"$(nproc)" "$b"
  fi
done

for f in "$XSEC" "$SLICE"; do
  [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done

# Verify the universe file referenced by the XSEC config exists
univ=$(awk '$1=="UnivFile"{print $2}' "$XSEC")
[[ -f "$univ" ]] || { echo "ERROR: UnivFile $univ not found — run run_universe_maker.sh first." >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
# UnfolderNuMI hardcodes its text/plot dumps to unfold_output relative to CWD
# (which is this script's dir), i.e. the repo root's unfold_output.
mkdir -p unfold_output

echo "[info] xsec config : $XSEC"
echo "[info] slice config: $SLICE"
echo "[info] univ file   : $univ"
echo "[info] output      : $OUT"
echo "[info] plots       : $PLOTS (UnfolderNuMI -> unfold_output)"

echo "── Unfolder (ROOT output) ──"
time Unfolder "$XSEC" "$SLICE" "$OUT"

if [[ "$PLOTS" == "1" ]]; then
  echo "── UnfolderNuMI (diagnostic plots) ──"
  UnfolderNuMI "$XSEC" "$SLICE" "$OUT" || echo "  (plots step failed — non-fatal)"
fi

echo
if [[ -f "$OUT" ]]; then
  echo "[done] cross-section file: $OUT ($(du -h "$OUT" | cut -f1))"
else
  echo "[warn] no ROOT output produced at $OUT"
fi
