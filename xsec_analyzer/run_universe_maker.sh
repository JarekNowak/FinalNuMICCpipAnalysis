#!/usr/bin/env bash
# run_universe_maker.sh
# Stage 2 of the xsec_analyzer chain: build the systematic-universe histograms
# from the processed xsec-ana-*.root files (configs/file_properties_numi.txt)
# using a bin-scheme / univmake config, writing one output ROOT file.
#
# univmake CLI:  univmake LIST_FILE BIN_CONFIG OUTPUT [FPM_CONFIG]
#   LIST_FILE : one input file per line (first whitespace token used).  We reuse
#               file_properties_numi.txt here — univmake reads only the path, and
#               passing it again as FPM_CONFIG gives the data POT/trigger info
#               needed for the total-event-count step (MCC9SystematicsCalculator).
#
# Usage:
#   ./run_universe_maker.sh                       # default pmu config
#   BIN_CONFIG=configs/foo.txt ./run_universe_maker.sh
#   OUT=/path/univ.root ./run_universe_maker.sh

set -euo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"

FPM="${FPM:-configs/file_properties_numi.txt}"
BIN_CONFIG="${BIN_CONFIG:-configs/ccpi_pmu_bin_config.txt}"
OUT="${OUT:-/data/uboone/processed/ccpi_Run1_pmu_univmake.root}"

set +u; source ./setup_xsec_analyzer.sh; set -u

command -v root &>/dev/null || { echo "ERROR: ROOT not in PATH" >&2; exit 1; }

if [[ ! -x bin/univmake ]]; then
  echo "[build] bin/univmake not found — running make ..."
  make -j"$(nproc)" bin/univmake
fi

[[ -f "$FPM" ]]        || { echo "ERROR: $FPM not found" >&2; exit 1; }
[[ -f "$BIN_CONFIG" ]] || { echo "ERROR: $BIN_CONFIG not found" >&2; exit 1; }

# Verify every processed input listed in the FPM file exists before starting
missing=0
while read -r path _rest; do
  [[ -z "${path:-}" || "$path" == \#* ]] && continue
  [[ -f "$path" ]] || { echo "  MISSING processed input: $path" >&2; missing=1; }
done < "$FPM"
[[ $missing -eq 0 ]] || { echo "ERROR: run ProcessNTuples (run_process_ntuples.sh) first." >&2; exit 1; }

echo "[info] FPM / list  : $FPM"
echo "[info] bin config  : $BIN_CONFIG"
echo "[info] output      : $OUT"
mkdir -p "$(dirname "$OUT")"

time univmake "$FPM" "$BIN_CONFIG" "$OUT" "$FPM"

echo
echo "[done] universe file: $OUT"
