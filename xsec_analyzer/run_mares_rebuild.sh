#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; mkdir -p "$LOG"
OBS=(pmu ppi costhmu costhpi thmupi)
for s in plus minus; do
  echo "########## mares $s univmake rebuild ##########"
  pids=()
  for o in "${OBS[@]}"; do
    FPM=configs/file_properties_numi_mares_${s}.txt \
    BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT=$PROC/ccpi_Run1_${o}_mares${s}_univmake.root \
    ./run_universe_maker.sh > "$LOG/mares${s}_${o}_univ.log" 2>&1 &
    pids+=($!)
  done
  fail=0; i=0
  for o in "${OBS[@]}"; do wait "${pids[$i]}" && echo "  univmake[$s/$o] OK" || { echo "  univmake[$s/$o] FAILED"; fail=1; }; i=$((i+1)); done
  echo "########## mares $s DONE (fail=$fail) ##########"
done
echo "########## MARES REBUILD ALLDONE ##########"
