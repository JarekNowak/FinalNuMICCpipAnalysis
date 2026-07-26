#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; mkdir -p "$LOG"
OBS=(pmu ppi costhmu costhpi thmupi)
echo "########## univmake rebuild with beamline (opt bins, parallel) ##########"
pids=()
for o in "${OBS[@]}"; do
  FPM=configs/file_properties_numi_beamline.txt \
  BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
  OUT=$PROC/ccpi_Run1_${o}_opt_univmake.root \
  ./run_universe_maker.sh > "$LOG/beamline_${o}_univ.log" 2>&1 &
  pids+=($!); echo "  launched univmake[$o] pid=${pids[-1]}"
done
fail=0; i=0
for o in "${OBS[@]}"; do
  if wait "${pids[$i]}"; then echo "  univmake[$o] OK"; else echo "  univmake[$o] FAILED"; fail=1; fi
  i=$((i+1))
done
echo "########## BEAMLINE REBUILD DONE (fail=$fail) ##########"
