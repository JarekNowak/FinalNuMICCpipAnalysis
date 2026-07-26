#!/usr/bin/env bash
set -uo pipefail
cd "$( dirname "${BASH_SOURCE[0]}" )"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$PWD"
pids=()
for o in pmu ppi costhpi thmupi; do
  FPM=configs/file_properties_numi_beamline.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
  OUT=/data/uboone/processed/ccpi_Run1_${o}_optbl_univmake.root ./run_universe_maker.sh > ../logs/bl4_${o}.log 2>&1 &
  pids+=($!)
done
fail=0; i=0; for o in pmu ppi costhpi thmupi; do wait "${pids[$i]}" && echo "  $o OK" || { echo "  $o FAILED"; fail=1; }; i=$((i+1)); done
echo "BEAMLINE4 DONE (fail=$fail)"
