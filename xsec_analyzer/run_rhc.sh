#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed; mkdir -p unfold_output
echo "[wait] RHC throw..."; for i in $(seq 1 60); do grep -q "RHC seed" ../logs/throw_rhc.log 2>/dev/null && break; sleep 20; done
grep -E "RHC seed" ../logs/throw_rhc.log
echo "[build] RHC univmakes..."
pids=()
for o in pmu ppi costhmu costhpi thmupi; do
  FPM=configs/file_properties_numi_rhc.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT=$PROC/ccpi_RHC_${o}_univmake.root ./run_universe_maker.sh > "../logs/rhc_${o}_univ.log" 2>&1 &
  pids+=($!)
done
fail=0; i=0; for o in pmu ppi costhmu costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || { echo "  univmake[$o] FAIL"; fail=1; }; i=$((i+1)); done
echo "########## RHC BUILD DONE (fail=$fail) ##########"
echo "===== RHC FLUX BREAKDOWN + CLOSURE ====="
for o in pmu ppi costhmu costhpi thmupi; do
  fb=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_rhc.txt "$o" 2>/dev/null | grep "flux RECO")
  cl=$(bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_rhc.txt configs/ccpi_${o}_slice_config_opt.txt /tmp/rhc_${o}.root 2>/dev/null | grep -iE "^truth: " | head -1)
  echo "--- $o ---"; echo "  $fb"; echo "  closure $cl"
done
echo "===== RHC ALLDONE ====="
