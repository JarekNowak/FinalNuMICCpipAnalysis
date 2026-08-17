#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; mkdir -p unfold_output
echo "[wait] throw..."; for i in $(seq 1 400); do grep -q "DONE throw_fakedata" ../logs/throw_fakedata.log 2>/dev/null && break; sleep 20; done
grep -q "DONE throw_fakedata" ../logs/throw_fakedata.log || { echo TIMEOUT_THROW; exit 1; }
grep -E "kept|throwing" ../logs/throw_fakedata.log
echo "[rebuild]"; ./run_mares_rebuild.sh 2>&1 | grep -E "univmake\[|DONE"
OBS=(pmu ppi costhmu costhpi thmupi)
echo "===== FLUX BREAKDOWN (Poisson-thrown, nue POT) : 2nd vs 1st deriv ====="
for s in plus minus; do for tag in "2nd:" "1st:_1d"; do d="${tag#*:}"; lbl="${tag%%:*}"
  echo "--- $s / $lbl-deriv ---"
  for o in "${OBS[@]}"; do /tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_mares${s}${d}.txt "$o" 2>/dev/null | grep "flux RECO"; done
done; done
echo "===== CLOSURE (unfolded vs A_C-smeared M_A^RES truth) : 2nd vs 1st deriv ====="
for s in plus minus; do for tag in "2nd:" "1st:_1d"; do d="${tag#*:}"; lbl="${tag%%:*}"
  echo "--- $s / $lbl-deriv ---"
  for o in "${OBS[@]}"; do
    r=$(bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_mares${s}${d}.txt configs/ccpi_${o}_slice_config_opt.txt /tmp/cl_${o}_${s}${d}.root 2>/dev/null | grep -A6 "Key: Fakedata" | grep -iE "truth:" | head -1)
    printf "  %-8s %s\n" "$o" "${r:-<none>}"
  done
done; done
echo "===== MARES V3 ALLDONE ====="
