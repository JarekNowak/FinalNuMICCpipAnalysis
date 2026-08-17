#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
echo "[wait] fake-data build (nue POT)..."
for i in $(seq 1 400); do grep -q "DONE build_fakedata" ../logs/build_fakedata.log 2>/dev/null && break; sleep 30; done
grep -q "DONE build_fakedata" ../logs/build_fakedata.log || { echo TIMEOUT_BUILD; exit 1; }
grep -E "rate change|wrote" ../logs/build_fakedata.log
echo "[rebuild] mares univmakes (nue POT)..."; ./run_mares_rebuild.sh 2>&1 | grep -E "univmake\[|DONE"
OBS=(pmu ppi costhmu costhpi thmupi)
echo "===== FLUX BREAKDOWN (nue POT 6.626e20) : 2nd vs 1st derivative ====="
for s in plus minus; do
  for tag in "2nd:" "1st:_1d"; do d="${tag#*:}"; lbl="${tag%%:*}"
    echo "--- M_A^RES $s / $lbl-deriv ---"
    for o in "${OBS[@]}"; do
      /tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_mares${s}${d}.txt "$o" 2>/dev/null | grep "flux RECO"
    done
  done
done
echo "===== CLOSURE (unfolded vs A_C-smeared M_A^RES truth) : 2nd vs 1st derivative ====="
for s in plus minus; do
  for tag in "2nd:" "1st:_1d"; do d="${tag#*:}"; lbl="${tag%%:*}"
    echo "--- M_A^RES $s / $lbl-deriv ---"
    for o in "${OBS[@]}"; do
      r=$(bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_mares${s}${d}.txt configs/ccpi_${o}_slice_config_opt.txt 2>/dev/null | grep -iE "Fakedata:" | head -1)
      printf "  %-8s %s\n" "$o" "${r:-<none>}"
    done
  done
done
echo "===== MARES V2 ALLDONE ====="
