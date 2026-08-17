#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
echo "===== M_A^RES FAKE-DATA CLOSURE (unfolded vs A_C-smeared shifted truth) ====="
for s in plus minus; do echo "--- M_A^RES $s ---"
  for o in pmu ppi costhmu costhpi thmupi; do
    r=$(bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_mares${s}.txt configs/ccpi_${o}_slice_config_opt.txt 2>/dev/null | grep -iE "Fakedata:" | head -1)
    printf "  %-8s %s\n" "$o" "${r:-<no closure line>}"
  done
done
echo "===== CLOSURE DONE ====="
