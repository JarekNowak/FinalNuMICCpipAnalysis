#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
REBUILD_PID="$1"
echo "[wait] rebuild (pid $REBUILD_PID)..."
while kill -0 "$REBUILD_PID" 2>/dev/null; do sleep 30; done
# verify all 10 univmakes present and non-trivial
n=0; for s in plus minus; do for o in pmu ppi costhmu costhpi thmupi; do
  f=/data/uboone/processed/ccpi_Run1_${o}_mares${s}_univmake.root
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0); [ "$sz" -gt 10000000 ] && n=$((n+1)) || echo "  MISSING/SMALL: ${o}_${s} ($sz)"
done; done
echo "univmakes OK: $n/10"
echo "===== FLUX BREAKDOWN (M_A^RES fake data, real FHC data POT 8.817e20) ====="
for s in plus minus; do echo "--- M_A^RES $s ---"
  for o in pmu ppi costhmu costhpi thmupi; do
    /tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_mares${s}.txt "$o" 2>/dev/null | grep "flux RECO"
  done
done
./run_mares_closure.sh
echo "===== MARES FINAL ALLDONE ====="
