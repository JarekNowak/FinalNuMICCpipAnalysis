#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
BLOG=../logs/build_fakedata.log
echo "[wait] fake-data build..."
for i in $(seq 1 400); do grep -q "DONE build_fakedata" "$BLOG" 2>/dev/null && break; sleep 30; done
grep -q "DONE build_fakedata" "$BLOG" || { echo "TIMEOUT build"; exit 1; }
grep -E "rate change|wrote" "$BLOG"
for s in plus minus; do
  f=/data/uboone/processed/xsec-ana-fakedata_MaRES_${s}_run15.root
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0); echo "  fake_$s = $((sz/1000000))MB"
  [ "$sz" -lt 1000000000 ] && { echo "ABORT: fake_$s too small"; exit 1; }
done
echo "[rebuild] mares univmakes..."; ./run_mares_rebuild.sh 2>&1 | grep -E "univmake\[|DONE"
echo "===== FLUX BREAKDOWN (high-stat MaRES fake data) ====="
for s in plus minus; do echo "--- M_A^RES $s ---"
  for o in pmu ppi costhmu costhpi thmupi; do
    /tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_mares${s}.txt "${o}" 2>/dev/null | grep "flux RECO"
  done
done
echo "===== MARES ALLDONE ====="
