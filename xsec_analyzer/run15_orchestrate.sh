#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
DRV=../logs/process_run15_driver.log
# 1) wait for processing
echo "[wait] processing..."
for i in $(seq 1 480); do grep -q "PROCESS RUN15 DONE" "$DRV" 2>/dev/null && break; sleep 30; done
grep -qE "PROCESS RUN15 DONE" "$DRV" || { echo "TIMEOUT waiting for processing"; exit 1; }
echo "[proc] $(grep -E 'OK|FAILED|DONE' "$DRV")"
for r in Run2_fhc Run4_fhc Run4c_fhc Run4d_fhc; do
  f=/data/uboone/processed/xsec-ana-${r}_new_numi_flux_fhc_pandora_ntuple.root
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0); echo "  $r size=$((sz/1000000))MB"
  [ "$sz" -lt 100000000 ] && { echo "ABORT: $r too small"; exit 1; }
done
# 2) normalise summed_pot -> Sigma
echo "[setpot]"; root.exe -l -b -q macros/set_summed_pot.C 2>&1 | grep -E "true_pot|SIGMA|summed_pot=|scale|DONE|ABORT"
# 3) rebuild univmakes on combined MC
echo "[rebuild]"; ./run_run15_rebuild.sh 2>&1 | grep -E "launched|univmake\[|DONE"
# 4) flux breakdown
echo "===== FLUX BREAKDOWN (run15, combined FHC 1,2,4,4c,4d) ====="
for o in pmu ppi costhmu costhpi thmupi; do
  /tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_run15.txt "$o" 2>/dev/null | grep "flux RECO"
done
# 5) full syst breakdown
echo "===== SYST BREAKDOWN (run15) ====="
for o in pmu ppi costhmu costhpi thmupi; do
  echo "--- $o ---"; /tmp/syst_breakdown configs/ccpi_xsec_config_numi_${o}_run15.txt "${o}_r15" 2>/dev/null \
    | grep -E "Total|Flux \(PPFX\)|Detector|MC \+|Cross section"
done
echo "===== RUN15 ALLDONE ====="
