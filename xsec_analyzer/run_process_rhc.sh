#!/usr/bin/env bash
# Phase 1 of RHC expansion: process the 5 RHC overlays (Run1,2,4a,4b,4c) through
# ProcessNTuples with the CC1mu1piXp selection (charge-blind, so nubar signal is
# included). Outputs xsec-ana-*_rhc_*.root to /data/uboone/processed.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
RAW=/data/uboone/new_numi_flux; OUT=/data/uboone/processed; LOG=../logs; SEL=CC1mu1piXp
RUNS=(Run1_rhc Run2_rhc Run4a_rhc Run4b_rhc Run4c_rhc)
pids=()
for r in "${RUNS[@]}"; do
  in="$RAW/${r}_new_numi_flux_rhc_pandora_ntuple.root"
  out="$OUT/xsec-ana-${r}_new_numi_flux_rhc_pandora_ntuple.root"
  ( ProcessNTuples "$in" numuMC "$SEL" "$out" > "$LOG/proc_${r}.log" 2>&1 && echo "  $r OK" || echo "  $r FAILED" ) &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "########## PROCESS RHC DONE (fail=$fail) ##########"
ls -la $OUT/xsec-ana-Run*_rhc_*.root 2>/dev/null | awk '{print $5/1e9" GB  "$NF}'
