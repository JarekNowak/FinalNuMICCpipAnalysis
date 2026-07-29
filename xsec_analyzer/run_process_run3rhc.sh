#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
RAW=/data/uboone/new_numi_flux; OUT=/data/uboone/processed; LOG=../logs
pids=()
for s in aa ab ac ad ae; do
  in="$RAW/Run3_rhc_new_numi_flux_rhc_pandora_ntuple_${s}.root"
  out="$OUT/xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_${s}.root"
  ( ProcessNTuples "$in" numuMC CC1mu1piXp "$out" > "$LOG/proc_run3rhc_${s}.log" 2>&1 && echo "  ${s} OK" || echo "  ${s} FAILED" ) &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "########## PROCESS RUN3RHC DONE (fail=$fail) ##########"
ls -la $OUT/xsec-ana-Run3_rhc_*_a*.root 2>/dev/null | awk '{print $5/1e9" GB  "$NF}'
