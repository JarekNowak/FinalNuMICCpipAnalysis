#!/usr/bin/env bash
# Process the additional FHC overlay runs (2,4,4c,4d) for the multi-run MC
# statistics boost. Run1 is already processed. Run5 excluded (corrupt SubRun POT).
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
RAW=/data/uboone/new_numi_flux; OUT=/data/uboone/processed; LOG=../logs
SEL=CC1mu1piXp
RUNS=(Run2_fhc Run4_fhc Run4c_fhc Run4d_fhc)
pids=()
for r in "${RUNS[@]}"; do
  in="$RAW/${r}_new_numi_flux_fhc_pandora_ntuple.root"
  out="$OUT/xsec-ana-${r}_new_numi_flux_fhc_pandora_ntuple.root"
  echo "[proc] $r -> $out"
  ( ProcessNTuples "$in" numuMC "$SEL" "$out" > "$LOG/proc_${r}.log" 2>&1 && echo "  $r OK" || echo "  $r FAILED" ) &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "########## PROCESS RUN15 DONE (fail=$fail) ##########"
ls -la $OUT/xsec-ana-Run{2,4,4c,4d}_fhc_new_numi_flux_fhc_pandora_ntuple.root 2>/dev/null | awk '{print $5/1e9" GB  "$NF}'
