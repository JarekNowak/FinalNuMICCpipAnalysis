#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
RAW=/data/uboone/EXT; OUT=/data/uboone/processed; LOG=../logs
pids=()
for s in run4a_rhc run4b_rhc run4c_fhc run4d_fhc run5_fhc; do
  in="$RAW/numi_pelee_ntuple_beam_off_${s}_ana.root"
  out="$OUT/xsec-ana-numi_pelee_ntuple_beam_off_${s}_ana.root"
  ( ProcessNTuples "$in" extBNB CC1mu1piXp "$out" > "$LOG/proc_ext_${s}.log" 2>&1 && echo "  ${s} OK" || echo "  ${s} FAILED" ) &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "########## PROCESS EXT DONE (fail=$fail) ##########"
