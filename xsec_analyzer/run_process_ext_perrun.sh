#!/usr/bin/env bash
# run_process_ext_perrun.sh -- process every per-run beam-off (EXT) sample with BOTH
# selections (CC1mu1piXp,CC1mu1pi1p) into one staging directory, so a single file per
# run period can serve the inclusive (processed/) and proton-tagged (processed/w/) trees.
# Uses the swtrig fail-open build (24c1682), which the CRT-era Run4/Run5 samples need.
# Nothing live is overwritten; promote with the *_perrun file_properties afterwards.
#   usage: ./run_process_ext_perrun.sh [NPAR]   (default 3)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
RAW=/data/uboone/EXT; STAGE=/data/uboone/processed/ext_perrun; LOG=../logs/ext_perrun
NPAR=${1:-3}; mkdir -p "$STAGE" "$LOG"; : > "$LOG/status.txt"
SAMPLES=( neutrinoselection_filt_run1_beamoff neutrinoselection_filt_run3b_beamoff
  numi_pelee_ntuple_beam_off_run4a_rhc_ana numi_pelee_ntuple_beam_off_run4b_rhc_ana
  numi_pelee_ntuple_beam_off_run4c_fhc_ana numi_pelee_ntuple_beam_off_run4d_fhc_ana
  numi_pelee_ntuple_beam_off_run5_fhc_ana )
one() {
  local s=$1 out="$STAGE/xsec-ana-$1.root"
  echo "[$(date +%H:%M)] START $s"
  nice -n 12 ProcessNTuples "$RAW/$s.root" extBNB CC1mu1piXp,CC1mu1pi1p "$out" > "$LOG/$s.log" 2>&1
  local rc=$?
  echo "[$(date +%H:%M)] END   $s rc=$rc ($(( $(stat -c%s "$out" 2>/dev/null || echo 0) / 1048576 ))MB)"
  echo "$( [ $rc -eq 0 ] && echo OK || echo FAIL ) $s" >> "$LOG/status.txt"
}
echo "==== EXT PER-RUN PROCESS START $(date) NPAR=$NPAR ===="
for s in "${SAMPLES[@]}"; do one "$s" & while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done; done
wait; echo "ALL_DONE" >> "$LOG/status.txt"; echo "==== EXT PER-RUN PROCESS DONE $(date) ===="
