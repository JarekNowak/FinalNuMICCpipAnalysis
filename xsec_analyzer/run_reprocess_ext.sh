#!/usr/bin/env bash
# run_reprocess_ext.sh — reprocess every beam-off (EXT) sample with the swtrig fail-open
# fix (commit 24c1682). The CRT-era Run4/Run5 productions carry a vestigial swtrig==0 that
# previously rejected 100% of their events, so those four selected zero cosmics; the
# pre-CRT samples are unaffected by the fix and are reprocessed only so the whole set is
# produced by one build.
#
# Writes to a STAGING directory. Nothing live is overwritten -- promote only after the
# per-sample cosmic rates have been checked against each other.
#   usage: ./run_reprocess_ext.sh [NPAR]      (default 3)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"

RAW=/data/uboone/EXT
STAGE=/data/uboone/processed/ext_reprocessed
LOG=../logs
NPAR=${1:-3}
mkdir -p "$STAGE"
rm -f "$LOG/reprocess_ext.status"

# run4c is already done (0 -> 3 selected) but is reprocessed here too so the whole set
# comes from one build; skip it by removing it from this list if time is short.
SAMPLES=(
  neutrinoselection_filt_run1_beamoff
  neutrinoselection_filt_run3b_beamoff
  numi_pelee_ntuple_beam_off_run4a_rhc_ana
  numi_pelee_ntuple_beam_off_run4b_rhc_ana
  numi_pelee_ntuple_beam_off_run4d_fhc_ana
  numi_pelee_ntuple_beam_off_run5_fhc_ana
  beamoff_run1Andrun3
)

one() {
  local s=$1
  local out="$STAGE/xsec-ana-${s}.root"
  if [ "$(stat -c%s "$out" 2>/dev/null || echo 0)" -gt 10000000 ]; then
    echo "[skip] $s (already staged)"; echo "DONE_$s" >> "$LOG/reprocess_ext.status"; return 0
  fi
  local list; list=$(mktemp "$STAGE/.list_${s}_XXXX.txt")
  echo "$RAW/${s}.root extBNB" > "$list"
  echo "[$(date +%H:%M)] START $s"
  OUT_DIR="$STAGE" SELECTION=CC1mu1piXp NTUPLE_LIST="$list" \
    nice -n 12 ./run_process_ntuples.sh > "$LOG/reprocess_ext_${s}.log" 2>&1
  local rc=$?
  rm -f "$list"
  echo "[$(date +%H:%M)] END   $s rc=$rc ($(( $(stat -c%s "$out" 2>/dev/null || echo 0) / 1048576 ))MB)"
  echo "DONE_$s" >> "$LOG/reprocess_ext.status"
}

echo "==== EXT REPROCESS START $(date) NPAR=$NPAR ===="
for s in "${SAMPLES[@]}"; do
  one "$s" &
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done
done
wait
echo "ALL_DONE" >> "$LOG/reprocess_ext.status"
echo "==== EXT REPROCESS DONE $(date) ===="
