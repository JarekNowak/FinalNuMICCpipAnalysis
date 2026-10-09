#!/bin/bash
#SBATCH --job-name=mp2picr
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=3000M
#SBATCH --time=08:00:00
#SBATCH --output=/data/uboone/processed/mp2pi_cr_skim/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mp2pi_cr_skim/logs/%x_%A_%a.out
# SLURM array: the two-pion selection in the opened control regions of the single-pion analysis
# (report/other_notes/multipion_note.tex, status). Processes ONE raw input (beam-on, overlay, beam-off or dirt) with
# the released single-pion selection and CC1mu2piBDTEvt (CC1mu2piBDT with its track dump), ALWAYS as the
# signal-region-stripped control-region skim (XSEC_CR_SKIM=1), into /data/uboone/processed/mp2pi_cr_skim.
# CC1mu2piBDT reports no control region of its own, so the skim writes only events in the control regions of the
# single-pion selection and drops every event in the signal region of either selection. Guards, on every file: the
# XSEC_CR_SKIM marker is present, no written event is selected by either selection, every written event is in a
# single-pion control region, the two-pion branches exist; otherwise the output is DELETED.
# Comparison: scripts/cr_mp2pi_check.py.
#   mkdir -p /data/uboone/processed/mp2pi_cr_skim/logs
#   sbatch --array=0-33 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/pid_diag_skim_manifest.list \
#          slurm/slurm_mp2pi_cr_skim.sbatch
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/mp2pi_cr_skim; MANIFEST="${MANIFEST:?}"
SELS=CC1mu1piXp,CC1mu2piBDTEvt
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"; F="$OUT/xsec-ana-$TAG.root"
export XSEC_CR_SKIM=1                      # never unset: beam-on files are among the inputs
ProcessNTuples "$RAW" "$FT" "$SELS" "$F" > "$OUT/$TAG.log" 2>&1; RC=$?
CHK=$(root.exe -l -b -q -e '
  auto f=TFile::Open("'"$F"'");
  if(!f||f->IsZombie()){ printf("BAD no-file\n"); return; }
  bool mark = f->Get("XSEC_CR_SKIM")!=nullptr;
  auto t=(TTree*)f->Get("stv_tree");
  if(!t){ printf("BAD no-tree\n"); return; }
  Long64_t sel = t->GetEntries("CC1mu1piXp_Selected||CC1mu2piBDTEvt_Selected||CC1mu2piBDTEvt_evt_region==0");
  Long64_t cr  = t->GetEntries("CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0||CC1mu1piXp_sb_cosmic");
  bool evt = t->GetBranch("CC1mu2piBDTEvt_evt_score")!=nullptr && t->GetBranch("CC1mu2piBDTEvt_etrk_pid_pi")!=nullptr;
  printf("%s marker=%d selected=%lld outside_cr=%lld evt=%d entries=%lld\n", (mark&&sel==0&&cr==t->GetEntries()&&evt)?"OK":"BAD",
         mark, sel, t->GetEntries()-cr, evt, t->GetEntries());
' 2>/dev/null | grep -E '^(OK|BAD)')
NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "GUARD $TAG :: $CHK write_errors=$NERR"
case "$CHK" in
  OK*) [ "$RC" = 0 ] && [ "$NERR" = 0 ] && { echo "DONE $TAG rc=0 guards=pass"; exit 0; }
       echo "FAIL $TAG rc=$RC write_errors=$NERR (guards passed)"; exit 1 ;;
  *)   echo "FAIL $TAG rc=$RC guards=FAILED -- deleting output"; rm -f "$F"; exit 2 ;;
esac
