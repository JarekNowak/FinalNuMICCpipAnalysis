#!/bin/bash
# SLURM array: PID diagnostic in the opened control regions (report/planning/multipion/PHASE2_SUMMARY.md).
# Processes ONE raw input (beam-on, overlay, beam-off or dirt) with the released single-pion
# selection and its PID-diagnostic twin CC1mu1piXpPIDDiag, ALWAYS as the signal-region-stripped
# control-region skim (XSEC_CR_SKIM=1), into /data/uboone/processed/pid_diag_skim. The guards of
# slurm_beamon_skim.sbatch apply to every file: the XSEC_CR_SKIM marker is present, no written
# event is selected by either selection, every written event is in a control region; otherwise the
# output is DELETED. Comparison: scripts/cr_pid_diag_check.py.
#   sbatch --array=0-33 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/pid_diag_skim_manifest.list \
#          slurm/slurm_pid_diag_skim.sbatch
#SBATCH --job-name=piddiag
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/pid_diag_skim; MANIFEST="${MANIFEST:?}"
SELS=CC1mu1piXp,CC1mu1piXpPIDDiag
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"; F="$OUT/xsec-ana-$TAG.root"
export XSEC_CR_SKIM=1                      # never unset: beam-on files are among the inputs
ProcessNTuples "$RAW" "$FT" "$SELS" "$F" > "$OUT/$TAG.log" 2>&1; RC=$?
# --- guards: marker present, no written event selected, every written event in a control region ---
CHK=$(root.exe -l -b -q -e '
  auto f=TFile::Open("'"$F"'");
  if(!f||f->IsZombie()){ printf("BAD no-file\n"); return; }
  bool mark = f->Get("XSEC_CR_SKIM")!=nullptr;
  auto t=(TTree*)f->Get("stv_tree");
  if(!t){ printf("BAD no-tree\n"); return; }
  Long64_t sel = t->GetEntries("CC1mu1piXp_Selected||CC1mu1piXpPIDDiag_Selected");
  Long64_t cr  = t->GetEntries("CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0||CC1mu1piXp_sb_cosmic");
  bool diag = t->GetBranch("CC1mu1piXpPIDDiag_pdiag_bragg_p")!=nullptr;
  printf("%s marker=%d selected=%lld outside_cr=%lld diag=%d\n", (mark&&sel==0&&cr==t->GetEntries()&&diag)?"OK":"BAD", mark, sel, t->GetEntries()-cr, diag);
' 2>/dev/null | grep -E '^(OK|BAD)')
NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "GUARD $TAG :: $CHK write_errors=$NERR"
case "$CHK" in
  OK*) [ "$RC" = 0 ] && [ "$NERR" = 0 ] && { echo "DONE $TAG rc=0 guards=pass"; exit 0; }
       echo "FAIL $TAG rc=$RC write_errors=$NERR (guards passed)"; exit 1 ;;
  *)   echo "FAIL $TAG rc=$RC guards=FAILED -- deleting output"; rm -f "$F"; exit 2 ;;
esac
