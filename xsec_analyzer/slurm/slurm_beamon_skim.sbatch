#!/bin/bash
# SLURM array task: process ONE real beam-on file into the signal-region-stripped
# control-region skim (ProcessNTuples with XSEC_CR_SKIM=1). This is the ONLY way a
# beam-on file is ever processed before the signal region is approved: an event is
# written only if it FAILS the signal-region selection and belongs to at least one
# control region, and the file carries the XSEC_CR_SKIM marker that UniverseMaker
# refuses. The job verifies both after writing and DELETES the output if either
# fails, so a full beam-on processed file cannot be produced by accident.
#SBATCH --job-name=bonskim
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/beamon_skim; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"; F="$OUT/xsec-ana-$TAG.root"
export XSEC_CR_SKIM=1                      # never unset: this is a real beam-on file
ProcessNTuples "$RAW" onBNB CC1mu1piXp "$F"; RC=$?
# --- guards: marker present AND no written event carries the signal-region flag ---
CHK=$(root.exe -l -b -q -e '
  auto f=TFile::Open("'"$F"'");
  if(!f||f->IsZombie()){ printf("BAD no-file\n"); return; }
  bool mark = f->Get("XSEC_CR_SKIM")!=nullptr;
  auto t=(TTree*)f->Get("stv_tree");
  if(!t){ printf("BAD no-tree\n"); return; }
  Long64_t sel = t->GetEntries("CC1mu1piXp_Selected");
  Long64_t cr  = t->GetEntries("CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0||CC1mu1piXp_sb_cosmic");
  printf("%s marker=%d selected=%lld outside_cr=%lld\n", (mark&&sel==0&&cr==t->GetEntries())?"OK":"BAD", mark, sel, t->GetEntries()-cr);
' 2>/dev/null | grep -E '^(OK|BAD)')
echo "GUARD $TAG :: $CHK"
case "$CHK" in
  OK*) [ "$RC" = 0 ] && { echo "DONE $TAG rc=0 guards=pass"; exit 0; }
       echo "FAIL $TAG rc=$RC (guards passed but processing returned non-zero)"; exit 1 ;;
  *)   echo "FAIL $TAG rc=$RC guards=FAILED -- deleting output"; rm -f "$F"; exit 2 ;;
esac
