#!/bin/bash
# SLURM array: Phase 2, items 1-3, of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md. Reprocess ONE raw
# input of the release list (the 28 overlays, beam-off and dirt of the Phase 0 manifest) with the
# event-level track dump (CC1mu2piEvt) next to the current two- and three-pion selections, their
# particle-classifier variants and the released single-pion selection, into
# /data/uboone/processed/mp_evt, for the track assignment and the event classifier
# (scripts/mp_evt_*.py). No beam-on file is touched.
#   sbatch --array=0-27 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_phase0_manifest.list \
#          slurm/slurm_mp_evt.sbatch
#SBATCH --job-name=mpevt
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/mp_evt; MANIFEST="${MANIFEST:?}"
SELS=CC1mu2piEvt,CC1mu2piNewPID,CC1mu3pi,CC1mu3piNewPID,CC1mu1piXp
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
OUTF="$OUT/xsec-ana-$TAG.root"
ProcessNTuples "$RAW" "$FT" "$SELS" "$OUTF" > "$OUT/$TAG.log" 2>&1; RC=$?
# checks: the dump branches exist, every input event was written, no basket write failed
NIN=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$RAW\");auto t=(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\");printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
CHK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%lld %d\n\",t?t->GetEntries():-1LL,(t&&t->GetBranch(\"CC1mu2piEvt_etrk_pid_pi\")&&t->GetBranch(\"CC1mu3pi_MC_Signal\"))?1:0);" 2>/dev/null | tail -1)
NOUT=${CHK% *}; BR=${CHK#* }; NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "DONE $TAG rc=$RC in=$NIN out=$NOUT dump_branches=$BR write_errors=$NERR"
[ "$RC" = 0 ] && [ "$BR" = 1 ] && [ "$NIN" = "$NOUT" ] && [ "$NERR" = 0 ]
