#!/bin/bash
# SLURM array: Phase 2 of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md. Reprocess ONE raw input of the
# held-out run periods of the particle classifier (FHC Run 5, RHC Run 3: overlays, beam-off, dirt)
# with the released selections and their particle-classifier variants in one pass, into
# /data/uboone/processed/mp_pid_eval, for an event-level comparison. No beam-on file is touched.
#   sbatch --array=0-9 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_pid_eval_manifest.list \
#          slurm/slurm_mp_pid_eval.sbatch
#SBATCH --job-name=mppideval
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/mp_pid_eval; MANIFEST="${MANIFEST:?}"
SELS=CC1mu1piXp,CC1mu1piXpNewPID,CC1mu2pi,CC1mu2piNewPID,CC1mu3pi,CC1mu3piNewPID
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
OUTF="$OUT/xsec-ana-$TAG.root"
ProcessNTuples "$RAW" "$FT" "$SELS" "$OUTF" > "$OUT/$TAG.log" 2>&1; RC=$?
NIN=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$RAW\");auto t=(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\");printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
CHK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%lld %d\n\",t?t->GetEntries():-1LL,(t&&t->GetBranch(\"CC1mu1piXpNewPID_Selected\")&&t->GetBranch(\"CC1mu2piNewPID_cutflow_bits\"))?1:0);" 2>/dev/null | tail -1)
NOUT=${CHK% *}; BR=${CHK#* }; NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "DONE $TAG rc=$RC in=$NIN out=$NOUT newpid_branches=$BR write_errors=$NERR"
[ "$RC" = 0 ] && [ "$BR" = 1 ] && [ "$NIN" = "$NOUT" ] && [ "$NERR" = 0 ]
