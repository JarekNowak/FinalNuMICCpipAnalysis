#!/bin/bash
# SLURM array: the two-pion production for the framework (Phase 4, item 3, of
# report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md). Reprocess ONE raw input of the release list (the 28 overlays,
# beam-off and dirt of the Phase 0 manifest) with the two-pion selection with the event classifier (CC1mu2piBDT) and
# the released single-pion selection (CC1mu1piXp, whose signal region the two-pion sidebands leave out), into
# /data/uboone/processed/mp2pi. No beam-on file is touched.
#   sbatch --array=0-27 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_phase0_manifest.list slurm/slurm_mp2pi.sbatch
#SBATCH --job-name=mp2pi
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT="${OUT:-/data/uboone/processed/mp2pi}"; MANIFEST="${MANIFEST:?}"
SELS=CC1mu2piBDT,CC1mu1piXp
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
unset CC1MU1PIXP_MP_PID_MODEL XSEC_CR_SKIM CC1MU1PIXP_USE_PION_BDT
mkdir -p "$OUT"
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
OUTF="$OUT/xsec-ana-$TAG.root"
ProcessNTuples "$RAW" "$FT" "$SELS" "$OUTF" > "$OUT/$TAG.log" 2>&1; RC=$?
NIN=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$RAW\");auto t=(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\");printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
CHK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%lld %d\n\",t?t->GetEntries():-1LL,(t&&t->GetBranch(\"CC1mu2piBDT_evt_region\")&&t->GetBranch(\"CC1mu1piXp_Selected\"))?1:0);" 2>/dev/null | tail -1)
NOUT=${CHK% *}; BR=${CHK#* }; NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "DONE $TAG rc=$RC in=$NIN out=$NOUT branches=$BR write_errors=$NERR"
[ "$RC" = 0 ] && [ "$BR" = 1 ] && [ "$NIN" = "$NOUT" ] && [ "$NERR" = 0 ]
