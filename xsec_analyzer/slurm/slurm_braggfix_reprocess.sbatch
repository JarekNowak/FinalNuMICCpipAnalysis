#!/bin/bash
# SLURM array: reprocess ONE unit of the release (overlay, dirt, detector variation or beam-off) with
# the selection without the Bragg-pion cut (commit 1fb6be9) into a STAGING directory
# (processed/incl_braggfix, processed/w_braggfix, processed/ext_perrun_braggfix); the live trees are
# not touched. Manifest slurm/braggfix_manifest.list: raw|filetype|selections|output. detVar
# part1+part2 inputs are merged first with an entry-count check. No beam-on file is processed.
#   sbatch --chdir=$PWD/slurm --array=0-84 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/braggfix_manifest.list \
#          slurm/slurm_braggfix_reprocess.sbatch
#SBATCH --job-name=braggfix
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT SEL OUTF <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${OUTF:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
unset XSEC_CR_SKIM CC1MU1PIXP_USE_PION_BDT
OUT=$(dirname "$OUTF"); TAG=$(basename "$OUTF" .root); mkdir -p "$OUT"
NENT(){ root.exe -l -b -q -e "Long64_t n=0;for(auto p:{$1}){auto f=TFile::Open(p);auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;if(t)n+=t->GetEntries();}printf(\"%lld\n\",n);" 2>/dev/null | tail -1; }
IN="$RAW"
if [[ "$RAW" == *"+"* ]]; then
  P1="${RAW%%+*}"; P2="${RAW##*+}"; IN="$OUT/merge_tmp_${TAG}.root"; rm -f "$IN"
  hadd -f "$IN" "$P1" "$P2" > "$OUT/${TAG}_hadd.log" 2>&1 || { echo "hadd failed"; rm -f "$IN"; exit 3; }
  NEXP=$(NENT "\"$P1\",\"$P2\""); NGET=$(NENT "\"$IN\"")
  [ "$NEXP" = "$NGET" ] || { echo "merge entry mismatch $NGET vs $NEXP"; rm -f "$IN"; exit 4; }
fi
NIN=$(NENT "\"$IN\"")
ProcessNTuples "$IN" "$FT" "$SEL" "$OUTF" > "$OUT/$TAG.log" 2>&1; RC=$?
[[ "$IN" == "$OUT/merge_tmp_"* ]] && rm -f "$IN"
FIRST=${SEL%%,*}
CHK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%lld %d\n\",t?t->GetEntries():-1LL,(t&&t->GetBranch(\"${FIRST}_Selected\"))?1:0);" 2>/dev/null | tail -1)
NOUT=${CHK% *}; BR=${CHK#* }; NERR=$(grep -c 'WriteTObject' "$OUT/$TAG.log")
echo "DONE $TAG rc=$RC in=$NIN out=$NOUT branch=$BR write_errors=$NERR"
[ "$RC" = 0 ] && [ "$BR" = 1 ] && [ "$NIN" = "$NOUT" ] && [ "$NERR" = 0 ]
