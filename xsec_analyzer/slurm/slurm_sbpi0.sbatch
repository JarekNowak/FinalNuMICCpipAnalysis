#!/bin/bash
# SLURM array: reprocess ONE raw input with the sideband-instrumented ProcessNTuples
# (sb_pi0_final / sb_nprimtrk / sb_nnonproton branches, 2026-09-06) into a SEPARATE tree
# /data/uboone/processed/sb_pi0 for the pre-ratification comparison of the two pi0-region
# definitions. Monte Carlo, dirt, detVar and beam-off only: no beam-on file is touched.
#SBATCH --job-name=sbpi0
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/sb_pi0; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
IN="$RAW"
if [[ "$RAW" == *"+"* ]]; then   # detVar part1+part2: merge, verify the entry count, never modify the parts
  P1="${RAW%%+*}"; P2="${RAW##*+}"; IN="$OUT/merge_tmp_${TAG}.root"; rm -f "$IN"
  hadd -f "$IN" "$P1" "$P2" > "$OUT/merge_${TAG}.log" 2>&1 || { echo "hadd failed"; exit 3; }
  NEXP=$(root.exe -l -b -q -e "Long64_t n=0;for(auto p:{\"$P1\",\"$P2\"}){auto f=TFile::Open(p);auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;if(t)n+=t->GetEntries();}printf(\"%lld\n\",n);" 2>/dev/null | tail -1)
  NGET=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$IN\");auto t=(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\");printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
  [ "$NEXP" = "$NGET" ] || { echo "merge entry mismatch $NEXP vs $NGET"; exit 4; }
fi
ProcessNTuples "$IN" "$FT" CC1mu1piXp "$OUT/xsec-ana-$TAG.root"; RC=$?
[[ "$IN" == "$OUT/merge_tmp_"* ]] && rm -f "$IN"
OK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUT/xsec-ana-$TAG.root\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%d\n\",(t&&t->GetBranch(\"CC1mu1piXp_sb_pi0_final\"))?1:0);" 2>/dev/null | tail -1)
echo "DONE $TAG rc=$RC sb_pi0_final=$OK"; [ "$RC" = 0 ] && [ "$OK" = 1 ]
