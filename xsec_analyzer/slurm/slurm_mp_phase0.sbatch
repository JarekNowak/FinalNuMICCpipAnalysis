#!/bin/bash
# SLURM array: Phase 0 of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md. Reprocess ONE raw input
# with the inclusive, two-pion and three-pion selections in a single pass
# (CC1mu1piXp,CC1mu2pi,CC1mu3pi), so per-event overlaps between the three are available,
# into a SEPARATE tree /data/uboone/processed/mp_phase0. Monte Carlo, dirt and beam-off
# only: no beam-on file is touched (the 1pi and 2pi signal regions are blind).
#   sbatch --array=0-27 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_phase0_manifest.list \
#          slurm/slurm_mp_phase0.sbatch
#SBATCH --job-name=mpph0
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/mp_phase0; MANIFEST="${MANIFEST:?}"
SELS=CC1mu1piXp,CC1mu2pi,CC1mu3pi
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$OUT"
case "$FT" in onBNB) echo "REFUSED: beam-on input $RAW (blind analysis)"; exit 2;; esac
OUTF="$OUT/xsec-ana-$TAG.root"
ProcessNTuples "$RAW" "$FT" "$SELS" "$OUTF"; RC=$?
# checks: the new multi-pion branches exist and every input event was written
NIN=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$RAW\");auto t=(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\");printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
CHK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%lld %d\n\",t?t->GetEntries():-1LL,(t&&t->GetBranch(\"CC1mu2pi_cutflow_bits\")&&t->GetBranch(\"CC1mu3pi_cutflow_bits\"))?1:0);" 2>/dev/null | tail -1)
NOUT=${CHK% *}; BR=${CHK#* }
echo "DONE $TAG rc=$RC in=$NIN out=$NOUT mp_branches=$BR"
[ "$RC" = 0 ] && [ "$BR" = 1 ] && [ "$NIN" = "$NOUT" ]
