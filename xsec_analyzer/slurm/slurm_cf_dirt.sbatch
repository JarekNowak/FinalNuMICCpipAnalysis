#!/bin/bash
# SLURM array: cut-flow inputs (h_cutflow_* counters, full CV weight) for the per-run dirt samples,
# one selection per pass so the counters are not overwritten: CC1mu1piXp -> processed/cf,
# CC1mu1pi1p -> processed/cf_1p. Dirt only. Manifest lines RAW|FT|TAG|OUTDIR|SELECTION.
#SBATCH --job-name=cfdirt
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW FT TAG OUT SEL <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${SEL:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
[ "$FT" = dirtMC ] || { echo "REFUSED: not a dirt input"; exit 2; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
OUTF="$OUT/xsec-ana-$TAG.root"
[ -e "$OUTF" ] && { echo "REFUSED: $OUTF exists"; exit 3; }
ProcessNTuples "$RAW" "$FT" "$SEL" "$OUTF"; RC=$?
OK=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");printf(\"%d\n\",(f&&!f->IsZombie()&&f->Get(\"h_cutflow_tot\"))?1:0);" 2>/dev/null | tail -1)
echo "DONE $TAG $SEL rc=$RC cutflow_hist=$OK"; [ "$RC" = 0 ] && [ "$OK" = 1 ]
