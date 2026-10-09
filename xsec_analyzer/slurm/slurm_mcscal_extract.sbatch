#!/bin/bash
#SBATCH --job-name=mcscalx
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=3500M
#SBATCH --time=10:00:00
#SBATCH --output=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
# SLURM array task: univmake + UnfolderNuMI for ONE MCS-calibrated result (line SLURM_ARRAY_TASK_ID of
# slurm/mcscal_manifest.list, "xsec bin slice"), as slurm/slurm_extract.sbatch but keeping the unfolded tables:
# /data/uboone/processed/mcscal/nominal_out/<tag>/{unfolded_signal.txt,cov_*.txt,run.log,xsec_<tag>.root}.
#   sbatch --array=0-18 --export=ALL,REPO=$PWD slurm/slurm_mcscal_extract.sbatch
set -u
REPO="${REPO:?}"; CAL=/data/uboone/processed/mcscal
read -r XC BIN SC <<< "$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$REPO/slurm/mcscal_manifest.list")"
[ -n "${SC:-}" ] || { echo "ERROR: no work on line $(( SLURM_ARRAY_TASK_ID + 1 ))" >&2; exit 1; }
TAG=$(basename "$XC" .txt)
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
export UNIV_THREADS="${SLURM_CPUS_PER_TASK:-4}"
U=$(awk '$1=="UnivFile"{print $2}' "$XC"); FP=$(awk '$1=="FPFile"{print $2}' "$XC")
for f in "$XC" "$BIN" "$SC" "$FP"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
awk '$1!~/^#/ && NF>=3 {print $1}' "$FP" | while read -r p; do [ -s "$p" ] || echo "MISSING $p"; done | grep -q . && { echo "ERROR: inputs missing" >&2; exit 1; }
echo "host=$(hostname) tag=$TAG univ=$U"
mkdir -p "$(dirname "$U")"; command rm -f -- "$U"
bin/univmake "$FP" "$BIN" "$U" "$FP" > "$(dirname "$U")/univ_${TAG}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $TAG" >&2; exit 1; }
D="$CAL/nominal_out/$TAG"; W="$CAL/nominal_out/work_$TAG"; command rm -rf -- "$W" "$D"; mkdir -p "$W/unfold_output" "$D"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$REPO/$XC" "$REPO/$SC" "$D/xsec_${TAG}.root" ) > "$D/run.log" 2>&1
rc=$?
for f in "$W"/unfold_output/mat_table_cov_*.txt; do [ -f "$f" ] && cp -p "$f" "$D/$(basename "$f" | sed 's/^mat_table_//')"; done
[ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$D/unfolded_signal.txt"
command rm -rf -- "$W"
[ $rc -eq 0 ] || { echo "FAILED unfold $TAG rc=$rc" >&2; exit 1; }
echo "OK $TAG cov=$(ls "$D"/cov_*.txt | wc -l) sigma_int=$(grep 'SYSTDUMP] sigma_int' "$D/run.log" | tail -1 | awk '{print $3}')"
