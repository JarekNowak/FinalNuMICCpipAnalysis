#!/bin/bash
#SBATCH --job-name=mcscalv
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=2500M
#SBATCH --time=03:00:00
#SBATCH --output=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
# SLURM array task: calibration uncertainty of ONE MCS-calibrated result for ONE variant map. Task i: product
# i / 7 of slurm/mcscal_manifest.list, variant i % 7 of VARIANTS. The calibrated universe file (slurm/slurm_extract.sbatch
# on slurm/mcscal_manifest.list) is copied; univmake runs on the pseudo-data calibrated with the variant map alone, their
# per-file directories are added to the copy (macros/merge_univ_dirs.C), and the extraction runs with the pseudo-data of
# the variant against the response calibrated with the nominal map. Output: /data/uboone/processed/mcscal/var/<variant>/<tag>/.
#   sbatch --array=0-132 --export=ALL,REPO=$PWD slurm/slurm_mcscal_var.sbatch
set -u
REPO="${REPO:?}"; CAL=/data/uboone/processed/mcscal
VARIANTS=(shrink coarse fine old_production new_production FHC RHC)
P=$(( SLURM_ARRAY_TASK_ID / ${#VARIANTS[@]} )); V=${VARIANTS[$(( SLURM_ARRAY_TASK_ID % ${#VARIANTS[@]} ))]}
read -r XC BIN SC <<< "$(sed -n "$(( P + 1 ))p" "$REPO/slurm/mcscal_manifest.list")"
[ -n "${SC:-}" ] || { echo "ERROR: no product $P" >&2; exit 1; }
TAG=$(basename "$XC" .txt)
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
echo "host=$(hostname) tag=$TAG variant=$V"
U=$(awk '$1=="UnivFile"{print $2}' "$XC"); FP=$(awk '$1=="FPFile"{print $2}' "$XC")
[ -s "$U" ] || { echo "FAIL no calibrated universes $U" >&2; exit 1; }
D="$CAL/var/$V/$TAG"; WK="$CAL/var/work"; mkdir -p "$D" "$WK"
UV="$WK/${TAG}_${V}_univmake.root"
cp -p "$U" "$UV.tmp" && mv "$UV.tmp" "$UV" || { echo "FAIL copy" >&2; exit 1; }
FPV="$WK/${TAG}_${V}_fp.txt"; sed "s|^$CAL/nominal/\(\(w/\)\{0,1\}xsec-ana-fakedata_\)|$CAL/$V/\1|" "$FP" > "$FPV"
NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FP" | wc -l); NV=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPV" | grep -c "^$CAL/$V/")
[ "$NON" -gt 0 ] && [ "$NON" = "$NV" ] || { echo "FAIL file list: $NON onBNB, $NV variant" >&2; exit 1; }
LIST="$WK/${TAG}_${V}_fakedata.list"; awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPV" > "$LIST"
SMALL="$WK/${TAG}_${V}_small_univmake.root"; command rm -f -- "$SMALL"
UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$LIST" "$BIN" "$SMALL" "$FPV" > "$WK/${TAG}_${V}_univmake.log" 2>&1
echo "univmake rc=$? (the closing total is expected to fail)"
root -l -b -q "macros/merge_univ_dirs.C(\"$SMALL\",\"$UV\")" > "$WK/${TAG}_${V}_merge.log" 2>&1
grep -q "^$NON per-file directories copied under 1 root" "$WK/${TAG}_${V}_merge.log" || { echo "FAIL merge" >&2; cat "$WK/${TAG}_${V}_merge.log"; exit 1; }
XV="$WK/${TAG}_${V}_xsec.txt"
sed -e "s|^UnivFile .*|UnivFile $UV|" -e "s|^FPFile .*|FPFile $FPV|" "$XC" > "$XV"
W="$WK/${TAG}_${V}_unfold"; command rm -rf -- "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XV" "$REPO/$SC" "$D/xsec_${TAG}.root" ) > "$D/run.log" 2>&1
rc=$?
for f in "$W"/unfold_output/mat_table_cov_*.txt; do [ -f "$f" ] && cp -p "$f" "$D/$(basename "$f" | sed 's/^mat_table_//')"; done
[ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$D/unfolded_signal.txt"
echo "rc=$rc cov=$(ls "$D"/cov_*.txt 2>/dev/null | wc -l) cache=$(grep -c 'Reusing cached universes' "$D/run.log") sigma_int=$(grep 'SYSTDUMP] sigma_int' "$D/run.log" | awk '{print $3}')"
command rm -rf -- "$W" "$UV" "$SMALL"
echo "DONE $TAG $V"
