#!/bin/bash
#SBATCH --job-name=mcsmap
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=2300M
#SBATCH --time=03:00:00
#SBATCH --output=/data/uboone/processed/mcsmap_study/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcsmap_study/logs/%x_%A_%a.out
# SLURM array task: impact of the MCS mis-modelling of exiting muons on ONE released result that uses the muon momentum
# (line SLURM_ARRAY_TASK_ID of slurm/mcsmap_products.list: tag|xsec config|slice config|bin config|family).
# The release universe file is copied twice. "base": the release extraction on the first copy (systematics as in
# version 1.8, DVfrac). "map": the pseudo-data with their MCS-measured muon momenta mapped to the control-region data
# (scripts/mcs_map_derive.py, macros/mcs_map_fakedata{,_1p}.C): univmake runs on the mapped pseudo-data files alone,
# their per-file directories replace the nominal ones in the second copy (macros/merge_univ_dirs.C), and the release
# extraction runs with the nominal response. Everything is written under /data/uboone/processed/mcsmap_study.
#   mkdir -p /data/uboone/processed/mcsmap_study/logs
#   sbatch --array=0-17 --export=ALL,REPO=$PWD slurm/slurm_mcsmap_impact.sbatch
set -u
REPO="${REPO:?}"; OUT=/data/uboone/processed/mcsmap_study; PROC=/data/uboone/processed
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r TAG XC SC BIN FAM <<< "$(sed -n "${line}p" "$REPO/slurm/mcsmap_products.list")"
[ -n "${FAM:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
mkdir -p "$OUT/univ" "$OUT/cfg" "$OUT/base" "$OUT/map" "$OUT/work"
echo "host=$(hostname) tag=$TAG"
U=$(awk '$1=="UnivFile"{print $2}' "$XC"); FP=$(awk '$1=="FPFile"{print $2}' "$XC")
UA="$OUT/univ/${TAG}_base_univmake.root"; UB="$OUT/univ/${TAG}_map_univmake.root"
cp -p "$U" "$UA.tmp" && mv "$UA.tmp" "$UA" && cp -p "$U" "$UB.tmp" && mv "$UB.tmp" "$UB" || { echo "FAIL copy"; exit 1; }
if [ "$FAM" = incl ]; then FROM="$PROC/xsec-ana-fakedata_"; TO="$OUT/fakedata/xsec-ana-fakedata_"
else FROM="$PROC/w/xsec-ana-fakedata_"; TO="$OUT/fakedata/w/xsec-ana-fakedata_"; fi
FPM="$OUT/cfg/${TAG}_fp_map.txt"; sed "s|^$FROM|$TO|" "$FP" > "$FPM"
NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FP" | wc -l); NMAP=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPM" | grep -c "^$TO")
[ "$NON" -gt 0 ] && [ "$NON" = "$NMAP" ] || { echo "FAIL file list: $NON onBNB, $NMAP mapped"; exit 1; }
LIST="$OUT/cfg/${TAG}_fakedata.list"; awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPM" > "$LIST"
SMALL="$OUT/work/${TAG}_fakedata_univmake.root"; rm -f "$SMALL"
UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$LIST" "$BIN" "$SMALL" "$FPM" > "$OUT/work/${TAG}_univmake.log" 2>&1
echo "univmake rc=$? (the closing total is expected to fail)"
root -l -b -q "macros/merge_univ_dirs.C(\"$SMALL\",\"$UB\")" > "$OUT/work/${TAG}_merge.log" 2>&1
grep -q "^$NON per-file directories copied under 1 root" "$OUT/work/${TAG}_merge.log" || { echo "FAIL merge"; cat "$OUT/work/${TAG}_merge.log"; exit 1; }
sed -e "s|^UnivFile .*|UnivFile $UA|" -e "s|^SystFile .*|SystFile configs/ccpi_systcalc_numi_dvfrac.conf|" "$XC" > "$OUT/cfg/${TAG}_base.txt"
sed -e "s|^UnivFile .*|UnivFile $UB|" -e "s|^SystFile .*|SystFile configs/ccpi_systcalc_numi_dvfrac.conf|" \
    -e "s|^FPFile .*|FPFile $FPM|" "$XC" > "$OUT/cfg/${TAG}_map.txt"
for mode in base map; do
  W="$OUT/work/${TAG}_$mode"; D="$OUT/$mode/$TAG"; rm -rf "$W" "$D"; mkdir -p "$W/unfold_output" "$D"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$OUT/cfg/${TAG}_$mode.txt" "$REPO/$SC" "$D/xsec_${TAG}.root" ) > "$D/run.log" 2>&1
  rc=$?
  for f in "$W"/unfold_output/mat_table_cov_*.txt; do [ -f "$f" ] && cp -p "$f" "$D/$(basename "$f" | sed 's/^mat_table_//')"; done
  [ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$D/unfolded_signal.txt"
  echo "$mode rc=$rc cov=$(ls "$D"/cov_*.txt 2>/dev/null | wc -l) cache=$(grep -c 'Reusing cached universes' "$D/run.log") sigma_int=$(grep 'SYSTDUMP] sigma_int' "$D/run.log" | awk '{print $3}')"
  rm -rf "$W"
done
command rm -f "$UA" "$UB"
echo "DONE $TAG"
