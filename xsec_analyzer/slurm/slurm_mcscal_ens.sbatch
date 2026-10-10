#!/bin/bash
#SBATCH --job-name=mcscale
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=2500M
#SBATCH --time=02:00:00
#SBATCH --output=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
# SLURM array task: one member of a release ensemble of a calibrated result (version 1.9). Task i: ensemble i / 100 of
# slurm/mcscal_ens_manifest.list ("member name, calibrated tag, release config, bin, slice"), member t = i % 100 + 1.
# The member's throw (calibrated copy, scripts/mcs_cal_setup.py validation) is the only data difference from the release,
# so the calibrated universes are copied, univmake runs on the throw alone and its per-file directories are added
# (macros/merge_univ_dirs.C). The member files keep the names of the release ensembles in /data/uboone/processed/ens
# (read by report/tools/ensemble_stat_pulls.py); those of version 1.8 are first moved to ens/backup_v1.8.
#   sbatch --array=0-499%40 --export=ALL,REPO=$PWD slurm/slurm_mcscal_ens.sbatch
set -u
REPO="${REPO:?}"; CAL=/data/uboone/processed/mcscal; E=/data/uboone/processed/ens
K=$(( SLURM_ARRAY_TASK_ID / 100 )); T=$(( SLURM_ARRAY_TASK_ID % 100 + 1 ))
read -r N TAG XC0 BIN SC <<< "$(sed -n "$(( K + 1 ))p" "$REPO/slurm/mcscal_ens_manifest.list")"
[ -n "${SC:-}" ] || { echo "ERROR: no ensemble $K" >&2; exit 1; }
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
echo "host=$(hostname) member=${N}_t$T tag=$TAG"
U0="$CAL/univ/${TAG}_univmake.root"; [ -s "$U0" ] || { echo "FAIL no calibrated universes $U0" >&2; exit 1; }
FPC=$(awk '$1=="FPFile"{print $2}' "$XC0")
B="$E/backup_v1.8"; mkdir -p "$B" "$E/work"
for f in "$E/univ_${N}_t${T}.root" "$E/univ_${N}_t${T}.log" "$E/xsec_${N}_t${T}.txt" "$E/fp_${N}_t${T}.txt" "$E/statcov_${N}_t${T}" \
         "$E/closure_hists_xsec_${N}_t${T}.root" "$E/closure_hists_all_xsec_${N}_t${T}.root" "$E/unfold_${N}_t${T}.log" "$E/xsec_${N}_t${T}.root"; do
  [ -e "$f" ] && [ ! -e "$B/$(basename "$f")" ] && mv "$f" "$B/"
done
FPM="$E/fp_${N}_t${T}.txt"
sed -E "s#^$CAL/nominal/(w/)?xsec-ana-fakedata_(fhc|rhc)_run([0-9])\.root#$CAL/ens/\1throws/fakedata_\2_run\3_t${T}.root#" "$FPC" > "$FPM.tmp" && mv "$FPM.tmp" "$FPM"
NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPM" | wc -l); NT=$(grep -c "_t${T}\.root" "$FPM")
[ "$NON" -gt 0 ] && [ "$NON" = "$NT" ] || { echo "FAIL fp rewrite: $NON onBNB, $NT throws" >&2; exit 1; }
LIST="$E/work/${N}_t${T}.list"; awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPM" > "$LIST"
while read -r f; do [ -s "$f" ] || { echo "FAIL missing $f" >&2; exit 1; }; done < "$LIST"
U="$E/univ_${N}_t${T}.root"; cp -p "$U0" "$U.tmp" && mv "$U.tmp" "$U" || { echo "FAIL copy" >&2; exit 1; }
SMALL="$E/work/small_${N}_t${T}.root"; command rm -f -- "$SMALL"
UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$LIST" "$BIN" "$SMALL" "$FPM" > "$E/univ_${N}_t${T}.log" 2>&1
root -l -b -q "macros/merge_univ_dirs.C(\"$SMALL\",\"$U\")" > "$E/work/merge_${N}_t${T}.log" 2>&1
grep -q "^$NON per-file directories copied under 1 root" "$E/work/merge_${N}_t${T}.log" || { echo "FAIL merge" >&2; cat "$E/work/merge_${N}_t${T}.log"; exit 1; }
XC="$E/xsec_${N}_t${T}.txt"; sed -E "s#^UnivFile .*#UnivFile $U#; s#^FPFile .*#FPFile $FPM#" "$XC0" > "$XC"
W="$E/statcov_${N}_t${T}"; command rm -rf -- "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SC" "$E/xsec_${N}_t${T}.root" ) > "$W/run.log" 2>&1
rc=$?
cp -p "$W/run.log" "$E/unfold_${N}_t${T}.log"
command rm -f -- "$SMALL" "$LIST" "$E/work/merge_${N}_t${T}.log"
[ $rc -eq 0 ] && [ -f "$W/unfold_output/mat_table_cov_DataStats.txt" ] || { echo "FAIL unfold rc=$rc" >&2; exit 1; }
echo "DONE ${N}_t$T $(grep -m1 '^truth: ' "$W/run.log")"
