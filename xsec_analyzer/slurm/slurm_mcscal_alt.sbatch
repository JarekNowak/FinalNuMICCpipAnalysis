#!/bin/bash
#SBATCH --job-name=mcscala
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=2500M
#SBATCH --time=03:00:00
#SBATCH --output=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
# SLURM array task: one alternative-model extraction of a calibrated result (version 1.9), line SLURM_ARRAY_TASK_ID of
# slurm/mcscal_alt_manifest.list (scripts/mcs_cal_setup.py validation): "alt config, bin, slice, calibrated tag, output".
# The alternative-model file list differs from the nominal one only in its pseudo-data, so the calibrated nominal
# universes are copied, univmake runs on the calibrated alternative pseudo-data alone and their per-file directories are
# added (macros/merge_univ_dirs.C). The extraction writes the version-1.8 output name in rebuild_alt (read by
# report/tools/altmodel_eval.py); the version-1.8 sidecars and log are first moved to rebuild_alt/backup_v1.8.
#   sbatch --array=0-13 --export=ALL,REPO=$PWD slurm/slurm_mcscal_alt.sbatch
set -u
REPO="${REPO:?}"; CAL=/data/uboone/processed/mcscal
read -r XC BIN SC TAG OUTROOT <<< "$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$REPO/slurm/mcscal_alt_manifest.list")"
[ -n "${OUTROOT:-}" ] || { echo "ERROR: no work on line $(( SLURM_ARRAY_TASK_ID + 1 ))" >&2; exit 1; }
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
U0="$CAL/univ/${TAG}_univmake.root"; UA=$(awk '$1=="UnivFile"{print $2}' "$XC"); FPA=$(awk '$1=="FPFile"{print $2}' "$XC")
STEM=$(basename "$OUTROOT" .root); STEM=${STEM#xsec_}; D=$(dirname "$OUTROOT"); B="$D/backup_v1.8"; WK="$CAL/alt/work"
echo "host=$(hostname) xc=$XC tag=$TAG stem=$STEM"
[ -s "$U0" ] || { echo "FAIL no calibrated universes $U0" >&2; exit 1; }
mkdir -p "$(dirname "$UA")" "$WK" "$B"
cp -p "$U0" "$UA.tmp" && mv "$UA.tmp" "$UA" || { echo "FAIL copy" >&2; exit 1; }
NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPA" | wc -l)
LIST="$WK/${STEM}_fakedata.list"; awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPA" > "$LIST"
while read -r f; do [ -s "$f" ] || { echo "FAIL missing $f" >&2; exit 1; }; done < "$LIST"
SMALL="$WK/${STEM}_small_univmake.root"; command rm -f -- "$SMALL"
UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$LIST" "$BIN" "$SMALL" "$FPA" > "$WK/${STEM}_univmake.log" 2>&1
echo "univmake rc=$? (the closing total is expected to fail)"
root -l -b -q "macros/merge_univ_dirs.C(\"$SMALL\",\"$UA\")" > "$WK/${STEM}_merge.log" 2>&1
grep -q "^$NON per-file directories copied under 1 root" "$WK/${STEM}_merge.log" || { echo "FAIL merge" >&2; cat "$WK/${STEM}_merge.log"; exit 1; }
for f in "$D/closure_hists_xsec_${STEM}.root" "$D/closure_hists_all_xsec_${STEM}.root" "$D/unfold_${STEM}.log" "$OUTROOT"; do
  [ -f "$f" ] && [ ! -e "$B/$(basename "$f")" ] && mv "$f" "$B/"
done
W="$WK/${STEM}_unfold"; command rm -rf -- "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$REPO/$XC" "$REPO/$SC" "$OUTROOT" ) > "$D/unfold_${STEM}.log" 2>&1
rc=$?
command rm -rf -- "$W" "$SMALL"
[ $rc -eq 0 ] && [ -s "$D/closure_hists_xsec_${STEM}.root" ] || { echo "FAIL unfold rc=$rc" >&2; exit 1; }
echo "DONE $STEM cache=$(grep -c 'Reusing cached universes' "$D/unfold_${STEM}.log") $(grep -m1 '^truth: ' "$D/unfold_${STEM}.log")"
