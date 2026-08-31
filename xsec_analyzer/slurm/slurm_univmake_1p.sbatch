#!/bin/bash
# SLURM array task: rebuild ONE proton-tagged univmake with the detector variations
# included. One manifest line per task: "<cfg> <obs>".
#
# univmake is internally multi-threaded and uses ~13 of 16 cores, so each task takes a
# whole node rather than sharing one.
#SBATCH --job-name=univ1p
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
MANIFEST="${MANIFEST:?MANIFEST not set (submit via submit_univmake_1p.sh)}"
REPO="${REPO:?REPO not set}"
PROC="${PROC:-/data/uboone/processed}"

line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r CFG OBS <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${CFG:-}" ] && [ -n "${OBS:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }

cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"

declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
OUT="$PROC/ccpi1p_${TAG[$CFG]}_${OBS}_univmake.root"
FPM="configs/file_properties_numi_${CFG}_w.txt"
# Bin-config naming: a coarsened observable "<base>2bin" lives in
# ccpi1p_<base>_bin_config_2bin.txt, not ccpi1p_<base>2bin_bin_config.txt.
case "$OBS" in
  *2bin) BIN="configs/ccpi1p_${OBS%2bin}_bin_config_2bin.txt" ;;
  *)     BIN="configs/ccpi1p_${OBS}_bin_config.txt" ;;
esac

for f in "$FPM" "$BIN"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

# Rebuild if the existing output predates the file_properties file. The whole point of
# this pass is that the config gained detVar entries, so a plain skip-if-exists would
# skip every task and silently do nothing.
if [ -f "$OUT" ] && [ "$OUT" -nt "$FPM" ]; then
  sz=$(stat -c%s "$OUT")
  if [ "$sz" -gt 30000000 ]; then
    echo "[skip] $CFG $OBS -- output newer than $FPM ($((sz/1000000))MB)"
    exit 0
  fi
fi

echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID cfg=$CFG obs=$OBS"
echo "  FPM=$FPM"
echo "  BIN=$BIN"
echo "  OUT=$OUT"
echo "  detVar entries in FPM: $(grep -vE '^[[:space:]]*#' "$FPM" | awk '$3 ~ /^detVar/' | wc -l)"

# Write to a temporary name and move into place, so an interrupted task cannot leave a
# truncated univmake that a later skip-if-exists would accept as complete.
TMP="${OUT}.partial.$$"
t0=$(date +%s)
"$REPO/bin/univmake" "$FPM" "$BIN" "$TMP" "$FPM"
rc=$?
t1=$(date +%s)

if [ $rc -ne 0 ] || [ ! -s "$TMP" ]; then
  echo "FAILED $CFG $OBS rc=$rc after $(( (t1-t0)/60 ))min" >&2
  rm -f "$TMP"
  exit 1
fi
sz=$(stat -c%s "$TMP")
if [ "$sz" -lt 30000000 ]; then
  echo "FAILED $CFG $OBS -- output only $((sz/1000000))MB, expected >30MB" >&2
  rm -f "$TMP"
  exit 1
fi
mv -f "$TMP" "$OUT"
echo "OK $CFG $OBS $((sz/1000000))MB in $(( (t1-t0)/60 ))min"
