#!/bin/bash
# SLURM array task: unfold ONE proton-tagged extraction from the detVar-carrying
# universes, and harvest its covariance and A_C matrices for the data release.
#
# UnfolderNuMI writes its matrix dumps to a FIXED relative path (unfold_output/...), so
# concurrent unfolds sharing a working directory would overwrite each other's matrices --
# silently, since each run would still find a plausible file. Every task therefore gets
# its own scratch directory with the repo symlinked in.
#SBATCH --job-name=unf1p
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
MANIFEST="${MANIFEST:?MANIFEST not set}"
REPO="${REPO:?REPO not set}"
PROC="${PROC:-/data/uboone/processed}"
RELEASE="${RELEASE:?RELEASE not set}"
SCRATCH="${SCRATCH:-/data/uboone/processed/unfold_scratch}"

line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r CFG OBS <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${CFG:-}" ] && [ -n "${OBS:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }

declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
T="${TAG[$CFG]}"
UNIV="$PROC/ccpi1p_${T}_${OBS}_univmake.root"
[ -s "$UNIV" ] || { echo "ERROR: univmake missing: $UNIV" >&2; exit 1; }

# The whole point of this pass is the detector systematic, so refuse to unfold a
# universe file that predates the config which introduced it.
FPM="$REPO/configs/file_properties_numi_${CFG}_w.txt"
if [ "$FPM" -nt "$UNIV" ]; then
  echo "ERROR: $UNIV is older than $FPM -- it has no detVar universes" >&2
  exit 1
fi

case "$OBS" in
  *2bin) SLICE="configs/ccpi1p_${OBS%2bin}_slice_config_2bin.txt" ;;
  *)     SLICE="configs/ccpi1p_${OBS}_slice_config.txt" ;;
esac
XSEC="configs/ccpi1p_xsec_config_numi_${OBS}_${CFG}.txt"

WORK="$SCRATCH/${T}_${OBS}.$SLURM_JOB_ID"
rm -rf "$WORK"; mkdir -p "$WORK/unfold_output" || exit 1
ln -s "$REPO/configs" "$WORK/configs"
ln -s "$REPO/bin"     "$WORK/bin"
ln -s "$REPO/lib"     "$WORK/lib"
cd "$WORK" || exit 1

[ -f "$XSEC" ]  || { echo "ERROR: missing $XSEC"  >&2; exit 1; }
[ -f "$SLICE" ] || { echo "ERROR: missing $SLICE" >&2; exit 1; }

export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"

echo "host=$(hostname) cfg=$CFG obs=$OBS"
echo "  univ=$UNIV ($(date -r "$UNIV" '+%m-%d %H:%M'))"
echo "  xsec=$XSEC"
echo "  slice=$SLICE"

t0=$(date +%s)
"$REPO/bin/UnfolderNuMI" "$XSEC" "$SLICE" "$PROC/xsec_ccpi1p_${T}_${OBS}.root"
rc=$?
t1=$(date +%s)
[ $rc -eq 0 ] || { echo "FAILED unfold $CFG $OBS rc=$rc" >&2; exit 1; }

sig=$(grep 'SYSTDUMP] sigma_int' slurm-*.out 2>/dev/null | tail -1 | awk '{print $3}')
dv=$(grep 'SYSTDUMP] detVar_total' "$SLURM_SUBMIT_DIR/unf1p_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" 2>/dev/null | tail -1 | awk '{print $3}')

# Harvest the matrices while they are still this task's own.
OUTD="$RELEASE/cov/1p_${T}_${OBS}"
mkdir -p "$OUTD"
n=0
for f in unfold_output/XsecUnits_mat_table_cov_*.txt; do
  [ -f "$f" ] || continue
  cp -p "$f" "$OUTD/$(basename "$f" | sed 's/^XsecUnits_mat_table_//')"
  n=$((n+1))
done
[ -f unfold_output/mat_table_add_smear.txt ] && cp -p unfold_output/mat_table_add_smear.txt "$OUTD/add_smear.txt"
echo "  harvested $n covariance matrices -> $OUTD"

cd / && rm -rf "$WORK"
echo "OK $CFG $OBS in $(( (t1-t0)/60 ))min  detVar_total=${dv:-?}%"
