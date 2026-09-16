#!/bin/bash
# SLURM array task: full extraction (univmake + UnfolderNuMI) of ONE candidate binning of the
# 0.50 migration-diagonal study (2026-09-16) into a SEPARATE tree, $PROC/rebuild_c50. Configs are
# written by scripts/make_bin_configs.py into configs/c50/. Nothing in the live tree is touched.
#SBATCH --job-name=c50
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; MANIFEST="${MANIFEST:?}"
RB="$PROC/rebuild_c50"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r PFX CFG OBS <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${OBS:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; export UNIV_THREADS=16
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL ); T="${TAG[$CFG]}"; mkdir -p "$RB"
if [ "$PFX" = ccpi1p ]; then FPM="configs/file_properties_numi_${CFG}_w.txt"; else FPM="configs/file_properties_numi_${CFG}.txt"; fi
BIN="configs/c50/${PFX}_${OBS}_bin_config_c50.txt"; SLICE="configs/c50/${PFX}_${OBS}_slice_config_c50.txt"
XC="configs/c50/${PFX}_xsec_config_numi_${OBS}_c50_${CFG}.txt"
for f in "$FPM" "$BIN" "$SLICE" "$XC"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
grep -q 'ext_run' "$FPM" || { echo "ERROR: $FPM is not the per-run EXT configuration" >&2; exit 1; }
grep -q "^UnivFile $RB/" "$XC" || { echo "ERROR: $XC does not point into $RB" >&2; exit 1; }
echo "host=$(hostname) pfx=$PFX cfg=$CFG obs=$OBS fpm=$FPM"
U="$RB/${PFX}_${T}_${OBS}_univmake.root"; rm -f "$U"
./bin/univmake "$FPM" "$BIN" "$U" "$FPM" > "$RB/univ_${PFX}_${T}_${OBS}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $PFX $CFG $OBS" >&2; exit 1; }
W="$RB/work_${PFX}_${T}_${OBS}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$REPO/$XC" "$REPO/$SLICE" "$RB/xsec_${PFX}_${T}_${OBS}.root" ) > "$RB/unfold_${PFX}_${T}_${OBS}.log" 2>&1
rc=$?; rm -rf "$W"
[ $rc -eq 0 ] || { echo "FAILED unfold $PFX $CFG $OBS rc=$rc" >&2; exit 1; }
sig=$(grep 'SYSTDUMP] sigma_int' "$RB/unfold_${PFX}_${T}_${OBS}.log" | tail -1 | awk '{print $3}')
echo "OK $PFX $CFG $OBS sigma_int=$sig"
