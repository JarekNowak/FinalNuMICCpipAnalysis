#!/bin/bash
# SLURM array task: universe file (univmake only) of ONE observable/configuration of the binnings adopted
# under the 0.50 migration criterion (2026-09-26), from the staged configs in configs/adopt050/, into
# $PROC/rebuild_adopt050/ under the release file name. promote_050.sh installs configs and files together
# and runs the release chain (which unfolds with the generator predictions). Nothing live is touched here.
#SBATCH --job-name=adopt050
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; MANIFEST="${MANIFEST:?}"
RB="$PROC/rebuild_adopt050"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r PFX CFG OBS BIN <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${BIN:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; export UNIV_THREADS="${UNIV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB ); T="${TAG[$CFG]}"; mkdir -p "$RB"
if [ "$PFX" = ccpi1p ]; then FPM="configs/file_properties_numi_${CFG}_w.txt"; else FPM="configs/file_properties_numi_${CFG}.txt"; fi
B="configs/adopt050/$BIN"
for f in "$FPM" "$B"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
grep -q 'ext_run' "$FPM" || { echo "ERROR: $FPM is not the per-run EXT configuration" >&2; exit 1; }
echo "host=$(hostname) pfx=$PFX cfg=$CFG obs=$OBS bin=$B fpm=$FPM"
U="$RB/${PFX}_${T}_${OBS}_univmake.root"; rm -f "$U"
./bin/univmake "$FPM" "$B" "$U" "$FPM" > "$RB/univ_${PFX}_${T}_${OBS}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $PFX $CFG $OBS" >&2; exit 1; }
echo "OK $PFX $CFG $OBS $(stat -c %s "$U")"
