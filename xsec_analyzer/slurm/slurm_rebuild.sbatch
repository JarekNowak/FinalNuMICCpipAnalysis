#!/bin/bash
# Clean-rebuild one extraction from the processed ntuples, into a SEPARATE tree.
#
# The point is to demonstrate that the release can be regenerated from raw inputs. It
# therefore must not touch the current univmakes, sidecars or release -- if it overwrote
# them, a failure halfway would destroy the thing it is trying to verify. Everything is
# built under $PROC/rebuild and compared afterwards.
#SBATCH --job-name=rebuild
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; MANIFEST="${MANIFEST:?}"
RB="$PROC/rebuild"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r FAM CFG OBS <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${OBS:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }

cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
T="${TAG[$CFG]}"
mkdir -p "$RB"

if [ "$FAM" = "1p" ]; then
  PFX=ccpi1p; FPM="configs/file_properties_numi_${CFG}_w.txt"
  case "$OBS" in
    *2bin) BIN="configs/${PFX}_${OBS%2bin}_bin_config_2bin.txt"
           SLICE="configs/${PFX}_${OBS%2bin}_slice_config_2bin.txt" ;;
    *)     BIN="configs/${PFX}_${OBS}_bin_config.txt"
           SLICE="configs/${PFX}_${OBS}_slice_config.txt" ;;
  esac
else
  PFX=ccpi; FPM="configs/file_properties_numi_${CFG}.txt"
  if [ "$OBS" = ppi2bin ]; then BIN="configs/ccpi_ppi_bin_config_2bin.txt"
                                SLICE="configs/ccpi_ppi_slice_config_2bin.txt"
  else BIN="configs/ccpi_${OBS}_bin_config_opt.txt"
       SLICE="configs/ccpi_${OBS}_slice_config_opt.txt"; fi
fi
XC="configs/${PFX}_xsec_config_numi_${OBS}_${CFG}.txt"
for f in "$FPM" "$BIN" "$SLICE" "$XC"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

echo "host=$(hostname) fam=$FAM cfg=$CFG obs=$OBS"
U="$RB/${PFX}_${T}_${OBS}_univmake.root"
rm -f "$U"                       # empty cache: no universes to reuse
./bin/univmake "$FPM" "$BIN" "$U" "$FPM" > "$RB/univ_${PFX}_${T}_${OBS}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $FAM $CFG $OBS" >&2; exit 1; }

XC2="$RB/xsec_${PFX}_${T}_${OBS}.txt"
sed -E "s#^UnivFile .*#UnivFile $U#" "$XC" > "$XC2"
W="$RB/work_${PFX}_${T}_${OBS}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC2" "$REPO/$SLICE" \
    "$RB/xsec_${PFX}_${T}_${OBS}.root" ) > "$RB/unfold_${PFX}_${T}_${OBS}.log" 2>&1
rc=$?
rm -rf "$W"
[ $rc -eq 0 ] || { echo "FAILED unfold $FAM $CFG $OBS rc=$rc" >&2; exit 1; }
sig=$(grep 'SYSTDUMP] sigma_int' "$RB/unfold_${PFX}_${T}_${OBS}.log" | tail -1 | awk '{print $3}')
echo "OK $FAM $CFG $OBS sigma_int=$sig"
