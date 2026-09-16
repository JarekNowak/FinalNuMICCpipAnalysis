#!/bin/bash
# SLURM array task: univmake + UnfolderNuMI for ONE (xsec config, bin config, slice config) triple
# from a manifest ("XSEC_CONFIG BIN_CONFIG SLICE_CONFIG", paths relative to the repo). The universe
# file and the file-properties list are taken FROM the xsec config (UnivFile / FPFile), so the
# manifest cannot point a build at the wrong place. Unfolding output goes next to the universe file
# as xsec_<univ basename minus _univmake>.root, with its log alongside.
#SBATCH --job-name=extract
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r XC BIN SLICE <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${SLICE:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; export UNIV_THREADS=16
for f in "$XC" "$BIN" "$SLICE"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
U=$(awk '$1=="UnivFile"{print $2}' "$XC"); FPM=$(awk '$1=="FPFile"{print $2}' "$XC")
[ -n "$U" ] && [ -f "$FPM" ] || { echo "ERROR: $XC lacks UnivFile/FPFile" >&2; exit 1; }
grep -q 'ext_run' "$FPM" || { echo "ERROR: $FPM is not the per-run EXT configuration" >&2; exit 1; }
D=$(dirname "$U"); B=$(basename "$U" _univmake.root); mkdir -p "$D"
echo "host=$(hostname) xc=$XC bin=$BIN slice=$SLICE univ=$U"
rm -f "$U"
./bin/univmake "$FPM" "$BIN" "$U" "$FPM" > "$D/univ_${B}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $XC" >&2; exit 1; }
W="$D/work_${B}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$REPO/$XC" "$REPO/$SLICE" "$D/xsec_${B}.root" ) > "$D/unfold_${B}.log" 2>&1
rc=$?; mkdir -p "$D/plots_${B}"; cp "$W"/unfold_output/*.pdf "$D/plots_${B}/" 2>/dev/null; rm -rf "$W"
[ $rc -eq 0 ] || { echo "FAILED unfold $XC rc=$rc" >&2; exit 1; }
echo "OK $XC sigma_int=$(grep 'SYSTDUMP] sigma_int' "$D/unfold_${B}.log" | tail -1 | awk '{print $3}')"
