#!/bin/bash
# One ensemble member: throw the fake data, rebuild the universes from it, unfold, and
# record the unfolded spectrum with its uncertainty and the A_C-smeared truth.
#
# Each member is a genuinely independent Poisson realisation of the data, put through the
# whole chain. That is the point: the closure test in the note uses ONE throw, so it can
# say the extraction reproduces its input but not whether the quoted uncertainty is the
# right size. Only an ensemble gives pull widths and interval coverage.
#SBATCH --job-name=ens
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; OBS="${OBS:-pmu}"
T=$(( SLURM_ARRAY_TASK_ID + 1 ))

cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"

echo "host=$(hostname) throw=$T obs=$OBS"

# --- 1. throw -----------------------------------------------------------------------
root -l -b -q "macros/throw_ensemble_fhc.C($T)" 2>&1 | grep POTSCALE
for r in 1 2 4 5; do
  [ -s "$PROC/ens/fakedata_fhc_run${r}_t${T}.root" ] || { echo "FAILED throw $T run $r" >&2; exit 1; }
done

# --- 2. a file_properties pointing at this member's data ----------------------------
FPM="$PROC/ens/fp_fhc_t${T}.txt"
sed -E "s#(xsec-ana-)fakedata_fhc_run([0-9])\.root#ens/fakedata_fhc_run\2_t${T}.root#" \
  configs/file_properties_numi_fhc5.txt > "$FPM"
n=$(grep -c "_t${T}\.root" "$FPM")
[ "$n" -eq 4 ] || { echo "FAILED fp rewrite for throw $T (matched $n of 4)" >&2; exit 1; }

# --- 3. universes from THIS throw ---------------------------------------------------
U="$PROC/ens/univ_fhc_${OBS}_t${T}.root"
rm -f "$U"
./bin/univmake "$FPM" "configs/ccpi_${OBS}_bin_config_opt.txt" "$U" "$FPM" > "$PROC/ens/univ_t${T}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake throw $T" >&2; exit 1; }

# --- 4. unfold ----------------------------------------------------------------------
XC="$PROC/ens/xsec_fhc_${OBS}_t${T}.txt"
sed -E "s#^UnivFile .*#UnivFile $U#; s#^FPFile .*#FPFile $FPM#" \
  "configs/ccpi_xsec_config_numi_${OBS}_fhc5.txt" > "$XC"
W="$PROC/ens/work_t${T}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/configs/ccpi_${OBS}_slice_config_opt.txt" \
    "$PROC/ens/xsec_t${T}.root" ) > "$PROC/ens/unfold_t${T}.log" 2>&1
rc=$?
SIDE="$PROC/closure_hists_xsec_t${T}.root"
[ $rc -eq 0 ] || { echo "FAILED unfold throw $T rc=$rc" >&2; rm -rf "$W"; exit 1; }
rm -rf "$W"
echo "OK throw $T sigma_int=$(grep 'SYSTDUMP] sigma_int' "$PROC/ens/unfold_t${T}.log" | tail -1 | awk '{print $3}')"
