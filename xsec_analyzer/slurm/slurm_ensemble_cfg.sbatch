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
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; OBS="${OBS:-pmu}"; CFG="${CFG:-fhc5}"
export UNIV_THREADS="${UNIV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
case "$CFG" in fhc5) MODES="fhc"; RUNS_fhc="1 2 4 5";; rhcfull) MODES="rhc"; RUNS_rhc="1 2 3 4";; comb) MODES="fhc rhc"; RUNS_fhc="1 2 4 5"; RUNS_rhc="1 2 3 4";; esac
T=$(( SLURM_ARRAY_TASK_ID + 1 ))

cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"

echo "host=$(hostname) throw=$T obs=$OBS"

# --- 1. throw -----------------------------------------------------------------------
for m in $MODES; do
  root -l -b -q "macros/throw_ensemble_${m}.C($T)" 2>&1 | grep POTSCALE
  eval runs=\$RUNS_$m
  for r in $runs; do [ -s "$PROC/ens/fakedata_${m}_run${r}_t${T}.root" ] || { echo "FAILED throw $T $m run $r" >&2; exit 1; }; done
done

# --- 2. a file_properties pointing at this member's data ----------------------------
FPM="$PROC/ens/fp_${CFG}_t${T}.txt"
sed -E "s#(xsec-ana-)fakedata_(fhc|rhc)_run([0-9])\.root#ens/fakedata_\2_run\3_t${T}.root#" \
  configs/file_properties_numi_${CFG}.txt > "$FPM"
n=$(grep -c "_t${T}\.root" "$FPM"); want=4; [ "$CFG" = comb ] && want=8
[ "$n" -eq $want ] || { echo "FAILED fp rewrite for throw $T (matched $n of $want)" >&2; exit 1; }

case "$OBS" in
  ppi2bin) BIN="configs/ccpi_ppi_bin_config_2bin.txt"; SLICE="configs/ccpi_ppi_slice_config_2bin.txt";;
  *)       BIN="configs/ccpi_${OBS}_bin_config_opt.txt"; SLICE="configs/ccpi_${OBS}_slice_config_opt.txt";;
esac
[ -f "$BIN" ] && [ -f "$SLICE" ] || { echo "FAILED: no bin/slice config for $OBS" >&2; exit 1; }
# --- 3. universes from THIS throw ---------------------------------------------------
U="$PROC/ens/univ_${CFG}_${OBS}_t${T}.root"
rm -f "$U"
./bin/univmake "$FPM" "$BIN" "$U" "$FPM" > "$PROC/ens/univ_${OBS}_t${T}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake throw $T" >&2; exit 1; }

# --- 4. unfold ----------------------------------------------------------------------
XC="$PROC/ens/xsec_${CFG}_${OBS}_t${T}.txt"
sed -E "s#^UnivFile .*#UnivFile $U#; s#^FPFile .*#FPFile $FPM#" \
  "configs/ccpi_xsec_config_numi_${OBS}_${CFG}.txt" > "$XC"
W="$PROC/ens/work_${OBS}_t${T}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SLICE" \
    "$PROC/ens/xsec_${CFG}_${OBS}_t${T}.root" ) > "$PROC/ens/unfold_${CFG}_${OBS}_t${T}.log" 2>&1
rc=$?
SIDE="$PROC/closure_hists_xsec_${OBS}_t${T}.root"  # sidecar name is set by the framework
[ $rc -eq 0 ] || { echo "FAILED unfold throw $T rc=$rc" >&2; rm -rf "$W"; exit 1; }
rm -rf "$W"
echo "OK throw $T sigma_int=$(grep 'SYSTDUMP] sigma_int' "$PROC/ens/unfold_${CFG}_${OBS}_t${T}.log" | tail -1 | awk '{print $3}')"
