#!/bin/bash
# One PROTON-TAGGED ensemble member: throw the w/-tree fake data, rebuild the CC1mu1pi1p universes
# from it, unfold, and record the unfolded spectrum with its uncertainty and the A_C-smeared truth.
# Counterpart of slurm_ensemble_cfg.sbatch for the ccpi1p family (configs/file_properties_numi_<cfg>_w.txt,
# configs/ccpi1p_*): outputs carry the tag 1p_<cfg> so ensemble_stat_pulls.evaluate('1p_<cfg>_<obs>')
# and ens_statcov.sh 1p_<cfg> <obs> work unchanged. Throw files live under w/ens/throws/ (created once, never replaced: macros/throw_guard.h).
#SBATCH --job-name=ens1p
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out

set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; OBS="${OBS:-Whad}"; CFG="${CFG:-fhc5}"
export UNIV_THREADS="${UNIV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
case "$CFG" in fhc5) MODES="fhc"; RUNS_fhc="1 2 4 5";; rhcfull) MODES="rhc"; RUNS_rhc="1 2 3 4";; comb) MODES="fhc rhc"; RUNS_fhc="1 2 4 5"; RUNS_rhc="1 2 3 4";; esac
T=$(( SLURM_ARRAY_TASK_ID + 1 )); TAG="1p_${CFG}"
# CONF: directory of the bin, slice and xsec configs (configs/adopt050 runs the staged binnings before promotion)
CONF="${CONF:-configs}"

cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
mkdir -p "$PROC/w/ens/throws" "$PROC/ens/throws"
echo "host=$(hostname) throw=$T fam=1p cfg=$CFG obs=$OBS"

# --- 1. throw (w/ tree) --------------------------------------------------------------
for m in $MODES; do
  root -l -b -q "macros/throw_ensemble_w.C($T,\"$m\")" 2>&1 | grep -E 'POTSCALE|FAILED'
  eval runs=\$RUNS_$m
  for r in $runs; do [ -s "$PROC/w/ens/throws/fakedata_${m}_run${r}_t${T}.root" ] || { echo "FAILED throw $T $m run $r" >&2; exit 1; }; done
done

# --- 2. a file_properties pointing at this member's data ----------------------------
FPM="$PROC/ens/fp_${TAG}_${OBS}_t${T}.txt"   # per observable, written atomically (see slurm_ensemble_cfg.sbatch)
sed -E "s#/w/xsec-ana-fakedata_(fhc|rhc)_run([0-9])\.root#/w/ens/throws/fakedata_\1_run\2_t${T}.root#" \
  configs/file_properties_numi_${CFG}_w.txt > "$FPM.tmp.$$" && mv -f "$FPM.tmp.$$" "$FPM"
n=$(grep -c "_t${T}\.root" "$FPM"); want=4; [ "$CFG" = comb ] && want=8
[ "$n" -eq $want ] || { echo "FAILED fp rewrite for throw $T (matched $n of $want)" >&2; exit 1; }

case "$OBS" in
  *2bin) BIN="$CONF/ccpi1p_${OBS%2bin}_bin_config_2bin.txt"; SLICE="$CONF/ccpi1p_${OBS%2bin}_slice_config_2bin.txt";;
  *3bin) BIN="$CONF/ccpi1p_${OBS%3bin}_bin_config_3bin.txt"; SLICE="$CONF/ccpi1p_${OBS%3bin}_slice_config_3bin.txt";;   # 2026-09-26
  *)     BIN="$CONF/ccpi1p_${OBS}_bin_config.txt"; SLICE="$CONF/ccpi1p_${OBS}_slice_config.txt";;
esac
XC0="$CONF/ccpi1p_xsec_config_numi_${OBS}_${CFG}.txt"
[ -f "$BIN" ] && [ -f "$SLICE" ] && [ -f "$XC0" ] || { echo "FAILED: missing bin/slice/xsec config for 1p $CFG $OBS in $CONF" >&2; exit 1; }
echo "configs: $BIN $SLICE $XC0"
# --- 3. universes from THIS throw ---------------------------------------------------
U="$PROC/ens/univ_${TAG}_${OBS}_t${T}.root"
rm -f "$U"
./bin/univmake "$FPM" "$BIN" "$U" "$FPM" > "$PROC/ens/univ_${TAG}_${OBS}_t${T}.log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake throw $T" >&2; exit 1; }

# --- 4. unfold ----------------------------------------------------------------------
XC="$PROC/ens/xsec_${TAG}_${OBS}_t${T}.txt"
sed -E "s#^UnivFile .*#UnivFile $U#; s#^FPFile .*#FPFile $FPM#" "$XC0" > "$XC"
W="$PROC/ens/work_${TAG}_${OBS}_t${T}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SLICE" \
    "$PROC/ens/xsec_${TAG}_${OBS}_t${T}.root" ) > "$PROC/ens/unfold_${TAG}_${OBS}_t${T}.log" 2>&1
rc=$?
[ $rc -eq 0 ] || { echo "FAILED unfold throw $T rc=$rc" >&2; rm -rf "$W"; exit 1; }
rm -rf "$W"
echo "OK throw $T sigma_int=$(grep 'SYSTDUMP] sigma_int' "$PROC/ens/unfold_${TAG}_${OBS}_t${T}.log" | tail -1 | awk '{print $3}')"
