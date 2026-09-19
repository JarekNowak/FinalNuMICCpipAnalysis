#!/usr/bin/env bash
# ens_statcov.sh CFG OBS [NMAX] -- re-run UnfolderNuMI on each ensemble member's universe file with a
# KEPT work directory, so the per-bin statistical covariance tables (mat_table_cov_DataStats.txt etc.)
# are available for the statistics-only pull test; the SLURM ensemble job deletes its work directory.
# Members already done are skipped, so the script can be re-run while the ensemble is still filling.
#   CFG = fhc5 | rhcfull | comb ; OBS = pmu | costhmu | ppi2bin | total ; NMAX default 100.
# (2026-09-19: generalised from the FHC-only form, which read xsec_fhc_OBS_tT.txt and wrote statcov_OBS_tT.)
set -u
CFG=${1:?cfg}; OBS=${2:?obs}; NMAX=${3:-100}
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; PROC=/data/uboone/processed
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
case "$OBS" in ppi2bin) SL=configs/ccpi_ppi_slice_config_2bin.txt;; total) SL=configs/ccpi_total_slice_config_opt.txt;; *) SL=configs/ccpi_${OBS}_slice_config_opt.txt;; esac
ok=0; skip=0; fail=0
for T in $(seq 1 "$NMAX"); do
  XC="$PROC/ens/xsec_${CFG}_${OBS}_t${T}.txt"; [ -f "$XC" ] || continue
  [ -s "$PROC/ens/closure_hists_xsec_${CFG}_${OBS}_t${T}.root" ] || continue   # member not finished
  W="$PROC/ens/statcov_${CFG}_${OBS}_t${T}"
  [ -f "$W/unfold_output/mat_table_cov_DataStats.txt" ] && { skip=$((skip+1)); continue; }
  rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SL" "$W/xsec_${CFG}_${OBS}_t${T}.root" ) > "$W/run.log" 2>&1
  if [ -f "$W/unfold_output/mat_table_cov_DataStats.txt" ]; then ok=$((ok+1)); else echo "FAIL $T"; fail=$((fail+1)); fi
done
echo "STATCOV_DONE $CFG $OBS: $ok new, $skip already done, $fail failed  $(date +%H:%M)"
