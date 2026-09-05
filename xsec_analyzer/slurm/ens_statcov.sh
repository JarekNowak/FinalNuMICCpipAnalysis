#!/usr/bin/env bash
# ens_statcov.sh OBS -- re-run UnfolderNuMI on each ensemble member's universe file with a KEPT work
# directory, so the per-bin statistical covariance tables (mat_table_cov_DataStats.txt etc.) are
# available for the statistics-only pull test; the SLURM ensemble job deletes its work directory.
set -u
OBS=${1:?obs}; REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; PROC=/data/uboone/processed
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
SL=configs/ccpi_${OBS}_slice_config_opt.txt
for T in $(seq 1 32); do
  XC="$PROC/ens/xsec_fhc_${OBS}_t${T}.txt"; [ -f "$XC" ] || { echo "no config for throw $T"; continue; }
  W="$PROC/ens/statcov_${OBS}_t${T}"; rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SL" "$W/xsec_${OBS}_t${T}.root" ) > "$W/run.log" 2>&1
  [ -f "$W/unfold_output/mat_table_cov_DataStats.txt" ] && echo "OK $T" || echo "FAIL $T"
done
echo "STATCOV_DONE $OBS"
