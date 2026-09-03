#!/bin/bash
# Local (sequential) counterpart of submit_unfold_1p.sh: re-unfold the 33 proton-tagged
# extractions and harvest covariance + A_C for the data release. Same per-extraction
# body as slurm_unfold_1p.sbatch; used after a code change that leaves the universes
# valid but refreshes the closure sidecars (e.g. the fake-data-truth weighting fix).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RELEASE="$(cd "$REPO/../report" && pwd)/data_release"
PROC=/data/uboone/processed; SCRATCH=$PROC/unfold_scratch
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
ok=0; fail=0
while read -r CFG OBS; do
  T=${TAG[$CFG]}
  case "$OBS" in *2bin) SLICE="configs/ccpi1p_${OBS%2bin}_slice_config_2bin.txt" ;; *) SLICE="configs/ccpi1p_${OBS}_slice_config.txt" ;; esac
  XSEC="configs/ccpi1p_xsec_config_numi_${OBS}_${CFG}.txt"
  W="$SCRATCH/1p_${T}_${OBS}.$$"; rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XSEC" "$SLICE" "$PROC/xsec_ccpi1p_${T}_${OBS}.root" ) > "$W/run.log" 2>&1
  rc=$?; n=$(ls "$W"/unfold_output/mat_table_cov_*.txt 2>/dev/null | wc -l)
  if [ $rc -ne 0 ] || [ "$n" -eq 0 ]; then echo "  FAILED $CFG $OBS (rc=$rc, $n)"; fail=$((fail+1)); rm -rf "$W"; continue; fi
  OUTD="$RELEASE/cov/1p_${T}_${OBS}"; mkdir -p "$OUTD"
  for f in "$W"/unfold_output/mat_table_cov_*.txt; do cp -p "$f" "$OUTD/$(basename "$f" | sed 's/^mat_table_//')"; done
  for extra in add_smear unfolding err_prop; do [ -f "$W/unfold_output/mat_table_${extra}.txt" ] && cp -p "$W/unfold_output/mat_table_${extra}.txt" "$OUTD/${extra}.txt"; done
  [ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$OUTD/unfolded_signal.txt"
  echo "  OK $CFG $OBS ($n cov matrices)"; ok=$((ok+1)); rm -rf "$W"
done < "$REPO/slurm/unfold_1p_manifest.list"
echo "  ---- $ok ok, $fail failed ----"
