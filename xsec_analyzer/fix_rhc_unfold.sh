#!/usr/bin/env bash
# fix_rhc_unfold.sh — re-unfold RHC W/TKI observables after the per-run POT
# normalization fix. Uses XSEC_FORCE_REBUILD so build_universes() recomputes the
# POT-summed universes with the corrected multi-file-per-run scaling instead of
# loading the stale cache. Observables passed as args (default: all 6). Serial to
# limit disk I/O contention with the running COMB pipeline.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
export XSEC_FORCE_REBUILD=1
PROC=/data/uboone/processed
OBS=("$@"); [ ${#OBS[@]} -eq 0 ] && OBS=(Wpipr Whad dpt dalphat dphit pn)
echo "==== FIX RHC UNFOLD (force-rebuild, per-run POT) START $(date) ===="
for k in "${OBS[@]}"; do
  uni=$PROC/ccpi1p_RHCFULL_${k}_univmake.root
  [ -s "$uni" ] || { echo "  !! $k univmake missing"; continue; }
  nice bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_rhcfull.txt configs/ccpi1p_${k}_slice_config.txt \
    $PROC/xsec_ccpi1p_RHCFULL_${k}.root > logs/w_rhcfix_${k}_unfold.log 2>&1 \
    && echo "  OK RHCFULL $k  sigma_int=$(grep 'SYSTDUMP] sigma_int' logs/w_rhcfix_${k}_unfold.log | tail -1 | awk '{print $3}')  rebuilt=$(grep -c 'PROCESSING universes' logs/w_rhcfix_${k}_unfold.log)" \
    || echo "  !! RHCFULL $k unfold FAILED (see logs/w_rhcfix_${k}_unfold.log)"
  echo "DONE_$k" >> logs/fix_rhc.status
done
echo "==== FIX RHC UNFOLD DONE $(date) ===="
echo "ALL_DONE" >> logs/fix_rhc.status
