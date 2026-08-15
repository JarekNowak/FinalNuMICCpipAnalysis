#!/usr/bin/env bash
# fix_revert_reunfold.sh — re-unfold every RHC + COMB result with the REVERTED (per-file
# file_pot) POT scaling. The f60547d "summed run_type_mc_pot" fix was a misdiagnosis that
# under-scaled multi-file RHC/COMB runs; reverted in SystematicsCalculator.cxx. Univmakes
# need NO rebuild (raw histograms; scaling is read-time). XSEC_FORCE_REBUILD=1 busts each
# univmake's stale cached total_ subdir so it re-sums with the corrected scaling.
#   A: inclusive RHC + COMB   (5 obs each)
#   B: W/TKI RHC              (11 obs)   [W/TKI COMB handled by the running reprocess]
# FHC immune (single file/run) -> skipped.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; export XSEC_FORCE_REBUILD=1
PROC=/data/uboone/processed; LOG=../logs
declare -A TAG=( [rhcfull]=RHCFULL [comb]=COMB )
INCL=(pmu ppi costhmu costhpi thmupi)
WTKI=(Wpipr Whad dpt dalphat dphit pn pmu ppi costhmu costhpi thmupi)
rm -f $LOG/revert_reunfold.status

echo "==== REVERT RE-UNFOLD START $(date) ===="
echo "==== A: inclusive RHC + COMB ===="
for cfg in rhcfull comb; do t=${TAG[$cfg]}
  for o in "${INCL[@]}"; do
    nice -n 15 bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_${cfg}.txt \
      configs/ccpi_${o}_slice_config_opt.txt $PROC/xsec_${t}_${o}.root \
      > $LOG/revert_incl_${cfg}_${o}.log 2>&1 \
      && echo "  OK incl $t $o  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/revert_incl_${cfg}_${o}.log|tail -1|awk '{print $3}')" \
      || echo "  !! incl $t $o FAILED"
    echo "DONE_incl_${cfg}_${o}" >> $LOG/revert_reunfold.status
  done
done
echo "==== B: W/TKI RHC ===="
for k in "${WTKI[@]}"; do
  u=$PROC/ccpi1p_RHCFULL_${k}_univmake.root
  [ "$(stat -c%s "$u" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! W/TKI RHC $k univmake missing"; continue; }
  nice -n 15 bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_rhcfull.txt \
    configs/ccpi1p_${k}_slice_config.txt $PROC/xsec_ccpi1p_RHCFULL_${k}.root \
    > $LOG/revert_wtki_rhc_${k}.log 2>&1 \
    && echo "  OK wtki RHCFULL $k  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/revert_wtki_rhc_${k}.log|tail -1|awk '{print $3}')" \
    || echo "  !! wtki RHCFULL $k FAILED"
  echo "DONE_wtki_rhc_${k}" >> $LOG/revert_reunfold.status
done
echo "ALL_DONE" >> $LOG/revert_reunfold.status
echo "==== REVERT RE-UNFOLD DONE $(date) ===="
