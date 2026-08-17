#!/usr/bin/env bash
# run_ppi2bin_rest.sh — RHC and combined p_pi at the adopted two-bin binning
# ([0.175,0.205] + everything above). FHC is already done. Generator overlays use
# the corrected (untruncated) FTE predictions.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs
declare -A TAG=( [rhcfull]=RHCFULL [comb]=COMB )
rm -f $LOG/ppi2bin_rest.status
one(){
  local cfg=$1 t=${TAG[$1]}
  local u=$PROC/ccpi_${t}_ppi2bin_univmake.root
  if [ "$(stat -c%s "$u" 2>/dev/null||echo 0)" -lt 30000000 ]; then
    echo "[$(date +%H:%M)] START univmake $t"
    nice -n 10 bin/univmake configs/file_properties_numi_${cfg}.txt configs/ccpi_ppi_bin_config_2bin.txt \
      "$u" > $LOG/ppi2bin_${cfg}_univ.log 2>&1
    echo "[$(date +%H:%M)] END   univmake $t ($(( $(stat -c%s "$u" 2>/dev/null||echo 0)/1048576 ))MB)"
  else echo "[skip] univmake $t exists"; fi
  nice -n 10 bin/UnfolderNuMI configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt \
    configs/ccpi_ppi_slice_config_2bin.txt $PROC/xsec_${t}_ppi2bin.root \
    > $LOG/ppi2bin_${cfg}_unfold.log 2>&1 \
    && echo "  OK $t sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/ppi2bin_${cfg}_unfold.log|tail -1|awk '{print $3}')" \
    || echo "  !! unfold $t FAILED"
  echo "DONE_${cfg}" >> $LOG/ppi2bin_rest.status
}
echo "==== ppi 2-bin RHC+COMB START $(date) ===="
one rhcfull & one comb &
wait
echo "ALL_DONE" >> $LOG/ppi2bin_rest.status
echo "==== ppi 2-bin RHC+COMB DONE $(date) ===="
