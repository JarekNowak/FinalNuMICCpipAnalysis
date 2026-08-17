#!/usr/bin/env bash
# w_batch_rhc_comb.sh — CC1mu1pi1p (W/TKI) univmake + unfold for RHC then combined,
# 2-wide within each config. Reads the w/-reprocessed files; detVar-free systcalc.
# comb reuses the FHC+RHC w/ files (14 numuMC) so it is the slow config.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
obs=(Wpipr Whad dpt dalphat dphit pn)

run_cfg() {  # cfg TAG
  local cfg=$1 TAG=$2
  univ_one() { local k=$1 uout=$PROC/ccpi1p_${TAG}_${k}_univmake.root
    local sz=$(stat -c%s "$uout" 2>/dev/null || echo 0)
    if [ "$sz" -gt 30000000 ]; then echo "[skip] $TAG $k (done $((sz/1000000))MB)"; return 0; fi
    echo "[$(date +%H:%M)] START univmake $TAG $k"
    FPM=configs/file_properties_numi_${cfg}_w.txt BIN_CONFIG=configs/ccpi1p_${k}_bin_config.txt \
      OUT=$uout ./run_universe_maker.sh > $LOG/w_${cfg}_${k}_univ.log 2>&1
    echo "[$(date +%H:%M)] END   univmake $TAG $k ($(($(stat -c%s $uout 2>/dev/null||echo 0)/1000000))MB)"; }
  echo "==== $TAG univmakes $(date) ===="
  for k in "${obs[@]}"; do univ_one "$k" & while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done; done
  wait
  for k in "${obs[@]}"; do
    uout=$PROC/ccpi1p_${TAG}_${k}_univmake.root
    [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $TAG $k univmake missing"; continue; }
    bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_${cfg}.txt configs/ccpi1p_${k}_slice_config.txt \
      $PROC/xsec_ccpi1p_${TAG}_${k}.root > $LOG/w_${cfg}_${k}_unfold.log 2>&1 \
      && echo "  OK $TAG $k  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/w_${cfg}_${k}_unfold.log | tail -1 | awk '{print $3}')" \
      || echo "  !! unfold $TAG $k FAILED"
  done
}
echo "==== W/TKI RHC+COMB BATCH START $(date) NPAR=$NPAR ===="
run_cfg rhcfull RHCFULL
run_cfg comb    COMB
echo "==== W/TKI RHC+COMB BATCH DONE $(date) ===="
