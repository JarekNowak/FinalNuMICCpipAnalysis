#!/usr/bin/env bash
# w_batch_fhc.sh — CC1mu1pi1p (W/TKI) FHC univmake + unfold batch, 2-wide. Reads the
# w/-reprocessed FHC files via file_properties_numi_fhc5_w.txt and the per-observable
# ccpi1p bin/slice/xsec configs. detVar-free systcalc -> flux/xsec/reint systematics.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
obs=(Wpipr Whad dpt dalphat dphit pn)

univ_one() {  # key
  local k=$1 uout=$PROC/ccpi1p_FHC5_${k}_univmake.root
  local sz=$(stat -c%s "$uout" 2>/dev/null || echo 0)
  if [ "$sz" -gt 30000000 ]; then echo "[skip] $k (done $((sz/1000000))MB)"; return 0; fi
  echo "[$(date +%H:%M)] START univmake $k"
  FPM=configs/file_properties_numi_fhc5_w.txt BIN_CONFIG=configs/ccpi1p_${k}_bin_config.txt \
    OUT=$uout ./run_universe_maker.sh > $LOG/w_fhc_${k}_univ.log 2>&1
  echo "[$(date +%H:%M)] END   univmake $k ($(($(stat -c%s $uout 2>/dev/null||echo 0)/1000000))MB)"
}
echo "==== W/TKI FHC BATCH START $(date) NPAR=$NPAR ===="
for k in "${obs[@]}"; do
  univ_one "$k" &
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done
done
wait
for k in "${obs[@]}"; do
  uout=$PROC/ccpi1p_FHC5_${k}_univmake.root
  [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $k univmake missing"; continue; }
  bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_fhc5.txt configs/ccpi1p_${k}_slice_config.txt \
    $PROC/xsec_ccpi1p_FHC5_${k}.root > $LOG/w_fhc_${k}_unfold.log 2>&1 \
    && echo "  OK $k  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/w_fhc_${k}_unfold.log | tail -1 | awk '{print $3}')" \
    || echo "  !! unfold $k FAILED"
done
echo "==== W/TKI FHC BATCH DONE $(date) ===="
