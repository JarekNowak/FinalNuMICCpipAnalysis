#!/usr/bin/env bash
# run_comb_wtki.sh — finish the combined proton-tagged (W/TKI) measurement.
# Four univmakes (Wpipr, Whad, dpt, dalphat) survived the earlier stopped reprocess and are
# unfolded separately; this rebuilds the seven that were truncated or never built and
# unfolds them. Combined runs read 14 MC files, so each univmake is slow: 2-wide.
# Uses the convention-agnostic MC POT normalisation (sum of DISTINCT per-file summed_pot),
# which is what makes the proton-tagged RHC/combined results come out positive.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
TODO=(dphit pn pmu ppi costhmu costhpi thmupi)
rm -f $LOG/comb_wtki.status

univ_one(){
  local k=$1 out=$PROC/ccpi1p_COMB_${k}_univmake.root
  if [ "$(stat -c%s "$out" 2>/dev/null||echo 0)" -gt 30000000 ]; then echo "[skip] $k (exists)"; return 0; fi
  echo "[$(date +%H:%M)] START univmake COMB $k"
  FPM=configs/file_properties_numi_comb_w.txt BIN_CONFIG=configs/ccpi1p_${k}_bin_config.txt \
    OUT="$out" nice -n 10 ./run_universe_maker.sh > $LOG/comb_wtki_${k}_univ.log 2>&1
  echo "[$(date +%H:%M)] END   univmake COMB $k ($(($(stat -c%s "$out" 2>/dev/null||echo 0)/1000000))MB)"
}

echo "==== COMB W/TKI START $(date) NPAR=$NPAR ===="
for k in "${TODO[@]}"; do
  univ_one "$k" &
  while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done
done
wait
echo "==== COMB W/TKI unfolds $(date) ===="
for k in "${TODO[@]}"; do
  u=$PROC/ccpi1p_COMB_${k}_univmake.root
  [ "$(stat -c%s "$u" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $k univmake missing"; continue; }
  XSEC_FORCE_REBUILD=1 nice -n 10 bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_comb.txt \
    configs/ccpi1p_${k}_slice_config.txt $PROC/xsec_ccpi1p_COMB_${k}.root \
    > $LOG/comb_wtki_${k}_unfold.log 2>&1 \
    && echo "  OK COMB $k sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/comb_wtki_${k}_unfold.log|tail -1|awk '{print $3}')" \
    || echo "  !! unfold $k FAILED"
  echo "DONE_${k}" >> $LOG/comb_wtki.status
done
echo "ALL_DONE" >> $LOG/comb_wtki.status
echo "==== COMB W/TKI DONE $(date) ===="
