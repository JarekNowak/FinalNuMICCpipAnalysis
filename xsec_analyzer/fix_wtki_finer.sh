#!/usr/bin/env bash
# fix_wtki_finer.sh — CC1mu1pi1p (proton-tagged) reprocess at finer binning + 5 added
# kinematic variables, for FHC/RHC/comb. Now that the per-run POT bug is fixed, the
# W/TKI binning is restored from 3 to 6 equal-population bins (matching the inclusive
# granularity), and the 5 inclusive kinematic observables (pmu,ppi,costhmu,costhpi,thmupi)
# are added to the proton-tagged selection. univmake bakes in the binning, so the stale
# 3-bin W/TKI univmakes are deleted first. Resumable (skips univmakes already >30MB).
# NOTE: the per-config TAG is resolved in the PARENT and passed as an argument, because
# an associative-array lookup inside a backgrounded (&) function subshell misresolves.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
WTKI=(Wpipr Whad dpt dalphat dphit pn)         # rebinned 3->6: delete stale univmakes
KIN=(pmu ppi costhmu costhpi thmupi)           # newly added kinematic observables
ALL=("${WTKI[@]}" "${KIN[@]}")
rm -f $LOG/wf.status

# delete stale/wrong univmakes for the rebinned W/TKI so they rebuild at 6-bin
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}; for k in "${WTKI[@]}"; do rm -f $PROC/ccpi1p_${t}_${k}_univmake.root; done; done

univ_one(){  # cfg t k uout   (all resolved in parent; no assoc-array access here)
  local cfg=$1 t=$2 k=$3 uout=$4
  [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 30000000 ] && { echo "[skip] $t $k"; return 0; }
  echo "[$(date +%H:%M)] START univmake $t $k -> $(basename $uout)"
  FPM=configs/file_properties_numi_${cfg}_w.txt BIN_CONFIG=configs/ccpi1p_${k}_bin_config.txt \
    OUT="$uout" ./run_universe_maker.sh > $LOG/wf_${cfg}_${k}_univ.log 2>&1
  echo "[$(date +%H:%M)] END   univmake $t $k ($(($(stat -c%s $uout 2>/dev/null||echo 0)/1000000))MB)"
}
echo "==== W/TKI FINER-BINNING REPROCESS START $(date) NPAR=$NPAR ===="
for cfg in fhc5 rhcfull comb; do
  t=${TAG[$cfg]}
  echo "==== $cfg ($t) univmakes $(date) ===="
  for k in "${ALL[@]}"; do
    univ_one "$cfg" "$t" "$k" "$PROC/ccpi1p_${t}_${k}_univmake.root" &
    while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done
  done
  wait
  echo "==== $cfg ($t) unfolds $(date) ===="
  for k in "${ALL[@]}"; do
    uout=$PROC/ccpi1p_${t}_${k}_univmake.root
    [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $t $k univmake missing"; echo "DONE_${cfg}_${k}" >> $LOG/wf.status; continue; }
    nice bin/UnfolderNuMI configs/ccpi1p_xsec_config_numi_${k}_${cfg}.txt configs/ccpi1p_${k}_slice_config.txt \
      $PROC/xsec_ccpi1p_${t}_${k}.root > $LOG/wf_${cfg}_${k}_unfold.log 2>&1 \
      && echo "  OK $t $k  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/wf_${cfg}_${k}_unfold.log|tail -1|awk '{print $3}')" \
      || echo "  !! unfold $t $k FAILED"
    echo "DONE_${cfg}_${k}" >> $LOG/wf.status
  done
done
echo "ALL_DONE" >> $LOG/wf.status
echo "==== W/TKI FINER-BINNING REPROCESS DONE $(date) ===="
