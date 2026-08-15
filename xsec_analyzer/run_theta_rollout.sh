#!/usr/bin/env bash
# run_theta_rollout.sh — full theta (theta_mu, theta_pi) rollout replacing cos(theta).
# Motivation: the Wiener-SVD additional-smearing matrix A_C suppresses the sharply
# forward-peaked cos(theta) distributions ~x0.33, crushing the (peaked) external
# generator predictions relative to the (smoother) CV/tune and making the angular
# panels look mis-normalised. In theta the sin(theta) Jacobian flattens the peak and
# the suppression drops to ~x0.76 -- comparable to the momentum observables.
# Bin edges are exactly acos() of the existing cos(theta) analysis edges (same phase
# space); generator FTE predictions were produced by exact bin reversal (sigma conserved).
# RHC runs FIRST (user request), then FHC theta_pi, then COMB.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
rm -f $LOG/theta_rollout.status

univ_one(){ # cfg TAG obs
  local cfg=$1 t=$2 o=$3
  local out=$PROC/ccpi_${t}_${o}_univmake.root
  if [ "$(stat -c%s "$out" 2>/dev/null||echo 0)" -gt 30000000 ]; then echo "[skip] $t $o (exists)"; return 0; fi
  echo "[$(date +%H:%M)] START univmake $t $o"
  nice -n 10 bin/univmake configs/file_properties_numi_${cfg}.txt configs/ccpi_${o}_bin_config_opt.txt \
    "$out" > $LOG/theta_${cfg}_${o}_univ.log 2>&1
  echo "[$(date +%H:%M)] END   univmake $t $o ($(($(stat -c%s "$out" 2>/dev/null||echo 0)/1000000))MB)"
}
unfold_one(){ # cfg TAG obs
  local cfg=$1 t=$2 o=$3
  local u=$PROC/ccpi_${t}_${o}_univmake.root
  [ "$(stat -c%s "$u" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $t $o univmake missing"; return 1; }
  XSEC_FORCE_REBUILD=1 nice -n 10 bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_${cfg}.txt \
    configs/ccpi_${o}_slice_config_opt.txt $PROC/xsec_${t}_${o}.root \
    > $LOG/theta_${cfg}_${o}_unfold.log 2>&1 \
    && echo "  OK $t $o  sigma_int=$(grep 'SYSTDUMP] sigma_int' $LOG/theta_${cfg}_${o}_unfold.log|tail -1|awk '{print $3}')" \
    || echo "  !! unfold $t $o FAILED"
  echo "DONE_${cfg}_${o}" >> $LOG/theta_rollout.status
}

echo "==== THETA ROLLOUT START $(date) NPAR=$NPAR ===="
for cfg in rhcfull fhc5 comb; do
  t=${TAG[$cfg]}
  echo "==== $cfg ($t) univmakes $(date) ===="
  for o in thetamu thetapi; do
    univ_one "$cfg" "$t" "$o" &
    while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done
  done
  wait
  echo "==== $cfg ($t) unfolds $(date) ===="
  for o in thetamu thetapi; do unfold_one "$cfg" "$t" "$o"; done
done
echo "ALL_DONE" >> $LOG/theta_rollout.status
echo "==== THETA ROLLOUT DONE $(date) ===="
