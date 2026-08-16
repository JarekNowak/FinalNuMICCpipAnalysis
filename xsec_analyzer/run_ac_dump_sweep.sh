#!/usr/bin/env bash
# run_ac_dump_sweep.sh — re-unfold every inclusive config/observable so the additional-
# smearing matrix A_C is dumped into the closure file (UnfolderNuMI now writes h_A_C
# there). Without A_C on disk there is no way to check how much shape information the
# regularisation kept: an A_C whose rows are all the same vector is an outer product
# 1.v^T, which maps EVERY model onto one curve, so agreement with the data is by
# construction rather than physical. FHC theta_mu turned out to be exactly that.
# theta_mu is already done separately; this covers the remaining five observables.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
OBS=(pmu ppi costhmu costhpi thmupi)
CFG=(fhc5:FHC5 rhcfull:RHCFULL comb:COMB)
rm -f $LOG/ac_dump.status

one(){
  local cfg=$1 tag=$2 o=$3
  local xs=configs/ccpi_xsec_config_numi_${o}_${cfg}.txt
  local sl=configs/ccpi_${o}_slice_config_opt.txt
  local out=$PROC/xsec_${tag}_${o}.root
  [ -f "$xs" ] && [ -f "$sl" ] || { echo "  !! missing config $tag $o"; return 0; }
  XSEC_FORCE_REBUILD=1 nice -n 12 bin/UnfolderNuMI "$xs" "$sl" "$out" \
    > $LOG/acdump_${tag}_${o}.log 2>&1 \
    && echo "  OK $tag $o" || echo "  !! FAILED $tag $o"
  echo "DONE_${tag}_${o}" >> $LOG/ac_dump.status
}

echo "==== A_C dump sweep START $(date) NPAR=$NPAR ===="
for c in "${CFG[@]}"; do
  for o in "${OBS[@]}"; do
    one "${c%%:*}" "${c##*:}" "$o" &
    while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done
  done
done
wait
echo "ALL_DONE" >> $LOG/ac_dump.status
echo "==== A_C dump sweep DONE $(date) ===="
