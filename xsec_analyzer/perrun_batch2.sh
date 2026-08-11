#!/usr/bin/env bash
# perrun_batch2.sh — per-run univmake+unfold batch, NPAR observables in parallel
# (default 2) to use the idle cores. Parallelism is WITHIN a config only (same dirt
# summed_pot); each config is a barrier (wait for all its univmakes, then re-set the
# dirt for the next config). Idempotent: skips any univmake whose output already
# exists (>50MB). Different observables write different output files, so there is no
# same-file write collision.  usage: ./perrun_batch2.sh [NPAR]
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:$LD_LIBRARY_PATH"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs
DIRT=$PROC/xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root
NPAR=${1:-2}

setpot_dirt() {
  root.exe -l -b -q -e "TFile*f=TFile::Open(\"$DIRT\",\"update\");TParameter<float> p(\"summed_pot\",(float)$1);p.Write(\"summed_pot\",TObject::kOverwrite);f->Close();" >/dev/null 2>&1
  echo "[$(date +%H:%M)] [dirt] summed_pot = $1"
}
univ_one() {  # cfg TAG obs
  local cfg=$1 TAG=$2 obs=$3 uout=$PROC/ccpi_${2}_${3}_univmake.root
  local sz=$(stat -c%s "$uout" 2>/dev/null || echo 0)
  if [ "$sz" -gt 50000000 ]; then echo "[skip] $TAG $obs (done $((sz/1000000))MB)"; return 0; fi
  echo "[$(date +%H:%M)] START univmake $TAG $obs"
  FPM=configs/file_properties_numi_${cfg}.txt BIN_CONFIG=configs/ccpi_${obs}_bin_config_opt.txt \
    OUT=$uout ./run_universe_maker.sh > $LOG/perrun_${cfg}_${obs}_univ.log 2>&1
  echo "[$(date +%H:%M)] END   univmake $TAG $obs ($(($(stat -c%s $uout 2>/dev/null||echo 0)/1000000))MB)"
}
run_config() {  # cfg TAG dirt obs...
  local cfg=$1 TAG=$2 dirt=$3; shift 3
  setpot_dirt "$dirt"
  for obs in "$@"; do
    univ_one "$cfg" "$TAG" "$obs" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done
  done
  wait
  for obs in "$@"; do
    local uout=$PROC/ccpi_${TAG}_${obs}_univmake.root
    [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 50000000 ] || { echo "  !! $TAG $obs univmake missing, skip unfold"; continue; }
    bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${obs}_${cfg}.txt configs/ccpi_${obs}_slice_config_opt.txt \
      $PROC/xsec_perrun_${cfg}_${obs}.root > $LOG/perrun_${cfg}_${obs}_unfold.log 2>&1 \
      && echo "  OK $TAG $obs" || echo "  !! unfold $TAG $obs FAILED"
  done
}
echo "==== 2-WIDE BATCH START $(date) NPAR=$NPAR ===="
run_config fhc5    FHC5    6.2046e20 ppi costhmu costhpi thmupi
run_config rhcfull RHCFULL 9.1429e19 pmu ppi costhmu costhpi thmupi
run_config comb    COMB    2.757e20  pmu ppi costhmu costhpi thmupi
echo "==== 2-WIDE BATCH DONE $(date) ===="
