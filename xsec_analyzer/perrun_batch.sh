#!/usr/bin/env bash
# perrun_batch.sh — sequential per-run univmake+unfold batch for the remaining
# FHC observables + all RHC + all combined. Runs one univmake at a time (each uses
# all cores; sequential avoids the concurrent-write crash). Dirt summed_pot is
# re-set per config (shared file, config-specific scaling). ~1.5-2h per univmake.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:$LD_LIBRARY_PATH"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed
LOG=../logs
DIRT=$PROC/xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root

setpot_dirt() {  # $1 = summed_pot value
  root.exe -l -b -q -e "TFile*f=TFile::Open(\"$DIRT\",\"update\");TParameter<float> p(\"summed_pot\",(float)$1);p.Write(\"summed_pot\",TObject::kOverwrite);f->Close();" >/dev/null 2>&1
  echo "[dirt] summed_pot = $1"
}

run_one() {  # $1 cfg(lower)  $2 TAG(upper)  $3 obs
  local cfg=$1 TAG=$2 obs=$3
  local uout=$PROC/ccpi_${TAG}_${obs}_univmake.root
  echo "[$(date +%H:%M)] univmake $TAG $obs ..."
  FPM=configs/file_properties_numi_${cfg}.txt BIN_CONFIG=configs/ccpi_${obs}_bin_config_opt.txt \
    OUT=$uout ./run_universe_maker.sh > $LOG/perrun_${cfg}_${obs}_univ.log 2>&1
  local sz=$(stat -c%s "$uout" 2>/dev/null || echo 0)
  if [ "$sz" -lt 50000000 ]; then echo "  !! univmake $TAG $obs FAILED (size $sz) — see log"; return 1; fi
  echo "[$(date +%H:%M)] unfold $TAG $obs (univmake $((sz/1000000))MB) ..."
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${obs}_${cfg}.txt configs/ccpi_${obs}_slice_config_opt.txt \
    $PROC/xsec_perrun_${cfg}_${obs}.root > $LOG/perrun_${cfg}_${obs}_unfold.log 2>&1 \
    && echo "  OK $TAG $obs" || echo "  !! unfold $TAG $obs FAILED"
}

echo "==== PER-RUN BATCH START $(date) ===="
# FHC remaining observables (pmu already done); dirt already 6.2046e20 but re-assert
setpot_dirt 6.2046e20
for o in ppi costhmu costhpi thmupi; do run_one fhc5 FHC5 $o; done
# RHC (all five)
setpot_dirt 9.1429e19
for o in pmu ppi costhmu costhpi thmupi; do run_one rhcfull RHCFULL $o; done
# Combined (all five)
setpot_dirt 2.757e20
for o in pmu ppi costhmu costhpi thmupi; do run_one comb COMB $o; done
echo "==== PER-RUN BATCH DONE $(date) ===="
