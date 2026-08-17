#!/usr/bin/env bash
# fix_inclusive_rebuild.sh — REBUILD the inclusive RHC + COMB univmakes so the fake data
# is re-thrown from the per-run-POT-corrected CV (the multi-file POT fix, commit f60547d).
# Root cause: the existing inclusive RHC/COMB univmakes baked in fake data thrown BEFORE
# the fix (Run3's 5 files over-counted ~2.3x for the RHC mixture), so the correct CV/tune
# prediction sits ~2x below the inflated data. FHC is single-file -> immune -> NOT rebuilt.
# After each univmake rebuilds, re-unfold (config already carries the A_C-smeared
# generators). Writes to a .new temp then moves on success so a crash can't corrupt a file.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs; NPAR=${1:-2}
declare -A TAG=( [rhcfull]=RHCFULL [comb]=COMB )
OBS=(pmu ppi costhmu costhpi thmupi)
rm -f $LOG/incl_rebuild.status

univ_one(){  # cfg TAG obs
  local cfg=$1 t=$2 o=$3
  local out=$PROC/ccpi_${t}_${o}_univmake.root
  echo "[$(date +%H:%M)] START univmake $t $o"
  FPM=configs/file_properties_numi_${cfg}.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT="${out}.new" nice -n 15 ./run_universe_maker.sh > $LOG/incl_rebuild_${cfg}_${o}_univ.log 2>&1
  if [ "$(stat -c%s "${out}.new" 2>/dev/null||echo 0)" -gt 30000000 ]; then
    mv -f "${out}.new" "$out"
    echo "[$(date +%H:%M)] END   univmake $t $o ($(($(stat -c%s $out)/1000000))MB)"
  else
    echo "[$(date +%H:%M)] FAIL  univmake $t $o (output too small)"; rm -f "${out}.new"
  fi
}

echo "==== INCLUSIVE RHC+COMB REBUILD START $(date) NPAR=$NPAR ===="
for cfg in rhcfull comb; do
  t=${TAG[$cfg]}
  echo "==== $cfg ($t) univmakes $(date) ===="
  for o in "${OBS[@]}"; do
    univ_one "$cfg" "$t" "$o" &
    while [ "$(jobs -rp|wc -l)" -ge "$NPAR" ]; do wait -n; done
  done
  wait
  echo "==== $cfg ($t) unfolds $(date) ===="
  for o in "${OBS[@]}"; do
    uout=$PROC/ccpi_${t}_${o}_univmake.root
    [ "$(stat -c%s "$uout" 2>/dev/null||echo 0)" -gt 30000000 ] || { echo "  !! $t $o univmake missing"; continue; }
    XSEC_FORCE_REBUILD=1 nice -n 15 bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_${cfg}.txt \
      configs/ccpi_${o}_slice_config_opt.txt $PROC/xsec_${t}_${o}.root \
      > $LOG/incl_rebuild_${cfg}_${o}_unfold.log 2>&1 \
      && echo "  OK $t $o  tune/data=$(grep -c NORMDBG $LOG/incl_rebuild_${cfg}_${o}_unfold.log)dbg" \
      || echo "  !! unfold $t $o FAILED"
    echo "DONE_${cfg}_${o}" >> $LOG/incl_rebuild.status
  done
done
echo "ALL_DONE" >> $LOG/incl_rebuild.status
echo "==== INCLUSIVE RHC+COMB REBUILD DONE $(date) ===="
