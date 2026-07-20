#!/usr/bin/env bash
# run_nuwro_observables.sh
# Drives the NuWro fake-data closure through univmake + Unfolder + UnfolderNuMI
# for the four "new" differential observables whose branches were added to
# CC1mu1piXp (pion momentum, muon cos-theta, pion cos-theta, mu-pi opening angle).
# Assumes ProcessNTuples has already been re-run so the processed stv_tree files
# carry the new branches.
#
# univmakes are launched in parallel (each uses ~4 threads); the unfolders run
# afterwards (fast).  Per-observable logs go to ../logs/nuwro_<obs>_*.log.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u

LOGDIR=../logs
# The redirects below fail immediately if this does not exist, and the script
# runs without `set -e`, so every stage would be skipped without a clear error.
mkdir -p "$LOGDIR"
OBS=( ppi costhmu costhpi thmupi )

echo "########## STAGE A: univmake (parallel) ##########"
pids=()
for o in "${OBS[@]}"; do
  FPM=configs/file_properties_numi_nuwro.txt \
  BIN_CONFIG=configs/ccpi_${o}_bin_config.txt \
  OUT=/data/uboone/processed/ccpi_Run1_${o}_nuwro_univmake.root \
  ./run_universe_maker.sh > "$LOGDIR/nuwro_${o}_univmake.log" 2>&1 &
  pids+=($!)
  echo "  launched univmake[$o] pid=${pids[-1]}"
done
fail=0
for i in "${!OBS[@]}"; do
  if wait "${pids[$i]}"; then echo "  univmake[${OBS[$i]}] OK"
  else echo "  univmake[${OBS[$i]}] FAILED (rc=$?)"; fail=1; fi
done

echo "########## STAGE B: Unfolder + UnfolderNuMI (serial) ##########"
for o in "${OBS[@]}"; do
  univ=/data/uboone/processed/ccpi_Run1_${o}_nuwro_univmake.root
  if [[ ! -s "$univ" ]]; then echo "  [skip] $o: no univ file"; continue; fi
  XSEC=configs/ccpi_xsec_config_numi_nuwro_${o}.txt \
  SLICE=configs/ccpi_${o}_slice_config.txt \
  OUT=unfold_output/ccpi_Run1_${o}_nuwro_xsec.root \
  ./run_unfolder.sh > "$LOGDIR/nuwro_${o}_unfolder.log" 2>&1 \
    && echo "  unfolder[$o] OK -> unfold_output/ccpi_Run1_${o}_nuwro_xsec.root" \
    || echo "  unfolder[$o] rc=$? (see $LOGDIR/nuwro_${o}_unfolder.log)"
done

echo "########## DONE (univmake fail=$fail) ##########"
