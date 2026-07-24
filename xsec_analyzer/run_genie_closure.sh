#!/usr/bin/env bash
# run_genie_closure.sh
# Rebuild the MAIN-chain univmake files with the GENIE fake data (detVar-CV,
# onBNB in configs/file_properties_numi.txt) and re-unfold all five observables
# to regenerate the closure plots.  Mirrors run_nuwro_observables.sh but for the
# default (non-nuwro) config set, and includes pmu.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u

LOGDIR=../logs
mkdir -p "$LOGDIR" unfold_output
OBS=( pmu costhmu costhpi ppi thmupi )

echo "########## STAGE A: univmake (parallel, GENIE fake data) ##########"
pids=()
for o in "${OBS[@]}"; do
  FPM=configs/file_properties_numi.txt \
  BIN_CONFIG=configs/ccpi_${o}_bin_config.txt \
  OUT=/data/uboone/processed/ccpi_Run1_${o}_univmake.root \
  ./run_universe_maker.sh > "$LOGDIR/genie_${o}_univmake.log" 2>&1 &
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
  univ=/data/uboone/processed/ccpi_Run1_${o}_univmake.root
  if [[ ! -s "$univ" ]]; then echo "  [skip] $o: no univ file"; continue; fi
  XSEC=configs/ccpi_xsec_config_numi_${o}.txt \
  SLICE=configs/ccpi_${o}_slice_config.txt \
  OUT=unfold_output/ccpi_Run1_${o}_xsec.root \
  ./run_unfolder.sh > "$LOGDIR/genie_${o}_unfolder.log" 2>&1 \
    && echo "  unfolder[$o] OK -> unfold_output/ccpi_Run1_${o}_xsec.root" \
    || echo "  unfolder[$o] rc=$? (see $LOGDIR/genie_${o}_unfolder.log)"
done

echo "########## DONE (univmake fail=$fail) ##########"
