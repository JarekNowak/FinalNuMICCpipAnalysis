#!/bin/bash
# reana_nuwro_1pinc.sh -- NuWro part of regen_1pinc.sh: compile and run the proton-tagged reader with the
# muon/pion observables (obs_1pinc.h) inside the SL7 container (same environment as reana_nuwro_ext.sh)
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
[[ -z "$NUWRO_FQ_DIR" ]] && { echo NUWRO_SETUP_FAIL; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d
g++ nuwro_1p_inc.cc -o nuwro_1p_inc -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" 2>&1 | head -15
[[ -x nuwro_1p_inc ]] || { echo NO_BINARY_1p_inc; exit 1; }
./nuwro_1p_inc ../newg4/out_numu.root    nuwro_1p_inc_numu.root    2>&1 | grep -i signal &
./nuwro_1p_inc ../newg4/out_numubar.root nuwro_1p_inc_numubar.root 2>&1 | grep -i signal &
wait
echo NUWRO_1PINC_DONE
