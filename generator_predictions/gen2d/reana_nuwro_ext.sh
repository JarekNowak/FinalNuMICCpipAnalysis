#!/bin/bash
# gen2d: compile and run the NuWro readers with the obs_ext.h fills (inside the SL7 container)
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
[[ -z "$NUWRO_FQ_DIR" ]] && { echo NUWRO_SETUP_FAIL; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d
for r in incl 1p; do
  g++ nuwro_$r.cc -o nuwro_$r -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" 2>&1 | head -15
  [[ -x nuwro_$r ]] || { echo NO_BINARY_$r; exit 1; }
done
./nuwro_incl ../newg4/out_numu.root    nuwro_ext_numu.root    2>&1 | grep -i signal
./nuwro_incl ../newg4/out_numubar.root nuwro_ext_numubar.root 2>&1 | grep -i signal
./nuwro_1p   ../newg4/out_numu.root    nuwro_1p_ext_numu.root    2>&1 | grep -i signal
./nuwro_1p   ../newg4/out_numubar.root nuwro_1p_ext_numubar.root 2>&1 | grep -i signal
echo NUWRO_EXT_DONE
