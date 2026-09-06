#!/bin/bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
[[ -z "$NUWRO_FQ_DIR" ]] && { echo NUWRO_SETUP_FAIL; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/costhpi5
g++ nuwro_cc1pi_c5.cc -o nuwro_cc1pi_c5 -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" 2>&1 | head -15
[[ -x nuwro_cc1pi_c5 ]] || { echo NO_BINARY; exit 1; }
./nuwro_cc1pi_c5 ../newg4/out_numu.root    nuwro_c5_numu.root    2>&1 | grep -i signal
./nuwro_cc1pi_c5 ../newg4/out_numubar.root nuwro_c5_numubar.root 2>&1 | grep -i signal
echo NUWRO_C5_DONE
