#!/bin/bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
[[ -z "$NUWRO_FQ_DIR" ]] && { echo "NUWRO_SETUP_FAIL"; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4
rm -f nuwro_cc1pi
echo "[compile]"; g++ ../analysis/nuwro_cc1pi.cc -o nuwro_cc1pi -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" 2>&1 | head -15
[[ -x nuwro_cc1pi ]] || { echo "NO_BINARY_AFTER_COMPILE"; exit 1; }
echo "[compiled OK, binary $(ls -la --time-style=+%H:%M nuwro_cc1pi | awk '{print $6}')]"
./nuwro_cc1pi out_numu.root    pred_numu.root    2>/dev/null | grep -i signal
./nuwro_cc1pi out_numubar.root pred_numubar.root 2>/dev/null | grep -i signal
root -l -b -q "combine_newg4.C(\"pred_numu.root\",\"pred_numubar.root\",\"nuwro_newg4_final.root\")" 2>/dev/null | grep -i "ppi"
root -l -b -q "make_fte.C(\"nuwro_newg4_final.root\",\"nuwro_newg4_fte.root\")" 2>/dev/null | grep -i wrote
echo NUWRO_REANA_DONE
