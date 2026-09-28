#!/bin/bash
# reana_neut_1pinc.sh -- NEUT part of regen_1pinc.sh, inside the SL7 container (same environment as reana_neut_ext.sh)
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib; export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d
rm -f neut_1p_inc_C.so
root -l -b -q _reana_1pinc.C 2>&1 | grep -iE "NEUT|signal|error:|undefined|wrote" | head
echo NEUT_1PINC_DONE
