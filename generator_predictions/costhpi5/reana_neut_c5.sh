#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib; export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/costhpi5
rm -f neut_cc1pi_c5_C.so neut_cc1pi_c5_C_ACLiC_dict_rdict.pcm
root -l -b -q _reana_c5.C 2>&1 | grep -iE "signal|error:|undefined|wrote" | head
echo NEUT_C5_DONE
