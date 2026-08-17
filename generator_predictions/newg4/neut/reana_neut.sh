#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib; export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/neut
rm -f neut_cc1pi_C.so neut_cc1pi_C_ACLiC_dict_rdict.pcm
cat > _reana.C <<'CE'
{ gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi.C+");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numu.root\",\"pred_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numubar.root\",\"pred_numubar.root\")"); }
CE
root -l -b -q _reana.C 2>&1 | grep -iE "NEUT:|signal|error:|undefined" | head
root -l -b -q "../combine_newg4.C(\"pred_numu.root\",\"pred_numubar.root\",\"neut_newg4_final.root\")" 2>&1 | grep -i integral
root -l -b -q "../make_fte.C(\"neut_newg4_final.root\",\"../neut_newg4_fte.root\")" 2>&1 | grep -i wrote
echo NEUT_REANA_DONE
