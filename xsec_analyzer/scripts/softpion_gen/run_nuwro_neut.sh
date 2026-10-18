#!/bin/bash
# run_nuwro_neut.sh <nuwro|neut> -- the NuWro and NEUT readers of the soft-pion study on the samples of
# generator_predictions/newg4, as generator_predictions/gen2d/reana_{nuwro,neut}_ext.sh:
#   apptainer exec -B /cvmfs -B /home -B /data /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest \
#       bash scripts/softpion_gen/run_nuwro_neut.sh <nuwro|neut>
G=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4; O=/data/uboone/processed/mp2pi_2f/softpion
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); W=$(mktemp -d); cd "$W"
if [ "$1" = nuwro ]; then
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
  setup nuwro v21_09_1 -q e20:prof; export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
  g++ -O2 "$D/nuwro_softpion.cc" -o nuwro_softpion -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" || exit 1
  for fl in numu numubar; do ./nuwro_softpion $G/out_$fl.root $O/gen_nuwro_$fl.json; done
else
  source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
  export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib; export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
  export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
  N=/home/t2k/nowak/generators/neut/src/neutclass; cp "$D/neut_softpion.C" .
  cat > run.C <<MACRO
{ gSystem->AddIncludePath("-I$N");
  gSystem->Load("$N/neutvtx.so"); gSystem->Load("$N/neutpart.so"); gSystem->Load("$N/neutfsipart.so");
  gSystem->Load("$N/neutfsivert.so"); gSystem->Load("$N/neutvect.so");
  gROOT->ProcessLine(".L neut_softpion.C+");
  gROOT->ProcessLine("neut_softpion(\\"$G/neut/neutvect_numu.root\\",\\"$O/gen_neut_numu.json\\")");
  gROOT->ProcessLine("neut_softpion(\\"$G/neut/neutvect_numubar.root\\",\\"$O/gen_neut_numubar.json\\")"); }
MACRO
  root -l -b -q run.C 2>&1 | grep -iE "events|error"
fi
command rm -rf -- "$W"
