#!/bin/bash
# run_nuwro_neut.sh <nuwro|neut> -- the NuWro and NEUT readers of the two-pion generator predictions on the samples
# of generator_predictions/newg4, as generator_predictions/gen2d/reana_{nuwro,neut}_ext.sh:
#   apptainer exec -B /cvmfs -B /home -B /data /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest \
#       bash scripts/mp2pi_gen/run_nuwro_neut.sh <nuwro|neut>
G=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4; O=/data/uboone/processed/mp2pi_2f/gen
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); W=$(mktemp -d); cd "$W"; mkdir -p "$O"
if [ "$1" = nuwro ]; then
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
  setup nuwro v21_09_1 -q e20:prof; export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
  g++ -O2 "$D/nuwro_2pi.cc" -o nuwro_2pi -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so" || exit 1
  for fl in numu numubar; do ./nuwro_2pi $G/out_$fl.root $O/gen2pi_nuwro_$fl.txt 2>&1 | grep -v "AutoLoading\|Missing FileEntry\|requested to autoload"; done
else
  source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
  export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib; export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
  export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
  N=/home/t2k/nowak/generators/neut/src/neutclass; cp "$D/neut_2pi.C" .
  cat > run.C <<MACRO
{ gSystem->AddIncludePath("-I$N");
  gSystem->Load("$N/neutvtx.so"); gSystem->Load("$N/neutpart.so"); gSystem->Load("$N/neutfsipart.so");
  gSystem->Load("$N/neutfsivert.so"); gSystem->Load("$N/neutvect.so");
  gROOT->ProcessLine(".L neut_2pi.C+");
  gROOT->ProcessLine("neut_2pi(\\"$G/neut/neutvect_numu.root\\",\\"$O/gen2pi_neut_numu.txt\\")");
  gROOT->ProcessLine("neut_2pi(\\"$G/neut/neutvect_numubar.root\\",\\"$O/gen2pi_neut_numubar.txt\\")"); }
MACRO
  root -l -b -q run.C 2>&1 | grep -iE "events|error"
fi
command rm -rf -- "$W"
