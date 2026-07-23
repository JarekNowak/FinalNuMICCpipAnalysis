HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/neut
cd "$HERE"; echo "[neut] start"
# ROOT 5.34.36 runtime (what neutroot2 + neutclass were built against)
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 2>/dev/null
export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib
export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
export LD_LIBRARY_PATH="/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib:$LD_LIBRARY_PATH"
NEUT=/home/t2k/nowak/generators/neut/src/neutsmpl/neutroot2
NC=/home/t2k/nowak/generators/neut/src/neutclass
echo "root $(root-config --version)"

# --- generate numu + numubar (flux-driven; card uses flux.root in this dir) ---
[[ -s neutvect_numu.root ]]    || { echo "[gen] numu";    $NEUT card_numu.card    neutvect_numu.root    > gen_numu.log 2>&1; }
[[ -s neutvect_numubar.root ]] || { echo "[gen] numubar"; $NEUT card_numubar.card neutvect_numubar.root > gen_numubar.log 2>&1; }

# --- analyze (ACLiC: compile with the neutclass headers/libs) ---
cat > _run_ana.C <<AEOF
{
  gSystem->AddIncludePath("-I$NC");
  gSystem->Load("$NC/neutvtx.so"); gSystem->Load("$NC/neutpart.so");
  gSystem->Load("$NC/neutfsipart.so"); gSystem->Load("$NC/neutfsivert.so");
  gSystem->Load("$NC/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi.C+");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numu.root\",\"pred_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numubar.root\",\"pred_numubar.root\")");
}
AEOF
root -l -b -q _run_ana.C 2>&1 | grep -iE "NEUT:|error"

# --- combine numu+numubar with authoritative flux fractions ---
root -l -b -q "../combine_newg4.C(\"pred_numu.root\",\"pred_numubar.root\",\"neut_newg4_final.root\")" 2>&1 | grep -i integral
echo "[neut] done"
