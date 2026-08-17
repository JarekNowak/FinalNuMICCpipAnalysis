source /home/t2k/nowak/generators/setupNuWro.sh >/dev/null 2>&1
export CERN=/cvmfs/sft.cern.ch/lcg/external/cernlib
export CERN_LEVEL=2006a/x86_64-slc6-gcc47-opt
export ROOTSYS=$(root-config --prefix)
GMPLIB=/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/lib
export LD_LIBRARY_PATH="$GMPLIB:$LD_LIBRARY_PATH"
NEUT=/home/t2k/nowak/generators/neut/src/neutsmpl/neutroot2
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/neut
echo "neutroot2: $(ls -la $NEUT 2>/dev/null | awk '{print $5}')"
sed 's/EVCT-NEVT  500000/EVCT-NEVT  2000/' card_numu.card > card_numu_test.card
echo "[neut] test run 2000 events ..."
$NEUT card_numu_test.card neutvect_test.root 2>&1 | tail -8
echo "=== output ==="; ls -la neutvect_test.root 2>/dev/null | awk '{print $5,$9}'
