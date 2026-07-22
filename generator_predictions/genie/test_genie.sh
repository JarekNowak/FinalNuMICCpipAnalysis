set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000
setup genie_phyopt v3_04_00 -q dkcharmtau
set -e
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
FLUX=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/flux/uboone_numi_flux_histograms.root
echo "=== small test: 5k numu on Ar-40, NuMI flux ==="
gevgen -n 5000 -p 14 -t 1000180400 -e 0.05,20 \
  -f "$FLUX,h_numu" --cross-sections "$GENIEXSECFILE" \
  --tune G18_10a_02_11a --seed 1 -o test_numu_ar 2>&1 | grep -iE "cross section|flux.avg|flux.weight|scale|nu_mu|total" | tail -8
echo "=== convert to gst ==="
gntpc -i test_numu_ar.0.ghep.root -f gst -o test_numu_ar.gst.root 2>&1 | tail -2
echo "=== gst branches ==="
root -l -b -q -e 'TFile f("test_numu_ar.gst.root"); TTree*t=(TTree*)f.Get("gst"); t->Print();' 2>&1 | grep -iE "neu|cc|nfpip|nfpim|nfpi0|nf |pdgf|pxf|Ef|nfp|wght|Ev" | head -20
