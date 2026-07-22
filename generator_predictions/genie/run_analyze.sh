set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
# sigma_tot_cc per Ar (flux-avg numu) = 13.31 [1e-38 cm2] -> pass in cm2: 13.31e-38
root -l -b -q 'analyze_gst.C("test_numu_ar.gst.root",13.31e-38,"out_genie_numu.root")' 2>&1 | grep -iE "CC1pi|sigma"
root -l -b -q -e 'TFile f("out_genie_numu.root"); TH1D*h=(TH1D*)f.Get("costhmu"); printf("GENIE numu CC1pi costhmu integral(width) per-nucleon=%.4e ; per-Ar=%.4e [1e-38]\n",h->Integral("width")/1e-38, 40*h->Integral("width")/1e-38);' 2>&1 | grep -iE "GENIE|integral"
