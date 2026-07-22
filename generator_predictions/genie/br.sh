set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e 'TFile f("test_numu_ar.gst.root");TTree*t=(TTree*)f.Get("gst");
for(auto b:{"pdgl","pxl","pyl","pzl","El","pdgf","pxf","nfpip","nfpim","nfpi0"}) printf("  %-6s %s\n",b,t->GetBranch(b)?t->GetBranch(b)->GetTitle():"MISSING");
' 2>&1 | grep -E "pdgl|pxl|pyl|pzl|El|pdgf|pxf|nfpi"
