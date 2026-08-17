set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000 >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4
echo "=== convert partial g_numu ==="
gntpc -i g_numu -f gst -o g_numu.gst.root 2>&1 | tail -2
root -l -b -q -e 'TFile f("g_numu.gst.root"); TTree*t=(TTree*)f.Get("gst"); printf("gst entries=%lld\n", t?t->GetEntries():-1);' 2>&1 | grep entries
