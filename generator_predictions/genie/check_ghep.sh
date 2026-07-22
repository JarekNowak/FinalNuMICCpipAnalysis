set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e '
TFile f("test_numu_ar");
TTree* h=(TTree*)f.Get("gtree");
if(h){ printf("gtree entries=%lld\n", h->GetEntries());
  // the NtpMCTreeHeader / gtree may store run config
}
f.GetListOfKeys()->Print();
' 2>&1 | grep -iE "KEY|gtree|entries|header|xsec"
