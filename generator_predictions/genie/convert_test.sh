set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000 >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
echo "=== convert test_numu_ar -> gst ==="
gntpc -i test_numu_ar -f gst -o test_numu_ar.gst.root 2>&1 | tail -3
echo "=== gst content ==="
root -l -b -q -e 'TFile f("test_numu_ar.gst.root"); TTree*t=(TTree*)f.Get("gst");
if(!t){printf("no gst\n");return;}
printf("entries=%lld\n",t->GetEntries());
// list key branches
for(auto b:{"neu","cc","nc","Ev","nfp","nfn","nfpip","nfpim","nfpi0","nf","pdgf","pxf","pyf","pzf","Ef","wght","xsec"}){
  printf("  %-8s : %s\n", b, t->GetBranch(b)?"yes":"NO");
}
// signal-ish counts: CC, exactly 1 charged pi (nfpip+nfpim==1), 0 pi0
int n=t->Draw("1","cc==1 && (nfpip+nfpim)==1 && nfpi0==0","goff");
printf("CC 1charged-pi 0pi0 events: %d of %lld\n", n, t->GetEntries());
' 2>&1 | grep -E "entries|:|events|no gst"
