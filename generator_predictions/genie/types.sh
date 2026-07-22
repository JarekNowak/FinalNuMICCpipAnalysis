set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e 'TFile f("test_numu_ar.gst.root");TTree*t=(TTree*)f.Get("gst");
for(auto b:{"cc","nf","nfpip","nfpim","nfpi0","nfkp","nfkm","nfk0","pdgf"}){TBranch*br=t->GetBranch(b); if(br)printf("  %-7s %s\n",b,br->GetLeaf(b)?br->GetLeaf(b)->GetTypeName():br->GetTitle());}
// print first signal-ish event nf and pdgf
t->Scan("nf:nfpip:nfpim:nfpi0:pdgf","cc&&(nfpip+nfpim)==1&&nfpi0==0","",3);
' 2>&1 | grep -E "cc |nf |nfpi|nfk|pdgf|\*|Bool|Int|Double" | head -20
