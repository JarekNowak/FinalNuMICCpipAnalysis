set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e '
TFile fs("splines_numu_ar.root");
TDirectory* d=(TDirectory*)fs.Get("nu_mu_Ar40");
// collect all CC graphs (name contains "cc", excludes "nc")
std::vector<TGraph*> ccg;
TIter it(d->GetListOfKeys()); TKey* k;
while((k=(TKey*)it())){ TString n=k->GetName();
  if(n.Contains("cc") && !n.Contains("nc")){ TGraph* g=(TGraph*)d->Get(n); if(g) ccg.push_back(g); } }
printf("summed %zu CC channel graphs\n", ccg.size());
auto sigTotCC=[&](double E){ double s=0; for(auto g:ccg){ double v=g->Eval(E); if(v>0) s+=v; } return s; };
printf("sig_tot_cc: 0.5GeV=%.4e 1GeV=%.4e 2GeV=%.4e [spline units, 1e-38 cm2 per Ar]\n", sigTotCC(0.5),sigTotCC(1.0),sigTotCC(2.0));
TFile ff("/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/flux/uboone_numi_flux_histograms.root");
TH1D* h=(TH1D*)ff.Get("h_numu");
double num=0,den=0;
for(int i=1;i<=h->GetNbinsX();++i){ double E=h->GetBinCenter(i); if(E<0.02) continue; double phi=h->GetBinContent(i)*h->GetBinWidth(i);
  num+=sigTotCC(E)*phi; den+=phi; }
printf("GENIE numu flux-avg tot_CC = %.4e [spline units] per Ar nucleus\n", num/den);
printf("  (GENIE splines are per-nucleus; expect ~1e1 in 1e-38 units => ~1e-37 cm2/Ar for total CC)\n");
' 2>&1 | grep -iE "summed|sig_tot|flux-avg|per Ar|spline"
