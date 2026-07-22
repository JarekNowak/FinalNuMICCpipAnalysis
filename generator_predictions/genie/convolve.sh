set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e '
TFile fs("splines_numu_ar.root");
TDirectory* d=(TDirectory*)fs.Get("nu_mu_Ar40");
if(!d){ printf("no dir\n"); fs.ls(); return; }
printf("=== graphs in nu_mu_Ar40 ===\n"); d->ls();
TGraph* g=(TGraph*)d->Get("tot_cc");
if(!g){ printf("no tot_cc\n"); return; }
TFile ff("/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/flux/uboone_numi_flux_histograms.root");
TH1D* h=(TH1D*)ff.Get("h_numu");
double num=0,den=0;
for(int i=1;i<=h->GetNbinsX();++i){ double E=h->GetBinCenter(i); double phi=h->GetBinContent(i)*h->GetBinWidth(i);
  double sig=g->Eval(E); num+=sig*phi; den+=phi; }
printf("GENIE numu tot_cc: sig(0.5GeV)=%.4e sig(1GeV)=%.4e sig(2GeV)=%.4e [spline units]\n",g->Eval(0.5),g->Eval(1.0),g->Eval(2.0));
printf("GENIE numu flux-avg tot_cc = %.4e [spline units]\n", num/den);
' 2>&1 | grep -iE "tot_cc|flux-avg|sig\(|graphs|no dir|no tot|KEY|TGraph"
