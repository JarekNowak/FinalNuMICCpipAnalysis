set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e '
TFile fs("splines_numu_ar.root");
TGraph* g=(TGraph*)((TDirectory*)fs.Get("nu_mu_Ar40"))->Get("tot_cc");
TFile ff("/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/flux/uboone_numi_flux_histograms.root");
TH1D* h=(TH1D*)ff.Get("h_numu");
double num=0,den=0;
for(int i=1;i<=h->GetNbinsX();++i){ double E=h->GetBinCenter(i); if(E<0.02)continue; double phi=h->GetBinContent(i)*h->GetBinWidth(i);
  num+=g->Eval(E)*phi; den+=phi; }
double favg=num/den;
printf("GENIE numu tot_cc: 0.44GeV=%.3e 1GeV=%.3e [1e-38 cm2/Ar]\n",g->Eval(0.44),g->Eval(1.0));
printf("GENIE numu flux-avg total CC = %.4e [1e-38 cm2/Ar] = %.4e per nucleon\n", favg, favg/40.);
printf("  (NuWro flux-avg total CC per nucleon was 0.30e-38; per Ar 12e-38)\n");
' 2>&1 | grep -iE "tot_cc|flux-avg|NuWro"
