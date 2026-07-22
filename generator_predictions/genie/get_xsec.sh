set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000 >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
echo "=== extract numu + Ar-40 total CC xsec spline -> ROOT ==="
gspl2root -f "$GENIEXSECFILE" -p 14 -t 1000180400 -o splines_numu_ar.root 2>&1 | tail -2
echo "=== convolve tot_cc with NuMI flux -> flux-avg sigma ==="
root -l -b -q -e '
TFile fs("splines_numu_ar.root");
TGraph* g=(TGraph*)fs.Get("nu_mu_Ar40/tot_cc");
if(!g){ fs.ls(); printf("tot_cc graph not found; check names above\n"); return; }
TFile ff("/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/flux/uboone_numi_flux_histograms.root");
TH1D* h=(TH1D*)ff.Get("h_numu");
double num=0, den=0;
for(int i=1;i<=h->GetNbinsX();++i){ double E=h->GetBinCenter(i), phi=h->GetBinContent(i)*h->GetBinWidth(i);
  double sig=g->Eval(E); num+=sig*phi; den+=phi; }
double favg=num/den;
printf("GENIE numu CC flux-avg sigma = %.4e (spline units, 1e-38 cm^2)\n", favg);
printf("  sigma at 1 GeV = %.4e ; at 0.5 GeV = %.4e\n", g->Eval(1.0), g->Eval(0.5));
printf("  per-nucleus if ~1e1 (x1e-38=1e-37), per-nucleon if ~1e-1 at these E\n");
' 2>&1 | grep -iE "sigma|per-nucleus|not found|tot_cc|nu_mu|KEY"
