// w_resolution.C -- W_pipr and W_had residual resolution against their bin widths.
// W_pipr: RMS 0.142, robust half-width 0.101, bias -0.111 GeV/c^2, versus 40 MeV/c^2 bins.
// W_pipr resolution against its bin widths: is the binning finer than the resolution?
void w_res() {
  TChain ch("stv_tree");
  for (auto r : {"Run1","Run2","Run4"})
    ch.Add(Form("/data/uboone/processed/w/xsec-ana-%s_fhc_new_numi_flux_fhc_pandora_ntuple.root", r));
  const char* S="CC1mu1pi1p";
  TString sel=Form("%s_Selected && %s_MC_Signal",S,S);
  for (const char* v : {"W_pipr","W_had"}) {
    TString res=Form("%s_%s_reco - %s_%s_true",S,v,S,v);
    ch.Draw(Form("%s>>h(200,-1.0,1.0)",res.Data()), sel, "goff");
    TH1D* h=(TH1D*)gDirectory->Get("h");
    printf("\n  %s residual (reco - true): mean %+.3f  RMS %.3f GeV/c^2  (%.0f events)\n",
           v, h->GetMean(), h->GetRMS(), h->GetEntries());
    // quantiles for a robust width
    double q[3]={0.16,0.5,0.84}, x[3];
    h->GetQuantiles(3,x,q);
    printf("    16/50/84%%: %+.3f %+.3f %+.3f   -> half-width %.3f\n", x[0],x[1],x[2],(x[2]-x[0])/2);
    delete h;
  }
  printf("\n  W_pipr analysis bin widths: 0.110 0.040 0.040 0.070 0.130 1.430 GeV/c^2\n");
  printf("  W_had  analysis bin widths: 0.780 0.300 0.130 0.100 0.130 1.300 GeV/c^2\n");
}
