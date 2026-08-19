// w_resolution_decomp.C -- what drives the W_pipr resolution. Answer: the pion momentum
// (r = +0.726) and nothing else. RMS 0.142 -> 0.061 if the pion is well reconstructed.
// What drives the W_pipr resolution? Correlate its residual against each input's residual,
// and measure how the W spread collapses when a given input happens to be well reconstructed.
void wdecomp() {
  TChain ch("stv_tree");
  for (auto r : {"Run1","Run2","Run4"})
    ch.Add(Form("/data/uboone/processed/w/xsec-ana-%s_fhc_new_numi_flux_fhc_pandora_ntuple.root", r));
  const char* S="CC1mu1pi1p";
  TString sel=Form("%s_Selected && %s_MC_Signal",S,S);
  // pion reco momentum as the ANALYSIS uses it: mu-hypothesis range momentum, mass-corrected
  TString ppi=Form("sqrt(pow(sqrt(pow(%s_candidate_pion_mom_reco,2)+0.011164)-0.10566+0.13957,2)-0.019480)",S);
  TString dW  =Form("(%s_W_pipr_reco-%s_W_pipr_true)",S,S);
  TString dpi =Form("(%s-%s_candidate_pion_mom_true)",ppi.Data(),S);
  TString dpr =Form("(%s_proton_mom_reco-%s_proton_mom_true)",S,S);
  TString dcpi=Form("(%s_candidate_pion_costh_reco-%s_candidate_pion_costh_true)",S,S);
  TString dcpr=Form("(%s_proton_costh_reco-%s_proton_costh_true)",S,S);

  struct V{const char* nm; TString ex;};
  V vs[4]={{"pion momentum",dpi},{"proton momentum",dpr},
           {"pion cos(theta)",dcpi},{"proton cos(theta)",dcpr}};
  printf("\n  correlation of the W_pipr residual with each input residual\n");
  for (auto& v : vs) {
    ch.Draw(Form("%s:%s>>h(100,-1,1,100,-1,1)",dW.Data(),v.ex.Data()), sel, "goff");
    TH2D* h=(TH2D*)gDirectory->Get("h");
    printf("    %-18s  r = %+.3f     input RMS %.3f\n", v.nm, h->GetCorrelationFactor(),
           h->ProjectionX()->GetRMS());
    delete h;
  }
  printf("\n  W_pipr residual RMS when one input is well reconstructed\n");
  ch.Draw(Form("%s>>h0(200,-1,1)",dW.Data()), sel, "goff");
  printf("    all events                      %.3f\n", ((TH1D*)gDirectory->Get("h0"))->GetRMS());
  struct C{const char* nm; TString cut;};
  C cs[3]={{"|dp_pi| < 0.05 GeV/c",  Form("fabs%s<0.05",dpi.Data())},
           {"|dp_p|  < 0.05 GeV/c",  Form("fabs%s<0.05",dpr.Data())},
           {"both momenta < 0.05",   Form("fabs%s<0.05 && fabs%s<0.05",dpi.Data(),dpr.Data())}};
  for (auto& c : cs) {
    ch.Draw(Form("%s>>h1(200,-1,1)",dW.Data()), sel+" && "+c.cut, "goff");
    TH1D* h=(TH1D*)gDirectory->Get("h1");
    printf("    %-30s  %.3f   (%.0f events)\n", c.nm, h->GetRMS(), h->GetEntries());
    delete h;
  }
}
