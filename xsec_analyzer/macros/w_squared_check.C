// w_squared_check.C -- W vs W^2 with mapped edges gives identical diagonals (83.2/77.9).
// A monotonic reparameterisation cannot change the response.
// Does measuring in W^2 instead of W change the response?
void w2() {
  TChain ch("stv_tree");
  for (auto r : {"Run1","Run2","Run4"})
    ch.Add(Form("/data/uboone/processed/w/xsec-ana-%s_fhc_new_numi_flux_fhc_pandora_ntuple.root", r));
  const char* S="CC1mu1pi1p";
  TString base=Form("%s_Selected && %s_MC_Signal",S,S);
  TString Wt=Form("%s_W_pipr_true",S), Wr=Form("%s_W_pipr_reco",S);
  TString W2t=Form("pow(%s_W_pipr_true,2)",S), W2r=Form("pow(%s_W_pipr_reco,2)",S);

  auto diag=[&](TString tv, TString rv, double cut){
    double n1=ch.GetEntries(Form("%s && %s < %g",base.Data(),tv.Data(),cut));
    double d1=ch.GetEntries(Form("%s && %s < %g && %s < %g",base.Data(),tv.Data(),cut,rv.Data(),cut));
    double n2=ch.GetEntries(Form("%s && %s >= %g",base.Data(),tv.Data(),cut));
    double d2=ch.GetEntries(Form("%s && %s >= %g && %s >= %g",base.Data(),tv.Data(),cut,rv.Data(),cut));
    printf("      diagonals %5.1f%% / %5.1f%%   N %5.0f / %5.0f\n",
           n1?100*d1/n1:0., n2?100*d2/n2:0., n1, n2);
  };
  printf("\n  (a) W, split at 1.15 GeV/c^2   [the adopted scheme]\n");   diag(Wt,Wr,1.15);
  printf("  (b) W^2, split at 1.15^2 = 1.3225 GeV^2/c^4  [same events]\n"); diag(W2t,W2r,1.3225);
  printf("\n  resolution in each variable:\n");
  ch.Draw(Form("%s-%s>>hA(200,-1,1)",Wr.Data(),Wt.Data()), base, "goff");
  TH1D*a=(TH1D*)gDirectory->Get("hA");
  ch.Draw(Form("%s-%s>>hB(200,-3,3)",W2r.Data(),W2t.Data()), base, "goff");
  TH1D*b=(TH1D*)gDirectory->Get("hB");
  printf("    sigma(W)   = %.3f GeV/c^2     bin width at 1.15 split: first bin 0.07\n", a->GetRMS());
  printf("    sigma(W^2) = %.3f GeV^2/c^4   first bin width in W^2:  %.3f\n", b->GetRMS(), 1.3225-1.08*1.08);
  printf("    ratio sigma/binwidth:  W %.2f   W^2 %.2f\n",
         a->GetRMS()/0.07, b->GetRMS()/(1.3225-1.08*1.08));
}
