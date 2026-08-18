void flipcheck() {
  TChain ch("stv_tree");
  for (auto r : {"Run1","Run2","Run4"})
    ch.Add(Form("/data/uboone/processed/xsec-ana-%s_fhc_new_numi_flux_fhc_pandora_ntuple.root", r));
  const char* S="CC1mu1piXp";
  TString sel = Form("%s_Selected && %s_MC_Signal", S, S);
  TString cr  = Form("%s_candidate_muon_costh_reco", S);
  TString ct  = Form("%s_candidate_muon_costh_true", S);

  long n    = ch.GetEntries(sel);
  long back = ch.GetEntries(Form("%s && %s < -0.9", sel.Data(), cr.Data()));
  long bk5  = ch.GetEntries(Form("%s && %s < -0.5", sel.Data(), cr.Data()));
  // "flipped": reco and true directions on opposite sides, badly so
  long flip = ch.GetEntries(Form("%s && %s < -0.5 && %s > 0.5", sel.Data(), cr.Data(), ct.Data()));
  long flip2= ch.GetEntries(Form("%s && %s < 0 && %s > 0", sel.Data(), cr.Data(), ct.Data()));
  long trueback = ch.GetEntries(Form("%s && %s < -0.5", sel.Data(), ct.Data()));

  printf("\n  CC1mu1piXp selected signal, FHC Runs 1/2/4: %ld events\n", n);
  printf("    reco cos(theta_mu) < -0.9        : %5ld  (%.2f%%)\n", back, 100.*back/n);
  printf("    reco cos(theta_mu) < -0.5        : %5ld  (%.2f%%)\n", bk5,  100.*bk5/n);
  printf("    TRUE cos(theta_mu) < -0.5        : %5ld  (%.2f%%)\n", trueback, 100.*trueback/n);
  printf("    hard flip (reco<-0.5, true>+0.5) : %5ld  (%.3f%%)\n", flip, 100.*flip/n);
  printf("    sign flip (reco<0,   true>0)     : %5ld  (%.2f%%)\n", flip2, 100.*flip2/n);
}
