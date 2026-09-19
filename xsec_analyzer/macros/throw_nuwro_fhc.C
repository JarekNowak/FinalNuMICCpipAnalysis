// throw_nuwro_fhc.C -- INDEPENDENT-GENERATOR fake data: a Poisson throw of the NuWro overlay
// (xsec-ana-numi_nuwro_overlay_pion_ntuples_run1_fhc.root, 6.650e20 POT, Run-1 FHC detector
// simulation), to be unfolded with the nominal GENIE response. This is the test the reviewer
// asked for and the model-distortion closure cannot replace: a different generator, an
// independent detector-level sample, signal and background both from NuWro.
//
// The overlay carries sentinel tune and PPFX central-value weights (-1), which the safe-weight
// rule maps to one, so every event is thrown at unit weight times the exposure scale: the NuWro
// rate is NuWro's own, with the raw flux. The sample is Run-1 only; it is thrown at the FULL
// corrected FHC exposure (7.766e20) so that the statistical power matches the FHC measurement,
// and paired with Run-1 simulation and beam-off only, whose trigger count is scaled in
// proportion. Truth branches pass through, so the framework's FakeData universe is NuWro truth.
//   root -l -b -q 'macros/throw_nuwro_fhc.C(1)'
void throw_nuwro_fhc(int seed = 1) {
  const char* P = "/data/uboone/processed/";
  const double dpot = 7.766e20, mcpot = 6.65024e20;
  double potscale = dpot / mcpot;
  gRandom->SetSeed(seed);
  TChain cin("stv_tree"); cin.Add(Form("%sxsec-ana-numi_nuwro_overlay_pion_ntuples_run1_fhc.root", P));
  cin.SetBranchStatus("*", 1);
  for (auto b : {"weight_All_UBGenie", "weight_ppfx_all", "weight_reint_all"}) cin.SetBranchStatus(b, 0);
  float tcv, pcv, nw;
  cin.SetBranchAddress("tuned_cv_weight", &tcv); cin.SetBranchAddress("ppfx_cv_weight", &pcv);
  cin.SetBranchAddress("normalisation_weight", &nw);
  TFile* out = new TFile(Form("%sxsec-ana-fakedata_nuwro_fhc_run1.root", P), "recreate"); out->SetCompressionLevel(1);
  TTree* ot = cin.CloneTree(0); float one = 1.0f;
  ot->SetBranchAddress("tuned_cv_weight", &one); ot->SetBranchAddress("ppfx_cv_weight", &one); ot->SetBranchAddress("normalisation_weight", &one);
  Long64_t N = cin.GetEntries(); long kept = 0; double sumw = 0;
  for (Long64_t i = 0; i < N; ++i) {
    cin.GetEntry(i);
    double cv = tcv * pcv * nw; if (!std::isfinite(cv) || cv < 0 || cv > 30) cv = 1.0;   // safe rule: sentinels -> 1
    sumw += cv;
    int nc = gRandom->Poisson(cv * potscale); for (int c = 0; c < nc; ++c) { one = 1.0f; ot->Fill(); ++kept; }
  }
  ot->Write("", TObject::kOverwrite);
  TParameter<float> sp("summed_pot", (float)dpot); sp.Write("summed_pot", TObject::kOverwrite);
  out->Close();
  printf("  nuwro fake data: POTSCALE=%.4f  mean weight %.3f  kept=%ld\n", potscale, sumw / N, kept);
}
