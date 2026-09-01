// throw_perrun_fhc.C(seed) — per-RUN Poisson CV fake data for FHC, so the fake data
// is structured like the (per-run) beam-on data. Each run i is thrown at its own
// data POT: POTSCALE_i = D_i / MC_POT_i(native). Writes one fake-data file per run,
// carrying summed_pot = D_i (matches the onBNB convention; framework reads onBNB POT
// from file_properties, but the parameter must be present). Ready for real data.
void throw_one_ens(const char* infile, const char* outfile, double dpot, double mcpot, int seed){
  const char* P="/data/uboone/processed/";
  double potscale = dpot/mcpot;
  gRandom->SetSeed(seed);
  TChain cin("stv_tree"); cin.Add(Form("%s%s",P,infile));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; cin.SetBranchAddress("tuned_cv_weight",&tcv);
  cin.SetBranchAddress("ppfx_cv_weight",&pcv); cin.SetBranchAddress("normalisation_weight",&nw);
  TFile* out=new TFile(Form("%s%s",P,outfile),"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one); ot->SetBranchAddress("ppfx_cv_weight",&one); ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i); double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0)continue;
    int nc=gRandom->Poisson(cv*potscale); for(int c=0;c<nc;c++){one=1.0f;ot->Fill();kept++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)dpot); sp.Write("summed_pot",TObject::kOverwrite);
  out->Close();
  printf("  %-55s POTSCALE=%.5f kept=%ld\n",outfile,potscale,kept);
}
void throw_ensemble_fhc(int throw_id=1){
  // Ensemble member: identical statistics to throw_perrun_fhc, writing to throw-tagged
  // names. The fixed names are the fake data the released result is built on; writing
  // there per throw would destroy the release inputs and serialise the ensemble.
  // throw_one_ens above is the original throw_one verbatim -- do not reimplement it.
  int s = 1000*throw_id;
  printf("FHC ensemble member %d (seed base %d):\n", throw_id, s);
  throw_one_ens("xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    Form("ens/fakedata_fhc_run1_t%d.root",throw_id), 3.283e20, 2.3282e21, s+0);
  throw_one_ens("xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    Form("ens/fakedata_fhc_run2_t%d.root",throw_id), 1.268e20, 2.4934e21, s+1);
  throw_one_ens("xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    Form("ens/fakedata_fhc_run4_t%d.root",throw_id), 2.075e20, 2.8335e21, s+2);
  throw_one_ens("xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root",
    Form("ens/fakedata_fhc_run5_t%d.root",throw_id), 2.231e20, 1.9300e21, s+3);
}
