// throw_cv_comb.C(seed) — Poisson-thrown CV fake data for the combined FHC+RHC
// (nu_mu + nubar_mu) measurement, at the combined data POT (FHC 6.626e20 + RHC
// 6.079e20 = 1.2705e21). Reads the 8-file combined overlay (3 FHC + 5 RHC).
void throw_cv_comb(int seed){
  const char* PROC="/data/uboone/processed/";
  const char* fhc[]={"Run1_fhc","Run2_fhc","Run4_fhc"};
  const char* rhc[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
  const double SIGPOT=1.7600077e22, DATAPOT=1.2705e21, POTSCALE=DATAPOT/SIGPOT;
  gRandom->SetSeed(seed);
  TString outname=Form("%sxsec-ana-fakedata_CV_comb.root",PROC);
  TChain cin("stv_tree");
  for(auto r:fhc) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,r));
  for(auto r:rhc) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",PROC,r));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; bool ismc;
  cin.SetBranchAddress("tuned_cv_weight",&tcv);cin.SetBranchAddress("ppfx_cv_weight",&pcv);
  cin.SetBranchAddress("normalisation_weight",&nw);cin.SetBranchAddress("is_mc",&ismc);
  TFile* out=new TFile(outname,"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0);
  float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one);
  ot->SetBranchAddress("ppfx_cv_weight",&one);
  ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i);
    double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0)continue;
    int nc=gRandom->Poisson(cv*POTSCALE);
    for(int c=0;c<nc;c++){ one=1.0f; ot->Fill(); kept++; } }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)DATAPOT); sp.Write("summed_pot",TObject::kOverwrite);
  out->Close();
  printf("COMB seed %d: kept %ld (POTSCALE=%.5f, from %lld MC)\n",seed,kept,POTSCALE,N);
}
