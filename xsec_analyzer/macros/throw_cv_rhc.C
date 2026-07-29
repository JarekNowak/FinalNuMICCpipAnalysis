// throw_cv_rhc.C(seed) — Poisson-thrown CV fake data for the RHC (nubar-dominant)
// measurement, at the RHC data POT (Run1+2+4 RHC ~ 6.079e20, excl Run3). Same
// mechanism as throw_cv.C: each MC event -> Poisson(cv*POTSCALE) weight-1 copies.
void throw_cv_rhc(int seed){
  const char* PROC="/data/uboone/processed/";
  const char* runs[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; const int NR=5;
  const double SIGPOT=9.944763e21, DATAPOT=6.079e20, POTSCALE=DATAPOT/SIGPOT;
  gRandom->SetSeed(seed);
  TString outname=Form("%sxsec-ana-fakedata_CV_rhc.root",PROC);
  TChain cin("stv_tree");
  for(int r=0;r<NR;r++) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",PROC,runs[r]));
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
  printf("RHC seed %d: kept %ld (POTSCALE=%.5f)\n",seed,kept,POTSCALE);
}
