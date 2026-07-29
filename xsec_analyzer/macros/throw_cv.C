// throw_cv.C(seed) — one Poisson-thrown CV pseudo-dataset (no M_A^RES reweight)
// at the nue-convention FHC data POT, for averaging the flux amplification over
// many independent throws. Each MC event -> Poisson(cv*POTSCALE) weight-1 copies.
void throw_cv(int seed){
  const char* PROC="/data/uboone/processed/";
  const char* runs[]={"Run1_fhc","Run2_fhc","Run4_fhc"}; const int NR=3;
  const double SIGPOT=7.655314e21, DATAPOT=6.626e20, POTSCALE=DATAPOT/SIGPOT;
  gRandom->SetSeed(seed);
  TString outname=Form("%sxsec-ana-fakedata_CV_run15.root",PROC);
  TChain cin("stv_tree");
  for(int r=0;r<NR;r++) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,runs[r]));
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
  printf("seed %d: kept %ld\n",seed,kept);
}
