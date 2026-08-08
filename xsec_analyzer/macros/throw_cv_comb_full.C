// throw_cv_comb_full.C(seed) — Poisson-thrown CV fake data for the FULL combined
// FHC+RHC measurement. PER-MODE POTSCALE so each mode is thrown at its own data
// exposure (matches the per-mode summed_pot in file_properties_numi_comb.txt):
//   FHC (chain trees 0-3):  POTSCALE = D_FHC/S_FHC = 0.092402
//   RHC (chain trees 4-13): POTSCALE = D_RHC/S_RHC = 0.071666
//   summed_pot(out) = D_comb = 1.99390e21
void throw_cv_comb_full(int seed){
  const char* PROC="/data/uboone/processed/";
  const char* fhc[]={"Run1_fhc","Run2_fhc","Run4_fhc"};              // + Run5 reweightedPPFX
  const char* rhc[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
  const char* run3[]={"aa","ab","ac","ad","ae"};
  const double PS_FHC=0.092402, PS_RHC=0.071666, DATAPOT=1.99390e21;
  gRandom->SetSeed(seed);
  TChain cin("stv_tree");
  for(auto r:fhc) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,r));
  cin.Add(Form("%sxsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root",PROC)); // FHC tree idx 3
  for(auto r:rhc) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",PROC,r));
  for(auto s:run3) cin.Add(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",PROC,s));
  const int NFHC=4; // trees 0..3 are FHC, 4..13 RHC
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; bool ismc;
  cin.SetBranchAddress("tuned_cv_weight",&tcv);cin.SetBranchAddress("ppfx_cv_weight",&pcv);
  cin.SetBranchAddress("normalisation_weight",&nw);cin.SetBranchAddress("is_mc",&ismc);
  TFile* out=new TFile(Form("%sxsec-ana-fakedata_CV_comb_full.root",PROC),"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one);ot->SetBranchAddress("ppfx_cv_weight",&one);ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long keptF=0,keptR=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i);
    double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0)continue;
    double ps = (cin.GetTreeNumber()<NFHC) ? PS_FHC : PS_RHC;
    int nc=gRandom->Poisson(cv*ps);
    for(int c=0;c<nc;c++){one=1.0f;ot->Fill(); if(cin.GetTreeNumber()<NFHC)keptF++; else keptR++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)DATAPOT); sp.Write("summed_pot",TObject::kOverwrite); out->Close();
  printf("COMBFULL seed %d: kept FHC=%ld RHC=%ld total=%ld (ratio F:R=%.3f, PS_FHC=%.5f PS_RHC=%.5f)\n",
         seed,keptF,keptR,keptF+keptR,(double)keptF/keptR,PS_FHC,PS_RHC);
}
