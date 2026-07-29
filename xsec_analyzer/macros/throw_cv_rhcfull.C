void throw_cv_rhcfull(int seed){
  const char* P="/data/uboone/processed/";
  const char* files[]={
    "xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root",
    "xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root",
    "xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root",
    "xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root",
    "xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root",
    "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root",
    "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root",
    "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root",
    "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root",
    "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root"};
  const double SIGPOT=1.546330e22, DATAPOT=1.1082e21, POTSCALE=DATAPOT/SIGPOT;
  gRandom->SetSeed(seed);
  TChain cin("stv_tree"); for(auto f:files) cin.Add(Form("%s%s",P,f));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; bool ismc;
  cin.SetBranchAddress("tuned_cv_weight",&tcv);cin.SetBranchAddress("ppfx_cv_weight",&pcv);
  cin.SetBranchAddress("normalisation_weight",&nw);cin.SetBranchAddress("is_mc",&ismc);
  TFile* out=new TFile(Form("%sxsec-ana-fakedata_CV_rhcfull.root",P),"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one);ot->SetBranchAddress("ppfx_cv_weight",&one);ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i); double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0)continue;
    int nc=gRandom->Poisson(cv*POTSCALE); for(int c=0;c<nc;c++){one=1.0f;ot->Fill();kept++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)DATAPOT); sp.Write("summed_pot",TObject::kOverwrite); out->Close();
  printf("RHCFULL seed %d: kept %ld (POTSCALE=%.5f)\n",seed,kept,POTSCALE);
}
