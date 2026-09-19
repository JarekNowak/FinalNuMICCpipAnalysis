// throw_ensemble_rhc.C(throw_id) -- RHC ensemble member: the group throw of throw_perrun_rhc.C,
// writing throw-tagged names under ens/ so the released fake data is never overwritten.
void throw_group_ens(std::vector<const char*> infiles, const char* outfile, double dpot, double mcpot, int seed){
  const char* P="/data/uboone/processed/";
  double potscale = dpot/mcpot; gRandom->SetSeed(seed);
  TChain cin("stv_tree"); for(auto f:infiles) cin.Add(Form("%s%s",P,f));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; cin.SetBranchAddress("tuned_cv_weight",&tcv);
  cin.SetBranchAddress("ppfx_cv_weight",&pcv); cin.SetBranchAddress("normalisation_weight",&nw);
  TFile* out=new TFile(Form("%s%s",P,outfile),"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one); ot->SetBranchAddress("ppfx_cv_weight",&one); ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i); double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0||cv>30) cv=1.0;
    int nc=gRandom->Poisson(cv*potscale); for(int c=0;c<nc;c++){one=1.0f;ot->Fill();kept++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)dpot); sp.Write("summed_pot",TObject::kOverwrite);
  out->Close(); printf("  %-45s POTSCALE=%.5f kept=%ld\n",outfile,potscale,kept);
}
void throw_ensemble_rhc(int throw_id=1){
  int s=1000*throw_id;
  throw_group_ens({"xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/fakedata_rhc_run1_t%d.root",throw_id),0.6053e20,8.9972e20,s+0);
  throw_group_ens({"xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/fakedata_rhc_run2_t%d.root",throw_id),2.591e20,5.7865e21,s+1);
  throw_group_ens({"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root"},Form("ens/fakedata_rhc_run3_t%d.root",throw_id),5.003e20,5.5185e21,s+2);
  throw_group_ens({"xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/fakedata_rhc_run4_t%d.root",throw_id),2.883e20,3.2586e21,s+3);
}
