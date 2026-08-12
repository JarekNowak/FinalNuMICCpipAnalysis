// throw_perrun_sb.C(seed) — per-RUN Poisson CV fake data for the SIDEBAND study,
// reading the sb/-reprocessed MC (which carries the CC1mu1piXp_sb_* flags) and writing
// the fake data into the SAME sb/ dir. Identical throw logic + POT values as
// throw_perrun_fhc.C / throw_perrun_rhc.C; CloneTree(0) carries the sb_ branches
// through so the fake "data" also has the sideband flags. Separate dir → the running
// comb batch's standard fake-data files are untouched.
static const char* SB="/data/uboone/processed/sb/";

void throw_group_sb(std::vector<const char*> infiles, const char* outfile, double dpot, double mcpot, int seed){
  double potscale = dpot/mcpot;
  gRandom->SetSeed(seed);
  TChain cin("stv_tree"); for(auto f:infiles) cin.Add(Form("%s%s",SB,f));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; cin.SetBranchAddress("tuned_cv_weight",&tcv);
  cin.SetBranchAddress("ppfx_cv_weight",&pcv); cin.SetBranchAddress("normalisation_weight",&nw);
  TFile* out=new TFile(Form("%s%s",SB,outfile),"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one); ot->SetBranchAddress("ppfx_cv_weight",&one); ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i); double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0)continue;
    int nc=gRandom->Poisson(cv*potscale); for(int c=0;c<nc;c++){one=1.0f;ot->Fill();kept++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)dpot); sp.Write("summed_pot",TObject::kOverwrite);
  out->Close();
  printf("  %-45s POTSCALE=%.5f kept=%ld\n",outfile,potscale,kept);
}

void throw_perrun_sb(int seed=1){
  printf("SIDEBAND per-run fake data (seed %d):\n",seed);
  // FHC (single-file runs)
  throw_group_sb({"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root"},"xsec-ana-fakedata_fhc_run1.root",3.283e20,2.3282e21,seed);
  throw_group_sb({"xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root"},"xsec-ana-fakedata_fhc_run2.root",1.268e20,2.4934e21,seed+1);
  throw_group_sb({"xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"},"xsec-ana-fakedata_fhc_run4.root",2.075e20,2.8335e21,seed+2);
  throw_group_sb({"xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root"},"xsec-ana-fakedata_fhc_run5.root",2.231e20,1.9300e21,seed+3);
  // RHC (Run1/2 single, Run3/4 grouped)
  throw_group_sb({"xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root"},"xsec-ana-fakedata_rhc_run1.root",0.6053e20,8.9972e20,seed+4);
  throw_group_sb({"xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root"},"xsec-ana-fakedata_rhc_run2.root",2.591e20,5.7865e21,seed+5);
  throw_group_sb({"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root"},"xsec-ana-fakedata_rhc_run3.root",5.003e20,5.5185e21,seed+6);
  throw_group_sb({"xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root"},"xsec-ana-fakedata_rhc_run4.root",2.883e20,3.2586e21,seed+7);
}
