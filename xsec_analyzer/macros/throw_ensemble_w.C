// throw_ensemble_w.C(throw_id, mode) -- proton-tagged (w/ tree) ensemble member: the group throw of
// throw_ensemble_{fhc,rhc}.C applied to /data/uboone/processed/w/, writing throw-tagged names under
// w/ens/ so the released CC1mu1pi1p fake data is never overwritten. Same runs, same data POT, same MC
// POT and the same seeds as the inclusive ensembles (seed base 1000*throw_id + run index).
//   root -l -b -q 'macros/throw_ensemble_w.C(7,"fhc")'   (mode = "fhc" | "rhc")
#include "throw_guard.h"   // create-once throw files (2026-09-26)
void throw_group_ens_w(std::vector<const char*> infiles, const char* outfile, double dpot, double mcpot, int seed){
  const char* P = gSystem->Getenv("THROW_W_DIR") ? gSystem->Getenv("THROW_W_DIR") : "/data/uboone/processed/w/";
  double potscale = dpot/mcpot; gRandom->SetSeed(seed);
  const TString fin=Form("%s%s",P,outfile); const std::vector<TString> inp=[&]{ std::vector<TString> v; for(auto f:infiles) v.push_back(Form("%s%s",P,f)); return v; }();
  if(throw_reusable(fin,inp,seed,potscale)){ printf("  %-45s POTSCALE=%.5f reused\n",outfile,potscale); return; }
  TChain cin("stv_tree"); for(auto f:infiles) cin.Add(Form("%s%s",P,f));
  cin.SetBranchStatus("*",1);
  for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);
  float tcv,pcv,nw; cin.SetBranchAddress("tuned_cv_weight",&tcv);
  cin.SetBranchAddress("ppfx_cv_weight",&pcv); cin.SetBranchAddress("normalisation_weight",&nw);
  // private temporary name + rename: concurrent ensembles sharing a throw id re-throw the same file
  // (the 2026-09-20 race); a reader must only ever see a complete file.
  TString ftmp=Form("%s.tmp%d",fin.Data(),gSystem->GetPid());
  TFile* out=new TFile(ftmp,"recreate"); out->SetCompressionLevel(1);
  TTree* ot=cin.CloneTree(0); float one=1.0f;
  ot->SetBranchAddress("tuned_cv_weight",&one); ot->SetBranchAddress("ppfx_cv_weight",&one); ot->SetBranchAddress("normalisation_weight",&one);
  Long64_t N=cin.GetEntries(); long kept=0;
  for(Long64_t i=0;i<N;i++){ cin.GetEntry(i); double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0||cv>30) cv=1.0; // safe-weight rule
    int nc=gRandom->Poisson(cv*potscale); for(int c=0;c<nc;c++){one=1.0f;ot->Fill();kept++;} }
  ot->Write("",TObject::kOverwrite);
  TParameter<float> sp("summed_pot",(float)dpot); sp.Write("summed_pot",TObject::kOverwrite);
  throw_stamp(seed,potscale); out->Close(); throw_install(ftmp,fin,inp,seed,potscale);
  printf("  %-45s POTSCALE=%.5f kept=%ld\n",outfile,potscale,kept);
}
void throw_ensemble_w(int throw_id=1, const char* mode="fhc"){
  int s=1000*throw_id; std::string m(mode);
  printf("W/TKI %s ensemble member %d (seed base %d):\n", mode, throw_id, s);
  if(m=="fhc"){
    throw_group_ens_w({"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root"},Form("ens/throws/fakedata_fhc_run1_t%d.root",throw_id),2.192e20,2.3282e21,s+0);
    throw_group_ens_w({"xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root"},Form("ens/throws/fakedata_fhc_run2_t%d.root",throw_id),1.268e20,2.4934e21,s+1);
    throw_group_ens_w({"xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"},Form("ens/throws/fakedata_fhc_run4_t%d.root",throw_id),2.075e20,2.8335e21,s+2);
    throw_group_ens_w({"xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root"},Form("ens/throws/fakedata_fhc_run5_t%d.root",throw_id),2.231e20,1.9300e21,s+3);
  } else if(m=="rhc"){
    throw_group_ens_w({"xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/throws/fakedata_rhc_run1_t%d.root",throw_id),0.6053e20,8.9972e20,s+0);
    throw_group_ens_w({"xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/throws/fakedata_rhc_run2_t%d.root",throw_id),2.591e20,5.7865e21,s+1);
    throw_group_ens_w({"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root","xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root"},Form("ens/throws/fakedata_rhc_run3_t%d.root",throw_id),5.003e20,5.5185e21,s+2);
    throw_group_ens_w({"xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root","xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root"},Form("ens/throws/fakedata_rhc_run4_t%d.root",throw_id),2.883e20,3.2586e21,s+3);
  } else { printf("unknown mode %s\n", mode); }
}
