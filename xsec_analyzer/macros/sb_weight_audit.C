// sb_weight_audit.C -- how many events, per region, have a CV weight product that the
// framework's safe_weight() replaces by unity (non-finite, negative, or >30), and the
// weighted sums under the two rules that were in use before 2026-09-03.
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void sb_weight_audit(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb/"; std::vector<Src> mc;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); }
  const int NR=5; const char* reg[NR]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic"}; const char* rn[NR]={"SR","CC0pi","multipi","pi0","cosmic"};
  const char* W="(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)";
  auto cnt=[&](const TString& expr,const TString& cut,double& acc){ for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str()); TH1D h("h","",1,-0.5,1.5); c.Draw(TString("0.5>>h"),expr+"*("+cut+")","goff"); acc+=h.Integral(0,2)*x.scale; } };
  printf("\n==== %s: CV-weight audit per region (POT-scaled)\n",mode);
  printf("%-8s %10s %10s %10s %10s | %12s %12s %12s\n","region","n_nonfin","n_neg","n_zero","n_gt30","sum safe","sum zero-rule","sum raw");
  for(int r=0;r<NR;r++){ double nf=0,ng=0,nz=0,nb=0,ss=0,sz=0,sr=0;
    cnt(TString("(!TMath::Finite")+W+")",reg[r],nf); cnt(TString("(TMath::Finite")+W+"&&"+W+"<0)",reg[r],ng); cnt(TString("(")+W+"==0)",reg[r],nz); cnt(TString("(TMath::Finite")+W+"&&"+W+">30)",reg[r],nb);
    cnt(TString("(TMath::Finite")+W+"&&"+W+">=0&&"+W+"<=30?"+W+":1)",reg[r],ss); cnt(TString("(TMath::Finite")+W+"&&"+W+">=0?"+W+":0)",reg[r],sz); cnt(TString("(TMath::Finite")+W+"?"+W+":0)",reg[r],sr);
    printf("%-8s %10.1f %10.1f %10.1f %10.1f | %12.1f %12.1f %12.1f\n",rn[r],nf,ng,nz,nb,ss,sz,sr); }
}
