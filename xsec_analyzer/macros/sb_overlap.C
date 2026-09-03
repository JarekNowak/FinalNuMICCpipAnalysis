// sb_overlap.C -- entry-level overlap between the signal region (Selected) and each control
// region, and between control regions, on the sb_-instrumented files (MC, CV-weighted and
// POT-scaled; plus the per-run fake data, unweighted). Each TTree entry is one event, so an
// entry passing two flags is one event in both regions.
//   root -l -b -q 'macros/sb_overlap.C("fhc")'   (rhc)
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void sb_overlap(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb/"; std::vector<Src> mc; std::vector<std::string> data;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); for(auto r:{"run1","run2","run4","run5"}) data.push_back(std::string(P)+"xsec-ana-fakedata_fhc_"+r+".root"); }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066});
    for(auto r:{"run1","run2","run3","run4"}) data.push_back(std::string(P)+"xsec-ana-fakedata_rhc_"+r+".root"); }
  // Framework safe_weight rule: non-finite, negative or >30 -> 1 (one authoritative rule, 2026-09-03)
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* reg[5]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic"};
  const char* nm[5]={"SR","CC0pi","multipi","pi0","cosmic"};
  auto wsum=[&](TChain& c,const TString& cut,bool w)->double{ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",(w?TString(CVW):TString("1"))+"*("+cut+")","goff"); return h.Integral(0,2); };
  double M[5][5]={{0}}, Dm[5][5]={{0}};
  for(int i=0;i<5;i++) for(int j=i;j<5;j++){ TString cut=TString(reg[i])+" && "+reg[j];
    for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str()); M[i][j]+=wsum(c,cut,true)*x.scale; }
    for(auto&d:data){ TChain c("stv_tree"); c.Add(d.c_str()); Dm[i][j]+=wsum(c,cut,false); } }
  printf("\n==== %s: region overlaps (MC POT-scaled CV-weighted | fake data raw). Diagonal = region size.\n",mode);
  printf("%-9s","");for(int j=0;j<5;j++)printf("%14s",nm[j]);printf("\n");
  for(int i=0;i<5;i++){ printf("%-9s",nm[i]); for(int j=0;j<5;j++){ if(j<i)printf("%14s",""); else printf("%8.1f|%5.0f",M[i][j],Dm[i][j]); } printf("\n"); }
}
