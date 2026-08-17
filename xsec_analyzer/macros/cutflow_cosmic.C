// cutflow_cosmic.C — extend the final-cut yields with two new cosmic-rejection
// stages (crtveto==0, then + nu_flashmatch_score<15), for FHC / RHC / combined,
// to show the purity gain. Efficiencies from the reprocessed stv_tree (which now
// carries crtveto + nu_flashmatch_score); applied to the final-cut POT-weighted
// yields (signal/bkg/EXT/dirt) from the cut-flow table.
#include <vector>
#include <string>
void cutflow_cosmic(){
  const char* P="/data/uboone/processed/";
  auto fhc=[&](TChain&c){ for(auto r:{"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}) c.Add(Form("%sxsec-ana-%s.root",P,r)); };
  auto rhc=[&](TChain&c){ for(auto r:{"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}) c.Add(Form("%sxsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",P,r)); for(auto s:{"aa","ab","ac","ad","ae"}) c.Add(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",P,s)); };
  struct Cfg{const char* name; int mode; double yS,yB,yE,yD;};
  std::vector<Cfg> cfgs={{"FHC",0,993,653,121,1.3},{"RHC",1,1084,747,152,1.0},{"Combined",2,2077,1400,273,1.3}};
  TString SEL="CC1mu1piXp_Selected", SIG=SEL+" && CC1mu1piXp_MC_Signal", BKG=SEL+" && !CC1mu1piXp_MC_Signal";
  std::vector<std::pair<std::string,std::string>> stages={
    {"final selection",""},
    {"+ crtveto==0"," && crtveto==0"},
    {"+ nu_flashmatch<15"," && crtveto==0 && nu_flashmatch_score<15"}};
  for(auto&cf:cfgs){
    TChain mc("stv_tree"),ext("stv_tree");
    if(cf.mode==0)fhc(mc); else if(cf.mode==1)rhc(mc); else {fhc(mc);rhc(mc);}
    ext.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
    double S0=mc.GetEntries(SIG),B0=mc.GetEntries(BKG),E0=ext.GetEntries(SEL);
    printf("\n==== %s ====   (final-cut baseline: S=%.0f B=%.0f EXT=%.0f)\n",cf.name,cf.yS,cf.yB,cf.yE);
    printf("  %-20s %8s %8s %8s %6s %8s %8s\n","stage","Signal","Bkg","EXT","Dirt","Pred","Purity");
    for(auto&st:stages){
      double eS=mc.GetEntries(SIG+st.second.c_str())/S0, eB=mc.GetEntries(BKG+st.second.c_str())/B0, eE=ext.GetEntries(SEL+st.second.c_str())/E0;
      double s=cf.yS*eS,b=cf.yB*eB,e=cf.yE*eE,pred=s+b+e+cf.yD,pur=s/pred;
      printf("  %-20s %8.1f %8.1f %8.1f %6.1f %8.1f %7.1f%%\n",st.first.c_str(),s,b,e,cf.yD,pred,100*pur);
    }
  }
}
