// ext_scan.C — quantitative cut scan for residual EXT removal. For each candidate
// cut, compute signal / nu-bkg / EXT efficiencies and the resulting FHC final-cut
// purity and S*purity figure of merit (baseline yields S=993, B=653, EXT=121, dirt=1.3).
#include <vector>
#include <string>
void ext_scan() {
  const char* P="/data/uboone/processed/";
  TChain mc("stv_tree"), ext("stv_tree");
  const char* mcf[]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
    "Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc",
    "Run1_rhc_new_numi_flux_rhc_pandora_ntuple","Run2_rhc_new_numi_flux_rhc_pandora_ntuple",
    "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple","Run4b_rhc_new_numi_flux_rhc_pandora_ntuple",
    "Run4c_rhc_new_numi_flux_rhc_pandora_ntuple"};
  for(auto f:mcf) mc.Add(Form("%sxsec-ana-%s.root",P,f));
  for(auto s:{"aa","ab","ac","ad","ae"}) mc.Add(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",P,s));
  ext.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
  TString SEL="CC1mu1piXp_Selected", SIG=SEL+" && CC1mu1piXp_MC_Signal", BKG=SEL+" && !CC1mu1piXp_MC_Signal";
  double S0=mc.GetEntries(SIG), B0=mc.GetEntries(BKG), E0=ext.GetEntries(SEL);
  // FHC final-cut yields (POT-weighted)
  const double yS=993.0, yB=653.0, yE=121.0, yD=1.3;
  double p0=yS/(yS+yB+yE+yD);
  printf("baseline: signal=%.0f nu-bkg=%.0f EXT=%.0f | FHC purity=%.1f%%  S*p=%.0f\n\n",S0,B0,E0,100*p0,yS*p0);
  printf("%-46s %6s %6s %6s | %6s %6s\n","cut","effS","effB","effE","purity","S*p");
  auto scan=[&](std::vector<std::string> cuts){
    for(auto&c:cuts){ TString cc=c.c_str();
      double eS=mc.GetEntries(SIG+" && "+cc)/S0, eB=mc.GetEntries(BKG+" && "+cc)/B0, eE=ext.GetEntries(SEL+" && "+cc)/E0;
      double s=yS*eS,b=yB*eB,e=yE*eE, pur=s/(s+b+e+yD);
      printf("%-46s %5.1f%% %5.1f%% %5.1f%% | %5.1f%% %6.0f\n",c.c_str(),100*eS,100*eB,100*eE,100*pur,s*pur); }
    printf("\n"); };
  printf("--- topological score ---\n");
  scan({"topological_score>0.70","topological_score>0.75","topological_score>0.80","topological_score>0.85","topological_score>0.90","topological_score>0.93","topological_score>0.95"});
  printf("--- reco vertex Y (cosmic entry from top/bottom) ---\n");
  scan({"reco_nu_vtx_sce_y>-80","abs(reco_nu_vtx_sce_y)<90","abs(reco_nu_vtx_sce_y)<100","reco_nu_vtx_sce_y<100 && reco_nu_vtx_sce_y>-100"});
  printf("--- muon momentum (reco) ---\n");
  scan({"CC1mu1piXp_candidate_muon_mom_reco<1.5","CC1mu1piXp_candidate_muon_mom_reco<1.0","CC1mu1piXp_candidate_muon_mom_reco<0.8"});
  printf("--- combined: topo>0.80 + muon costh>-0.25 ---\n");
  scan({"topological_score>0.80 && CC1mu1piXp_candidate_muon_costh_reco>-0.25","topological_score>0.85 && CC1mu1piXp_candidate_muon_mom_reco<1.5"});
}
