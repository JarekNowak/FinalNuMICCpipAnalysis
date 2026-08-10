// ext_cosmic.C — re-measure the residual EXT using the flash/CRT cosmic taggers now
// carried through ProcessNTuples (crtveto, crthitpe, nu_flashmatch_score,
// best_cosmic_flashmatch_score, bdt_cosmic). For each cosmic cut, report signal /
// nu-bkg / EXT efficiency and the resulting FHC final-cut purity. crtveto is data-only
// (0 in MC), so crtveto==0 keeps all MC signal while removing CRT-tagged cosmics.
#include <string>
void ext_cosmic(){
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
  const double yS=993,yB=653,yE=121,yD=1.3; double p0=yS/(yS+yB+yE+yD);
  printf("baseline: signal=%.0f nu-bkg=%.0f EXT=%.0f | FHC purity=%.1f%%\n",S0,B0,E0,100*p0);
  // characterise EXT crtveto rate
  printf("EXT with crtveto==1 (CRT-tagged) at final selection: %lld / %.0f = %.1f%%\n\n",
    ext.GetEntries(SEL+" && crtveto==1"), E0, 100.*ext.GetEntries(SEL+" && crtveto==1")/E0);
  printf("%-34s %6s %6s %6s | %6s\n","cosmic cut (keep)","effS","effB","effE","purity");
  auto row=[&](const char* c){
    double eS=mc.GetEntries(SIG+" && "+TString(c))/S0, eB=mc.GetEntries(BKG+" && "+TString(c))/B0, eE=ext.GetEntries(SEL+" && "+TString(c))/E0;
    double s=yS*eS,b=yB*eB,e=yE*eE,pur=s/(s+b+e+yD);
    printf("%-34s %5.1f%% %5.1f%% %5.1f%% | %5.1f%%\n",c,100*eS,100*eB,100*eE,100*pur);
  };
  printf("-- CRT (data-only tag) --\n"); row("crtveto==0"); row("crtveto==0 && crthitpe<100");
  printf("-- neutrino flash-match score --\n"); row("nu_flashmatch_score<20"); row("nu_flashmatch_score<15"); row("nu_flashmatch_score<10");
  printf("-- best cosmic flash-match score --\n"); row("best_cosmic_flashmatch_score>5"); row("best_cosmic_flashmatch_score>8"); row("best_cosmic_flashmatch_score>10");
  printf("-- cosmic BDT --\n"); row("bdt_cosmic>0.2"); row("bdt_cosmic>0.4"); row("bdt_cosmic>0.5");
  printf("-- combined --\n"); row("crtveto==0 && bdt_cosmic>0.4");
}
