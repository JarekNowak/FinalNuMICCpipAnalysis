// ext_study.C — plot the selection/cosmic-discriminating variables at the FINAL
// selection (CC1mu1piXp_Selected) for Signal / nu-background / EXT, to look for a
// variable that isolates the residual beam-off (cosmic) background. Uses TTree::Draw
// directly on the processed ntuples (no instrumentation needed), so any branch can
// be explored. Signal/bkg from numuMC (FHC+RHC), EXT from the beam-off sample.
//   usage: root -l -b -q macros/ext_study.C
#include <vector>
#include <string>
void ext_study() {
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

  gStyle->SetOptStat(0);
  const TString SEL="CC1mu1piXp_Selected";       // final selection flag
  const TString SIG=SEL+" && CC1mu1piXp_MC_Signal";
  const TString BKG=SEL+" && !CC1mu1piXp_MC_Signal";
  int cSig=TColor::GetColor("#0072B2"), cBkg=TColor::GetColor("#E69F00"), cExt=TColor::GetColor("#999999");

  struct V{const char* expr; const char* name; const char* xt; int nb; double lo,hi;};
  std::vector<V> vars = {
    {"CC1mu1piXp_candidate_muon_costh_reco","mucosth","muon cos#theta (reco)",20,-1,1},
    {"CC1mu1piXp_candidate_pion_costh_reco","picosth","pion cos#theta (reco)",20,-1,1},
    {"contained_fraction","contfrac","contained fraction",20,0,1.001},
    {"topological_score","topo","topological score",20,0,1},
  };
  for(auto&v:vars){
    TH1D hs("hs","",v.nb,v.lo,v.hi), hb("hb","",v.nb,v.lo,v.hi), he("he","",v.nb,v.lo,v.hi);
    mc.Draw(Form("%s>>hs",v.expr), SIG,"goff");
    mc.Draw(Form("%s>>hb",v.expr), BKG,"goff");
    ext.Draw(Form("%s>>he",v.expr), SEL,"goff");
    printf("[%s] raw entries: signal=%.0f  nu-bkg=%.0f  EXT=%.0f\n",v.name,hs.Integral(),hb.Integral(),he.Integral());
    // area-normalise for shape comparison
    for(TH1D* h:{&hs,&hb,&he}) if(h->Integral()>0) h->Scale(1.0/h->Integral());
    hs.SetLineColor(cSig); hb.SetLineColor(cBkg); he.SetLineColor(cExt);
    hs.SetLineWidth(3); hb.SetLineWidth(2); he.SetLineWidth(3); he.SetLineStyle(1);
    hs.SetFillColorAlpha(cSig,0.15); he.SetFillColorAlpha(cExt,0.25);
    double mx=std::max({hs.GetMaximum(),hb.GetMaximum(),he.GetMaximum()});
    hs.SetMaximum(mx*1.3); hs.SetTitle(Form(";%s;area-normalised",v.xt));
    TCanvas c("c","",800,600); hs.Draw("hist"); hb.Draw("hist same"); he.Draw("hist same");
    TLegend lg(0.62,0.74,0.88,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0);
    lg.AddEntry(&hs,"Signal","l"); lg.AddEntry(&hb,"#nu background","l"); lg.AddEntry(&he,"EXT (cosmic)","l"); lg.Draw();
    TString out=Form("unfold_output/ext_%s.pdf",v.name); c.SaveAs(out); printf("  wrote %s\n",out.Data());
  }
}
