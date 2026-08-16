// cutflow_yields_1p.C — per-cut yield table + stacked figure for the proton-tagged
// CC1mu1pi1p subsample, full exposure, FHC/RHC/combined. Reads the 10-stage
// h_cutflow_tot/h_cutflow_sig counters from the w/-reprocessed files (filled by the
// inherited CC1mu1piXp cut sequence, with the CC1mu1pi1p signal definition) and appends
// an 11th stage, "IdentProton", from the CC1mu1pi1p_Selected flag in the tree -- the
// leading-proton identification that defines the subsample. Signal = numuMC sig
// (POT-scaled); Bkg = numuMC (tot-sig); EXT/Dirt from their own counters + scales.
//   usage: root -l -b -q 'macros/cutflow_yields_1p.C("fhc")'   // or "rhc","comb"
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void cutflow_yields_1p(const char* mode="fhc") {
  const char* P="/data/uboone/processed/w/";
  const char* stage[11]={"None","InFiducialVol","Topological","MuonCandidate",
    "ContainedPion","MuonIn3Planes","PionIn3Planes","ShowerCut","OpeningAngle",
    "Nonprotons","IdentProton"};
  std::vector<Src> mc, ext, dirt;
  auto fhcMC=[&](){
    const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
      "Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
      "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
    double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root", sc[i]}); };
  auto rhcMC=[&](){
    const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
    double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root", sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"})
      mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root", 0.09066}); };
  if(std::string(mode)=="fhc") fhcMC();
  else if(std::string(mode)=="rhc") rhcMC();
  else { fhcMC(); rhcMC(); }
  // NUMI_EXT_OCC: 2% NuMI beam-occupancy factor the framework applies on top of the
  // trigger ratio (SystematicsCalculator.cxx); it was missing here, over-counting EXT 2%.
  const double NUMI_EXT_OCC = 0.98;
  double sc_ext = NUMI_EXT_OCC * ( (std::string(mode)=="fhc") ? 5.9313 : (std::string(mode)=="rhc") ? 6.1584 : 12.0898 );
  ext.push_back({std::string(P)+"xsec-ana-beamoff_run1Andrun3.root", sc_ext});
  // "comb" sums both modes' exposures, as sc_ext does; the old two-branch ternary sent
  // comb down the FHC branch and scaled combined dirt to FHC exposure alone.
  double sc_dirt = ( (std::string(mode)=="fhc") ? 0.092402
                   : (std::string(mode)=="rhc") ? 0.071666
                   : (0.092402 + 0.071666) ) * 0.65;
  dirt.push_back({std::string(P)+"xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root", sc_dirt});

  double sig[11]={0}, tot[11]={0}, ex[11]={0}, dt[11]={0};
  // stages 1-10 from the histogram counters
  auto addH=[&](std::vector<Src>&v, double* sACC, const char* hname){
    for(auto&s:v){ TFile* f=TFile::Open(s.file.c_str()); if(!f||f->IsZombie())continue;
      TH1D* h=(TH1D*)f->Get(hname); if(h) for(int c=0;c<10;c++) sACC[c]+=h->GetBinContent(c+1)*s.scale; f->Close(); } };
  addH(mc, tot, "h_cutflow_tot"); addH(mc, sig, "h_cutflow_sig");
  addH(ext, ex, "h_cutflow_tot"); addH(dirt, dt, "h_cutflow_tot");
  // stage 11 (identified proton) from the tree flags
  auto addP=[&](std::vector<Src>&v, double* totA, double* sigA){
    for(auto&s:v){ TChain c("stv_tree"); c.Add(s.file.c_str());
      totA[10]+=c.GetEntries("CC1mu1pi1p_Selected")*s.scale;
      if(sigA) sigA[10]+=c.GetEntries("CC1mu1pi1p_Selected && CC1mu1pi1p_MC_Signal")*s.scale; } };
  addP(mc, tot, sig); addP(ext, ex, nullptr); addP(dirt, dt, nullptr);

  printf("\n=== %s CC1mu1pi1p cut-flow yields (full exposure) ===\n", mode);
  printf("%-15s %10s %10s %10s %10s %10s\n","Cut","Signal","Bkg","EXT","Dirt","Pred");
  for(int c=0;c<11;c++){ double bkg=tot[c]-sig[c], pred=sig[c]+bkg+ex[c]+dt[c];
    printf("%-15s %10.1f %10.1f %10.1f %10.1f %10.1f\n",stage[c],sig[c],bkg,ex[c],dt[c],pred); }

  gStyle->SetOptStat(0);
  TH1D* hs=new TH1D("hs_sig",Form("%s CC1mu1pi1p cut-flow;;events",mode),11,0,11);
  TH1D* hb=new TH1D("hb_bkg","",11,0,11), *he=new TH1D("he_ext","",11,0,11), *hd=new TH1D("hd_dirt","",11,0,11);
  for(int c=0;c<11;c++){ hs->SetBinContent(c+1,sig[c]); hb->SetBinContent(c+1,tot[c]-sig[c]); he->SetBinContent(c+1,ex[c]); hd->SetBinContent(c+1,dt[c]); hs->GetXaxis()->SetBinLabel(c+1,stage[c]); }
  hs->SetFillColor(TColor::GetColor("#0072B2")); hb->SetFillColor(TColor::GetColor("#E69F00"));
  he->SetFillColor(TColor::GetColor("#999999")); hd->SetFillColor(TColor::GetColor("#009E73"));
  THStack* st=new THStack("cf",Form("%s CC1mu1pi1p cut-flow;;events",mode)); st->Add(hd); st->Add(he); st->Add(hb); st->Add(hs);
  TCanvas c1("c1","",1000,600); c1.SetLogy(); c1.SetBottomMargin(0.22); st->Draw("hist"); st->GetXaxis()->LabelsOption("v");
  TLegend lg(0.7,0.72,0.88,0.88); lg.AddEntry(hs,"Signal","f"); lg.AddEntry(hb,"Bkg","f"); lg.AddEntry(he,"EXT","f"); lg.AddEntry(hd,"Dirt","f"); lg.Draw();
  TString out=Form("unfold_output/cutflow_yields_1p_%s.pdf",mode); c1.SaveAs(out); printf("wrote %s\n",out.Data());
}
