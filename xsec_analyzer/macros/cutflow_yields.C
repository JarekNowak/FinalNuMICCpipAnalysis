// cutflow_yields.C — per-cut yield table + stacked figure (Fig 1 / Table 7 style),
// at FULL exposure, for FHC / RHC / combined. Reads the h_cutflow_tot and
// h_cutflow_sig counters written by the extended CC1mu1piXp selection into every
// processed file. Signal = sum of numuMC h_cutflow_sig (POT-scaled); Bkg = numuMC
// (tot - sig); EXT and Dirt from their own h_cutflow_tot with their scale factors.
// Data (blind) = Pred by construction, so the plot shows the MC-prediction stack.
//   RUN AFTER the reprocess that includes numuMC + EXT + dirt with the counters.
//   usage: root -l -b -q 'macros/cutflow_yields.C("fhc")'   // or "rhc","comb"
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void cutflow_yields(const char* mode="fhc") {
  const char* P="/data/uboone/processed/";
  const char* stage[10]={"None","InFiducialVol","Topological","MuonCandidate",
    "ContainedPion","MuonIn3Planes","PionIn3Planes","ShowerCut","OpeningAngle","Nonprotons"};
  // ---- per-config numuMC files + per-mode POT scale (Data/SummedMC) ----
  std::vector<Src> mc, ext, dirt;
  double sc_fhc=0.092402, sc_rhc=0.071666;             // pooled per-mode scales
  auto fhcMC=[&](){ for(auto r:{"Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
      "Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
      "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"})
      mc.push_back({std::string(P)+"xsec-ana-"+r+".root", sc_fhc}); };
  auto rhcMC=[&](){ for(auto r:{"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"})
      mc.push_back({std::string(P)+"xsec-ana-"+r+"_new_numi_flux_rhc_pandora_ntuple.root", sc_rhc});
      for(auto s:{"aa","ab","ac","ad","ae"})
      mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root", sc_rhc}); };
  if(std::string(mode)=="fhc") fhcMC();
  else if(std::string(mode)=="rhc") rhcMC();
  else { fhcMC(); rhcMC(); }
  // EXT + dirt: scale factors are placeholders to be set from the final trigger
  // accounting; left here so the driver is complete once the per-run EXT is wired.
  double sc_dirt = (std::string(mode)=="rhc") ? 0.0717*0.65 : 0.0924*0.65;
  dirt.push_back({std::string(P)+"xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root", sc_dirt});

  double sig[10]={0}, tot[10]={0}, ex[10]={0}, dt[10]={0};
  auto add=[&](std::vector<Src>&v, double* sACC, const char* hname){
    for(auto&s:v){ TFile* f=TFile::Open(s.file.c_str()); if(!f||f->IsZombie())continue;
      TH1D* h=(TH1D*)f->Get(hname); if(h) for(int c=0;c<10;c++) sACC[c]+=h->GetBinContent(c+1)*s.scale; f->Close(); } };
  add(mc, tot, "h_cutflow_tot"); add(mc, sig, "h_cutflow_sig");
  add(ext, ex, "h_cutflow_tot"); add(dirt, dt, "h_cutflow_tot");

  printf("\n=== %s cut-flow yields (full exposure) ===\n", mode);
  printf("%-15s %10s %10s %10s %10s %10s\n","Cut","Signal","Bkg","EXT","Dirt","Pred");
  for(int c=0;c<10;c++){ double bkg=tot[c]-sig[c], pred=sig[c]+bkg+ex[c]+dt[c];
    printf("%-15s %10.1f %10.1f %10.1f %10.1f %10.1f\n",stage[c],sig[c],bkg,ex[c],dt[c],pred); }

  // stacked figure
  gStyle->SetOptStat(0);
  TH1D* hs=new TH1D("hs_sig",Form("%s cut-flow;;events",mode),10,0,10);
  TH1D* hb=new TH1D("hb_bkg","",10,0,10), *he=new TH1D("he_ext","",10,0,10), *hd=new TH1D("hd_dirt","",10,0,10);
  for(int c=0;c<10;c++){ hs->SetBinContent(c+1,sig[c]); hb->SetBinContent(c+1,tot[c]-sig[c]); he->SetBinContent(c+1,ex[c]); hd->SetBinContent(c+1,dt[c]); hs->GetXaxis()->SetBinLabel(c+1,stage[c]); }
  // Okabe-Ito colorblind-safe fills: signal blue, bkg orange, EXT grey, dirt green
  hs->SetFillColor(TColor::GetColor("#0072B2")); hb->SetFillColor(TColor::GetColor("#E69F00"));
  he->SetFillColor(TColor::GetColor("#999999")); hd->SetFillColor(TColor::GetColor("#009E73"));
  THStack* st=new THStack("cf",Form("%s cut-flow;;events",mode)); st->Add(hd); st->Add(he); st->Add(hb); st->Add(hs);
  TCanvas c1("c1","",1000,600); c1.SetLogy(); c1.SetBottomMargin(0.22); st->Draw("hist"); st->GetXaxis()->LabelsOption("v");
  TLegend lg(0.7,0.72,0.88,0.88); lg.AddEntry(hs,"Signal","f"); lg.AddEntry(hb,"Bkg","f"); lg.AddEntry(he,"EXT","f"); lg.AddEntry(hd,"Dirt","f"); lg.Draw();
  TString out=Form("unfold_output/cutflow_yields_%s.pdf",mode); c1.SaveAs(out); printf("wrote %s\n",out.Data());
}
