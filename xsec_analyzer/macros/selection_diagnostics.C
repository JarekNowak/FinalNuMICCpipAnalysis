// selection_diagnostics.C — N-1, final-cut candidate distributions, and background
// decomposition per configuration (FHC/RHC/combined), at full exposure, from the
// diagnostic histograms written by the extended CC1mu1piXp selection
// (h_nm1_topo/oa_{sig,bkg}, h_fin_{mupid,pipid,mulen,pilen}_{sig,bkg}, h_bkgcat).
// Signal + numuMC-background POT-scaled; EXT/dirt background from their own files.
//   RUN AFTER the reprocess that includes numuMC + EXT + dirt with the counters.
//   usage: root -l -b -q 'macros/selection_diagnostics.C("fhc")'  // or "rhc","comb"
#include <vector>
#include <string>
struct Src { std::string file; double scale; int isMC; };
void selection_diagnostics(const char* mode="fhc") {
  const char* P="/data/uboone/processed/";
  std::vector<Src> S;
  double s_fhc=0.092402, s_rhc=0.071666;
  auto fhc=[&](){ for(auto r:{"Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
    "Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
    "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}) S.push_back({std::string(P)+"xsec-ana-"+r+".root",s_fhc,1}); };
  auto rhc=[&](){ for(auto r:{"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"})
      S.push_back({std::string(P)+"xsec-ana-"+r+"_new_numi_flux_rhc_pandora_ntuple.root",s_rhc,1});
    for(auto s:{"aa","ab","ac","ad","ae"}) S.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",s_rhc,1}); };
  if(std::string(mode)=="fhc") fhc(); else if(std::string(mode)=="rhc") rhc(); else { fhc(); rhc(); }

  gStyle->SetOptStat(0);
  // stacked N-1 (signal + numuMC bkg) — EXT/dirt could be added when reprocessed
  auto stackNm1=[&](const char* base,const char* xt,double thr,const char* out){
    TH1D *hs=nullptr,*hb=nullptr;
    for(auto&s:S){ TFile*f=TFile::Open(s.file.c_str()); if(!f||f->IsZombie())continue;
      TH1D*a=(TH1D*)f->Get(Form("%s_sig",base)); TH1D*b=(TH1D*)f->Get(Form("%s_bkg",base));
      if(a){ if(!hs){hs=(TH1D*)a->Clone("hs");hs->SetDirectory(0);hs->Reset();} hs->Add(a,s.scale);}
      if(b){ if(!hb){hb=(TH1D*)b->Clone("hb");hb->SetDirectory(0);hb->Reset();} hb->Add(b,s.scale);} f->Close(); }
    if(!hs)return;
    hs->SetFillColor(TColor::GetColor("#0072B2")); hb->SetFillColor(TColor::GetColor("#E69F00")); // Okabe-Ito
    THStack st("st",Form("%s (%s);%s;events",base,mode,xt)); st.Add(hb); st.Add(hs);
    TCanvas c("c","",800,600); st.Draw("hist");
    TLine ln(thr,0,thr,hs->GetMaximum()*1.1); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw();
    TLegend lg(0.7,0.75,0.88,0.88); lg.AddEntry(hs,"Signal","f"); lg.AddEntry(hb,"Bkg","f"); lg.Draw();
    TString o=Form("unfold_output/%s_%s.pdf",out,mode); c.SaveAs(o); printf("wrote %s\n",o.Data());
  };
  stackNm1("h_nm1_topo","topological score",0.67,"nm1_topo");
  stackNm1("h_nm1_oa","#theta_{#mu#pi} [rad]",2.6,"nm1_oa");
  stackNm1("h_fin_mupid","muon LLR PID",0.2,"final_mupid");
  stackNm1("h_fin_pipid","pion LLR PID",0.1,"final_pipid");
  stackNm1("h_fin_mulen","muon length [cm]",10,"final_mulen");
  stackNm1("h_fin_pilen","pion length [cm]",20,"final_pilen");

  // background decomposition (POT-weighted) by category at the final cut
  TH1D* hcat=nullptr;
  for(auto&s:S){ TFile*f=TFile::Open(s.file.c_str()); if(!f||f->IsZombie())continue;
    TH1D*h=(TH1D*)f->Get("h_bkgcat"); if(h){ if(!hcat){hcat=(TH1D*)h->Clone("hcat");hcat->SetDirectory(0);hcat->Reset();} hcat->Add(h,s.scale);} f->Close(); }
  if(hcat){ printf("\n=== %s background decomposition (category : POT-weighted events) ===\n",mode);
    for(int b=1;b<=hcat->GetNbinsX();++b) if(hcat->GetBinContent(b)>0.05) printf("  cat %2d : %8.1f\n",b-1,hcat->GetBinContent(b)); }
}
