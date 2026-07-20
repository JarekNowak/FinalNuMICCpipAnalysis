// Regenerate the selection figures for the analysis note on a consistent naming
// scheme and style, from the current selection outputs:
//   histograms.root        : <Variable>_<CutStage>_<Sample> (Sample = Signal,
//                            Background, OOFV, EXT, Data, AllMC)
//   phase1_diagnostics.root : ph1_<var>_pi / ph1_<var>_pr (truth-particle PID)
// Produces (into report/figures/, overwriting):
//   nm1_<x>   : N-1 distribution (stacked MC truth + data), threshold line
//   final_<x> : final-cut distribution (stacked MC truth + data)
//   pid_<x>   : pion-candidate PID for true pi+- (red) vs true protons (blue)
// Run from the repo root:  root -l -b -q report/make_selection_figures.C
#include <vector>
#include <string>

static const char* FIG = "report/figures/";

// Stacked data/MC comparison: Signal+Background+OOFV+EXT stacked, Data overlaid.
void drawStack(TFile* f, const char* var, const char* cut, const char* xlab,
               const char* out, double thr=-1e30) {
  auto get=[&](const char* s){ return (TH1F*)f->Get(Form("%s_%s_%s",var,cut,s)); };
  TH1F* sig=get("Signal"); TH1F* bkg=get("Background");
  TH1F* oof=get("OOFV");   TH1F* ext=get("EXT");   TH1F* dat=get("Data");
  if(!sig||!bkg||!ext||!dat){ printf("  MISSING hists for %s (%s_%s)\n",out,var,cut); return; }
  sig->SetFillColor(kAzure+1); bkg->SetFillColor(kOrange+1);
  if(oof) oof->SetFillColor(kGreen+2); ext->SetFillColor(kGray+1);
  for(auto h:{sig,bkg,oof,ext}) if(h){h->SetLineColor(kBlack);h->SetLineWidth(1);}
  THStack* hs=new THStack("hs",Form(";%s;Events (3.283#times10^{20} POT)",xlab));
  hs->Add(sig); hs->Add(bkg); if(oof) hs->Add(oof); hs->Add(ext);
  dat->SetMarkerStyle(20); dat->SetMarkerSize(0.85); dat->SetLineColor(kBlack);
  TCanvas c("c","",720,540);
  double ymax=std::max(hs->GetMaximum(),dat->GetMaximum())*1.35;
  hs->SetMaximum(ymax); hs->SetMinimum(0); hs->Draw("hist"); dat->Draw("E1 same");
  TLegend* l=new TLegend(0.60,0.60,0.89,0.89);
  l->AddEntry(dat,"NuWro fake data","lep");
  l->AddEntry(sig,"Signal (CC#pi^{+})","f");
  l->AddEntry(bkg,"Beam background","f");
  if(oof) l->AddEntry(oof,"Out-of-FV","f");
  l->AddEntry(ext,"EXT (cosmic)","f");
  l->SetBorderSize(0); l->SetFillStyle(0); l->Draw();
  if(thr>-1e29){ TLine* ln=new TLine(thr,0,thr,ymax*0.78); ln->SetLineColor(kRed+1);
    ln->SetLineStyle(2); ln->SetLineWidth(2); ln->Draw();
    TLatex tx; tx.SetTextColor(kRed+1); tx.SetTextSize(0.035);
    tx.DrawLatex(thr,ymax*0.82,"cut"); }
  c.SaveAs(Form("%s%s.png",FIG,out));
}

// Truth-particle PID overlay (area-normalised): true pi+- vs true protons.
void drawPID(TFile* f, const char* base, const char* xlab, const char* out, double thr=-1e30) {
  TH1D* hpi=(TH1D*)f->Get(Form("ph1_%s_pi",base));
  TH1D* hpr=(TH1D*)f->Get(Form("ph1_%s_pr",base));
  if(!hpi||!hpr){ printf("  MISSING phase1 hists for %s\n",out); return; }
  hpi=(TH1D*)hpi->Clone(); hpr=(TH1D*)hpr->Clone();
  if(hpi->Integral()>0) hpi->Scale(1./hpi->Integral());
  if(hpr->Integral()>0) hpr->Scale(1./hpr->Integral());
  hpi->SetLineColor(kRed+1);  hpi->SetLineWidth(2); hpi->SetFillColorAlpha(kRed-9,0.4);
  hpr->SetLineColor(kAzure+2);hpr->SetLineWidth(2); hpr->SetFillColorAlpha(kAzure-9,0.4);
  hpi->SetTitle(Form(";%s;Area-normalised",xlab));
  double ymax=std::max(hpi->GetMaximum(),hpr->GetMaximum())*1.3;
  hpi->SetMaximum(ymax); hpi->SetMinimum(0); hpi->SetStats(0);
  TCanvas c("c","",720,540);
  hpi->Draw("hist"); hpr->Draw("hist same");
  TLegend* l=new TLegend(0.62,0.72,0.89,0.89);
  l->AddEntry(hpi,"true #pi^{#pm}","f"); l->AddEntry(hpr,"true proton","f");
  l->SetBorderSize(0); l->SetFillStyle(0); l->Draw();
  if(thr>-1e29){ TLine* ln=new TLine(thr,0,thr,ymax*0.85); ln->SetLineColor(kBlack);
    ln->SetLineStyle(2); ln->SetLineWidth(2); ln->Draw();
    TLatex tx; tx.SetTextSize(0.034); tx.DrawLatex(thr,ymax*0.88,"cut"); }
  c.SaveAs(Form("%s%s.png",FIG,out));
}

void make_selection_figures() {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TFile* fh = TFile::Open("histograms.root");

  // ---- N-1 distributions (variable at the stage just before its cut) ----
  drawStack(fh,"TopologicalScore","InFiducialVol","Topological score","nm1_topo",0.67);
  drawStack(fh,"MuonPID","Topological","Muon LLR PID score","nm1_mupid",0.2);
  drawStack(fh,"PionPID","MuonCandidate","Pion LLR PID score","nm1_pipid",0.1);
  drawStack(fh,"MuPiOpeningangle","ShowerCut","#mu#minus#pi opening angle [rad]","nm1_oa",2.6);
  drawStack(fh,"MuonTrackLength","Topological","Muon track length [cm]","nm1_mulen",10);
  drawStack(fh,"PionTrackLength","MuonCandidate","Pion track length [cm]","nm1_pilen",20);

  // ---- final-cut distributions (at Nonprotons) ----
  const char* C="Nonprotons";
  drawStack(fh,"MuonTrkMCSMom",C,"Muon momentum [GeV/c]","final_mumom");
  drawStack(fh,"MuonCosTheta", C,"cos#theta_{#mu}","final_mucosth");
  drawStack(fh,"MuonTrackLength",C,"Muon track length [cm]","final_mulen");
  drawStack(fh,"MuonPID",      C,"Muon LLR PID score","final_mupid");
  drawStack(fh,"PionTrkMuonMom",C,"Pion momentum [GeV/c]","final_pimom");
  drawStack(fh,"PionCosTheta", C,"cos#theta_{#pi}","final_picosth");
  drawStack(fh,"PionTrackLength",C,"Pion track length [cm]","final_pilen");
  drawStack(fh,"PionPID",      C,"Pion LLR PID score","final_pipid");
  drawStack(fh,"MuPiOpeningangle",C,"#mu#minus#pi opening angle [rad]","final_oa");
  drawStack(fh,"NPrimaryTracks",C,"N primary tracks","final_ntracks");
  fh->Close();

  // ---- pion-candidate PID (truth-particle separation) ----
  TFile* fp = TFile::Open("phase1_diagnostics.root");
  drawPID(fp,"bpion","Bragg-pion score","pid_bpion",0.08);
  drawPID(fp,"chipr","#chi^{2}_{proton}","pid_chipr");
  fp->Close();

  printf("Wrote selection figures to %s\n",FIG);
}
