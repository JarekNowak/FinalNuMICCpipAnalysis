// detvar_plots.C — overlay all nine detector-variation samples (CV + 8
// variations) for one observable, showing both the reco spectrum (top pad) and
// the true signal distribution (bottom pad). Each sample is POT-scaled to the
// CV exposure so the curves are directly comparable. The true pad drops the
// final catch-all (non-signal) true bin.
//
// Source: the unweighted_0_reco / unweighted_0_true histograms stored per input
// file in the univmake output (ccpi_Run1_<obs>_univmake.root).
//
// Run: root -l -b -q 'detvar_plots.C("costhmu")'
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

void detvar_plots(const char* obs) {
  struct DV { const char* tag; const char* label; double pot; int color; int style; };
  const double CVPOT = 7.630913e20;
  std::vector<DV> dv = {
    {"CV_slim",            "CV",             7.630913e20, kBlack,     1},
    {"LYDown_slim",        "LY Down",        7.053740e20, kRed+1,     2},
    {"LYRayleigh_slim",    "LY Rayleigh",    7.623732e20, kBlue,      2},
    {"Recombination_slim", "Recombination",  7.537606e20, kGreen+2,   2},
    {"SCE_slim",           "SCE",            7.532321e20, kMagenta+1, 2},
    {"WireModThetaXZ_slim","WireMod #thetaXZ",7.628598e20,kOrange+7, 7},
    {"WireModThetaYZ_slim","WireMod #thetaYZ",7.653130e20,kCyan+2,   7},
    {"WireModX_slim",      "WireMod X",      7.600698e20, kViolet,    7},
    {"WireModYZ_slim",     "WireMod YZ",     7.509717e20, kAzure+1,   7},
  };

  TFile* f = TFile::Open(Form("/data/uboone/processed/ccpi_Run1_%s_univmake.root", obs));
  TDirectoryFile* top = (TDirectoryFile*)f->Get(Form("ccpi_CC1mu1piXp_%s_1D", obs));
  if (!top) { printf("no top dir for %s\n", obs); return; }

  std::vector<TH1D*> reco, tru;
  for (auto& d : dv) {
    TDirectoryFile* sub = (TDirectoryFile*)top->Get(
      Form("+data+uboone+processed+xsec-ana-numi_detvars_pion_ntuples_run1_fhc_%s.root", d.tag));
    if (!sub) { printf("missing subdir %s\n", d.tag); return; }
    TH1D* hr = (TH1D*)sub->Get("unweighted_0_reco");
    TH1D* ht = (TH1D*)sub->Get("unweighted_0_true");
    hr = (TH1D*)hr->Clone(Form("reco_%s", d.tag));
    ht = (TH1D*)ht->Clone(Form("true_%s", d.tag));
    hr->SetDirectory(0); ht->SetDirectory(0);
    double s = CVPOT / d.pot;                 // POT-scale to CV
    hr->Scale(s); ht->Scale(s);
    // drop the final catch-all (non-signal) true bin so it doesn't swamp the pad
    int nt = ht->GetNbinsX(); ht->SetBinContent(nt,0); ht->SetBinError(nt,0);
    for (TH1D* h : {hr, ht}) {
      h->SetStats(0); h->SetLineColor(d.color); h->SetLineStyle(d.style);
      h->SetLineWidth(d.color==kBlack ? 3 : 2); h->SetTitle("");
    }
    reco.push_back(hr); tru.push_back(ht);
  }

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas(Form("dv_%s", obs), obs, 1500, 1100);
  TPad* p1 = new TPad("p1","",0.0,0.52,0.78,1.0); p1->SetBottomMargin(0.11); p1->SetLeftMargin(0.13); p1->Draw();
  TPad* p2 = new TPad("p2","",0.0,0.02,0.78,0.52); p2->SetBottomMargin(0.13); p2->SetLeftMargin(0.13); p2->Draw();
  TPad* pl = new TPad("pl","",0.78,0.02,1.0,1.0); pl->Draw();

  // reco pad
  p1->cd();
  double rmax=0; for (auto h: reco) rmax = std::max(rmax, h->GetMaximum());
  reco[0]->GetYaxis()->SetRangeUser(0., rmax*1.25);
  reco[0]->GetYaxis()->SetTitle("selected events (reco), POT-scaled");
  reco[0]->GetXaxis()->SetTitle(Form("%s reco bin", obs));
  reco[0]->Draw("hist");
  for (size_t i=1;i<reco.size();++i) reco[i]->Draw("hist same");
  reco[0]->Draw("hist same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.05);
  t.DrawLatex(0.16,0.93, Form("Detector variations: %s (reco)", obs));

  // true pad (drop the final catch-all true bin)
  p2->cd();
  int ntb = tru[0]->GetNbinsX();
  double tmax=0; for (auto h: tru){ for(int b=1;b<ntb;++b) tmax=std::max(tmax,h->GetBinContent(b)); }
  tru[0]->GetYaxis()->SetRangeUser(0., tmax*1.25);
  tru[0]->GetXaxis()->SetRange(1, ntb-1);
  tru[0]->GetYaxis()->SetTitle("true signal events, POT-scaled");
  tru[0]->GetXaxis()->SetTitle(Form("%s true bin (signal)", obs));
  tru[0]->Draw("hist");
  for (size_t i=1;i<tru.size();++i) tru[i]->Draw("hist same");
  tru[0]->Draw("hist same");
  t.DrawLatex(0.16,0.93, Form("Detector variations: %s (true signal)", obs));

  // legend
  pl->cd();
  TLegend* lg = new TLegend(0.02,0.3,0.98,0.72);
  lg->SetBorderSize(0); lg->SetTextSize(0.05);
  for (size_t i=0;i<dv.size();++i) lg->AddEntry(reco[i], dv[i].label, "l");
  lg->Draw();

  c->SaveAs(Form("unfold_output/detvar_%s.pdf", obs));
  printf("wrote unfold_output/detvar_%s.pdf\n", obs);
}
