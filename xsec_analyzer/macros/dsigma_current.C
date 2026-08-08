// dsigma_current.C — differential cross-section result figures for the CURRENT
// Poisson-thrown central-value fake data. Reads the UnfolderNuMI closure dumps
// closure_hists_xsec_<CFG>_<obs>.root (h_unfolded_nuwro = unfolded fake data,
// h_fakedata_truth = A_C-smeared truth, h_genie_tune = MicroBooNE tune), overlays
// the four generator FTE predictions, and writes one multi-panel figure per config.
//   usage: root -l -b -q 'macros/dsigma_current.C("FHC5","newg4")'  (or "RHCFULL","rhc" / "COMB","comb")
void dsigma_current(const char* cfg = "FHC5", const char* gtag = "newg4") {
  const char* PROC = "/data/uboone/processed/";
  const char* GP   = "../generator_predictions/newg4/";
  const char* obs[5]    = {"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* obsX[5]   = {"p_{#mu} [GeV/c]","p_{#pi} [GeV/c]","cos#theta_{#mu}",
                           "cos#theta_{#pi}","#theta_{#mu#pi} [rad]"};
  const char* gens[4]   = {"genie","gibuu","neut","nuwro"};
  int gcol[4]           = {kRed+1, kGreen+2, kMagenta+1, kBlue+1};
  gStyle->SetOptStat(0);
  TCanvas c(Form("ds_%s",cfg), "", 1600, 900);
  c.Divide(3, 2);
  std::vector<TObject*> keep;
  for (int o = 0; o < 5; ++o) {
    c.cd(o+1);
    TFile* f = TFile::Open(Form("%sclosure_hists_xsec_%s_%s.root", PROC, cfg, obs[o]));
    if (!f || f->IsZombie()) { printf("  missing closure %s %s\n", cfg, obs[o]); continue; }
    TH1D* hunf = (TH1D*)f->Get("h_unfolded_nuwro");   // unfolded fake data
    TH1D* htru = (TH1D*)f->Get("h_fakedata_truth");   // A_C-smeared truth
    TH1D* htun = (TH1D*)f->Get("h_genie_tune");       // uB tune
    if (!hunf) continue;
    hunf = (TH1D*)hunf->Clone(); hunf->SetDirectory(0);
    hunf->SetTitle(Form("%s;%s;d#sigma/dx [10^{-38} cm^{2}/Ar]", obs[o], obsX[o]));
    hunf->SetMarkerStyle(20); hunf->SetMarkerSize(0.8); hunf->SetLineColor(kBlack); hunf->SetMarkerColor(kBlack);
    hunf->SetMinimum(0);
    hunf->Draw("E1");
    TLegend* lg = new TLegend(0.45,0.62,0.88,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    lg->AddEntry(hunf, "Unfolded fake data", "lep");
    if (htru) { htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); htru->SetLineColor(kGray+2); htru->SetLineWidth(2); htru->SetLineStyle(2); htru->Draw("hist same"); lg->AddEntry(htru,"Fake-data truth (A_{C})","l"); keep.push_back(htru); }
    if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); htun->SetLineColor(kBlack); htun->SetLineWidth(2); htun->Draw("hist same"); lg->AddEntry(htun,"MicroBooNE tune","l"); keep.push_back(htun); }
    // overlay the four generator FTE predictions (per-bin dsigma*width -> divide by width)
    for (int g = 0; g < 4; ++g) {
      TFile* fg = TFile::Open(Form("%s%s_%s_fte.root", GP, gens[g], gtag));
      if (!fg || fg->IsZombie()) continue;
      TH1D* hfte = (TH1D*)fg->Get(Form("%s_fte", obs[o]));
      if (hfte && hfte->GetNbinsX()==hunf->GetNbinsX()) {
        TH1D* hg = (TH1D*)hunf->Clone(Form("g_%s_%d_%d",cfg,o,g)); hg->SetDirectory(0); hg->Reset();
        for (int b=1;b<=hg->GetNbinsX();++b) hg->SetBinContent(b, hfte->GetBinContent(b)/hg->GetBinWidth(b));
        hg->SetLineColor(gcol[g]); hg->SetLineWidth(2); hg->SetMarkerSize(0); hg->Draw("hist same");
        lg->AddEntry(hg, gens[g], "l"); keep.push_back(hg);
      }
      fg->Close();
    }
    hunf->Draw("E1 same");
    lg->Draw(); keep.push_back(lg); keep.push_back(hunf);
  }
  TString out = Form("unfold_output/dsigma_%s.pdf", cfg);
  c.SaveAs(out);
  printf("wrote %s\n", out.Data());
}
