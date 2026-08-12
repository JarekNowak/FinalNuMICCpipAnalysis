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
  // Okabe-Ito colorblind-safe palette + distinct line styles (redundant encoding,
  // so the four generators are separable in grayscale and for all colour-vision types)
  int gcol[4]   = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"),
                    TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00") }; // blue,green,purple,vermillion
  int gstyle[4] = { 1, 2, 7, 9 };
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
    // Build the generator overlays FIRST (per-bin dsigma*width -> divide by width) so
    // the y-axis can be scaled to include them (otherwise curves above the data clip
    // outside the frame). Warn if any generator is dropped for a bin-count mismatch.
    std::vector<TH1D*> gh; std::vector<int> gidx;
    for (int g = 0; g < 4; ++g) {
      TFile* fg = TFile::Open(Form("%s%s_%s_fte.root", GP, gens[g], gtag));
      if (!fg || fg->IsZombie()) { printf("  [warn] %s %s: no FTE file for %s\n",cfg,obs[o],gens[g]); continue; }
      TH1D* hfte = (TH1D*)fg->Get(Form("%s_fte", obs[o]));
      if (hfte && hfte->GetNbinsX()==hunf->GetNbinsX()) {
        TH1D* hg = (TH1D*)hunf->Clone(Form("g_%s_%d_%d",cfg,o,g)); hg->SetDirectory(0); hg->Reset();
        for (int b=1;b<=hg->GetNbinsX();++b) hg->SetBinContent(b, hfte->GetBinContent(b)/hg->GetBinWidth(b));
        hg->SetLineColor(gcol[g]); hg->SetLineStyle(gstyle[g]); hg->SetLineWidth(2); hg->SetMarkerSize(0);
        gh.push_back(hg); gidx.push_back(g); keep.push_back(hg);
      } else {
        printf("  [warn] %s %s: %s FTE bins=%d != data bins=%d -> DROPPED\n",
               cfg,obs[o],gens[g], hfte?hfte->GetNbinsX():-1, hunf->GetNbinsX());
      }
      fg->Close();
    }
    if (htru) { htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); keep.push_back(htru); }
    if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
    // y-axis max over data(+error), truth, tune and every generator curve
    double ymax = 0.;
    for (int b=1;b<=hunf->GetNbinsX();++b) ymax = std::max(ymax, hunf->GetBinContent(b)+hunf->GetBinError(b));
    if (htru) for (int b=1;b<=htru->GetNbinsX();++b) ymax = std::max(ymax, htru->GetBinContent(b));
    if (htun) for (int b=1;b<=htun->GetNbinsX();++b) ymax = std::max(ymax, htun->GetBinContent(b));
    for (auto hg : gh) for (int b=1;b<=hg->GetNbinsX();++b) ymax = std::max(ymax, hg->GetBinContent(b));
    hunf->SetMinimum(0); hunf->SetMaximum(1.30*ymax);
    hunf->Draw("E1");
    TLegend* lg = new TLegend(0.45,0.60,0.88,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
    lg->SetNColumns(2);
    lg->AddEntry(hunf, "Unfolded fake data", "lep");
    if (htru) { htru->SetLineColor(kGray+2); htru->SetLineWidth(2); htru->SetLineStyle(2); htru->Draw("hist same"); lg->AddEntry(htru,"Truth (A_{C})","l"); }
    if (htun) { htun->SetLineColor(kBlack); htun->SetLineWidth(2); htun->Draw("hist same"); lg->AddEntry(htun,"uB tune","l"); }
    for (size_t k=0;k<gh.size();++k) { gh[k]->Draw("hist same"); lg->AddEntry(gh[k], gens[gidx[k]], "l"); }
    hunf->Draw("E1 same");
    lg->Draw(); keep.push_back(lg); keep.push_back(hunf);
  }
  TString out = Form("unfold_output/dsigma_%s.pdf", cfg);
  c.SaveAs(out);
  printf("wrote %s\n", out.Data());
}
