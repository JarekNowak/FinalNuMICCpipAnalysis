// dsigma_ccpi1p.C — differential cross-section montage for the proton-tagged
// CC1mu1pi1p (W/TKI) observables. Reads the UnfolderNuMI closure dumps
// closure_hists_xsec_ccpi1p_<CFG>_<obs>.root (h_unfolded_nuwro = unfolded fake data,
// h_fakedata_truth = A_C-smeared truth, h_genie_tune = MicroBooNE tune) and draws one
// 6-panel figure per config. No external generator predictions exist for these
// observables yet, so only data/truth/tune are shown. Auto-scales the y-axis to
// include every curve (no clipping).
//   usage: root -l -b -q 'macros/dsigma_ccpi1p.C("FHC5")'  (or RHCFULL / COMB)
#include <vector>
void dsigma_ccpi1p(const char* cfg = "FHC5") {
  const char* PROC = "/data/uboone/processed/";
  const char* obs[6]  = {"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  const char* obsX[6] = {"W_{#pi p} [GeV/c^{2}]","W_{had} [GeV/c^{2}]",
                         "#deltap_{T} [GeV/c]","#delta#alpha_{T} [deg]",
                         "#delta#phi_{T} [deg]","p_{n} [GeV/c]"};
  gStyle->SetOptStat(0);
  TCanvas c(Form("ds1p_%s",cfg), "", 1500, 900);
  c.Divide(3, 2);
  std::vector<TObject*> keep;
  for (int o = 0; o < 6; ++o) {
    c.cd(o+1);
    TFile* f = TFile::Open(Form("%sclosure_hists_xsec_ccpi1p_%s_%s.root", PROC, cfg, obs[o]));
    if (!f || f->IsZombie()) { printf("  missing closure %s %s\n", cfg, obs[o]); continue; }
    TH1D* hunf = (TH1D*)f->Get("h_unfolded_nuwro");
    TH1D* htru = (TH1D*)f->Get("h_fakedata_truth");
    TH1D* htun = (TH1D*)f->Get("h_genie_tune");
    if (!hunf) { f->Close(); continue; }
    hunf = (TH1D*)hunf->Clone(); hunf->SetDirectory(0); keep.push_back(hunf);
    if (htru) { htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); keep.push_back(htru); }
    if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
    f->Close();
    hunf->SetTitle(Form("%s;%s;d#sigma/dx [10^{-38} cm^{2}/Ar]", obs[o], obsX[o]));
    hunf->SetMarkerStyle(20); hunf->SetMarkerSize(0.8);
    hunf->SetLineColor(kBlack); hunf->SetMarkerColor(kBlack);
    // y-axis max over data(+error), truth and tune
    double ymax = 0.;
    for (int b=1;b<=hunf->GetNbinsX();++b) ymax = std::max(ymax, hunf->GetBinContent(b)+hunf->GetBinError(b));
    if (htru) for (int b=1;b<=htru->GetNbinsX();++b) ymax = std::max(ymax, htru->GetBinContent(b));
    if (htun) for (int b=1;b<=htun->GetNbinsX();++b) ymax = std::max(ymax, htun->GetBinContent(b));
    hunf->SetMinimum(0); hunf->SetMaximum(1.30*ymax);
    hunf->Draw("E1");
    TLegend* lg = new TLegend(0.40,0.72,0.88,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
    lg->AddEntry(hunf, "Unfolded fake data", "lep");
    if (htru) { htru->SetLineColor(kGray+2); htru->SetLineWidth(2); htru->SetLineStyle(2); htru->Draw("hist same"); lg->AddEntry(htru,"Fake-data truth (A_{C})","l"); }
    if (htun) { htun->SetLineColor(TColor::GetColor("#0072B2")); htun->SetLineWidth(2); htun->Draw("hist same"); lg->AddEntry(htun,"MicroBooNE tune","l"); }
    hunf->Draw("E1 same");
    lg->Draw(); keep.push_back(lg);
  }
  TString out = Form("unfold_output/dsigma_ccpi1p_%s.pdf", cfg);
  c.SaveAs(out);
  printf("wrote %s\n", out.Data());
}
