// Framework-only differential cross-section figures (XSecAnalyzer / Version B).
// Reads the fresh Wiener-SVD unfolder outputs
//   xsec_analyzer/unfold_output/ccpi_Run1_<obs>_xsec.root
// and plots d(sigma)/dx [1e-38 cm^2/Ar per unit-x] with the FULL
// (systematic + statistical) uncertainty band, for the measured phase space
//   p_mu>0.15, p_pi>0.175 GeV/c, theta_mupi<2.6 rad.
//
// Convention (validated): the stored XsecUnits/<obs>/<obs>_total histogram is
// sigma-PER-BIN (= events/conversion_factor, conv const = 1374.331), NOT dsigma/dx.
//   sigma_bin[1e-38]   = 0.1 * content
//   dsigma/dx[1e-38]   = 0.1 * content / width
//   sigma_total[1e-38] = 0.1 * sum(content)
// The slice configs were not updated to the restricted phase space, so the
// physical measured bin edges are taken from the bin configs (below), not from
// the stored histogram axis (whose first p_mu/p_pi edge and last theta edge are
// stale).  Run from repo root:
//   root -l -b -q report/make_framework_figures.C

#include <vector>

static void one_obs(const char* obs, const char* sub, const char* xtitle,
                    const char* ytitle, std::vector<double> edges) {
  TString fn = Form("xsec_analyzer/unfold_output/ccpi_Run1_%s_xsec.root", obs);
  TFile* f = TFile::Open(fn);
  if (!f || f->IsZombie()) { printf("  %-8s : MISSING %s\n", obs, fn.Data()); return; }
  TH1D*     cv  = (TH1D*)    f->Get(Form("XsecUnits/%s/%s_total", sub, sub));
  TMatrixD* cov = (TMatrixD*)f->Get("XsecUnits/Covariances/total");
  if (!cv || !cov) { printf("  %-8s : missing CV/cov\n", obs); return; }
  int n = cv->GetNbinsX();
  if ((int)edges.size() != n+1) { printf("  %-8s : edge/bin mismatch (%d edges, %d bins)\n",
                                         obs, (int)edges.size(), n); return; }

  TH1D* h = new TH1D(Form("dsig_%s", obs), "", n, edges.data());
  double sig_tot = 0.0;
  for (int i = 1; i <= n; ++i) {
    double w   = h->GetBinWidth(i);
    double s   = 0.1 * cv->GetBinContent(i);            // sigma in bin [1e-38]
    double err = 0.1 * std::sqrt((*cov)(i-1, i-1));     // abs uncert on sigma_bin
    sig_tot += s;
    h->SetBinContent(i, s / w);
    h->SetBinError  (i, err / w);
  }

  gStyle->SetOptStat(0);
  TCanvas c("c", "", 760, 560);
  c.SetLeftMargin(0.15); c.SetBottomMargin(0.13);
  h->SetTitle(Form(";%s;%s", xtitle, ytitle));
  h->GetYaxis()->SetTitleOffset(1.45);
  h->SetMinimum(0.0);
  // uncertainty band (filled) drawn under the points
  TH1D* band = (TH1D*)h->Clone(Form("band_%s", obs));
  band->SetFillColorAlpha(kAzure-9, 0.8); band->SetFillStyle(1000);
  band->SetLineColor(kAzure-9); band->SetMarkerStyle(0);
  band->Draw("E2");
  h->SetLineColor(kBlue+2); h->SetLineWidth(2);
  h->SetMarkerStyle(20); h->SetMarkerColor(kBlue+2); h->SetMarkerSize(0.9);
  h->Draw("E1 X0 same");

  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.040);
  tl.DrawLatex(0.17, 0.92, "MicroBooNE NuMI FHC  CC1#mu1#pi^{#pm}Xp  (NuWro fake data)");
  tl.SetTextSize(0.038);
  tl.DrawLatex(0.55, 0.84, Form("#sigma_{tot} = %.3f #times10^{-38} cm^{2}/Ar", sig_tot));

  c.SaveAs(Form("report/figures/fw_xsec_%s.png", obs));
  printf("  %-8s : nbins=%2d  sigma_tot = %.4f x1e-38 cm^2/Ar -> fw_xsec_%s.png\n",
         obs, n, sig_tot, obs);
  f->Close();
}

void make_framework_figures() {
  gStyle->SetPalette(kBird);
  one_obs("pmu", "p_{#mu}",
          "true p_{#mu} [GeV/c]", "d#sigma/dp_{#mu} [10^{-38} cm^{2}/Ar/(GeV/c)]",
          {0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,
           1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.30,2.60,3.00});
  one_obs("ppi", "p_{#pi}",
          "true p_{#pi} [GeV/c]", "d#sigma/dp_{#pi} [10^{-38} cm^{2}/Ar/(GeV/c)]",
          {0.175,0.20,0.30,0.40,0.50,1.00});
  one_obs("costhmu", "cos#theta_{#mu}",
          "true cos#theta_{#mu}", "d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/Ar]",
          {-1.0,-0.5,0.0,0.2,0.4,0.55,0.65,0.75,0.82,0.88,0.93,0.97,1.0});
  one_obs("costhpi", "cos#theta_{#pi}",
          "true cos#theta_{#pi}", "d#sigma/dcos#theta_{#pi} [10^{-38} cm^{2}/Ar]",
          {-1.0,-0.7,-0.4,-0.2,0.0,0.2,0.35,0.5,0.65,0.78,0.88,0.95,1.0});
  one_obs("thmupi", "#theta_{#mu#pi}",
          "true #theta_{#mu#pi} [rad]", "d#sigma/d#theta_{#mu#pi} [10^{-38} cm^{2}/Ar/rad]",
          {0.0,0.314,0.628,0.942,1.257,1.571,1.885,2.199,2.513,2.600});
  printf("done -> report/figures/fw_xsec_*.png\n");
}
