// Reconstruction-resolution plots, one per measured observable:
//   momenta : relative (reco-true)/true  vs true
//   angles  : absolute (reco-true)        vs true   (relative diverges as cos->0 / theta->0)
// Source: framework GENIE-overlay processed tree, selected signal (new phase space
// pmu>0.15, ppi>0.175, theta_mupi<2.6 is encoded in CC1mu1piXp_MC_Signal).
// Writes PNGs into report/figures/.  Run from repo root:
//   root -l -b -q report/make_resolution_figures.C
void one_panel(TTree* t, const char* cut, const char* var_true, const char* metric_expr,
               const char* title, const char* xlab, const char* ylab,
               int nx, double xlo, double xhi, double ylo, double yhi, const char* out) {
  TString draw = Form("%s : %s >> h2(%d,%g,%g,80,%g,%g)", metric_expr, var_true, nx, xlo, xhi, ylo, yhi);
  t->Draw(draw, cut, "goff");
  TH2D* h2 = (TH2D*)gDirectory->Get("h2");
  if (!h2 || h2->GetEntries()==0) { printf("  %-14s : no entries\n", out); return; }
  // overall resolution = mean/RMS of the metric over the whole sample
  TH1D* hy = h2->ProjectionY("hy");
  double mean = hy->GetMean(), rms = hy->GetRMS();
  h2->SetTitle(Form("%s  (mean=%.3f, RMS=%.3f);%s;%s", title, mean, rms, xlab, ylab));
  TProfile* prof = h2->ProfileX("prof", 1, -1, "");          // mean metric per true bin
  TCanvas c("c","",760,560); c.SetRightMargin(0.14);
  gStyle->SetOptStat(0);
  h2->Draw("colz");
  TLine* z = new TLine(xlo,0,xhi,0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw();
  prof->SetLineColor(kRed+1); prof->SetMarkerColor(kRed+1); prof->SetMarkerStyle(20);
  prof->SetMarkerSize(0.8); prof->SetLineWidth(2); prof->Draw("E1 same");
  c.SaveAs(Form("report/figures/%s.png", out));
  printf("  %-14s mean=%+.3f RMS=%.3f  (%s)\n", out, mean, rms, title);
}

void make_resolution_figures() {
  gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(60);
  TFile* f = TFile::Open("/data/uboone/processed/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t = (TTree*)f->Get("stv_tree");
  const char* SEL = "CC1mu1piXp_Selected && CC1mu1piXp_MC_Signal";
  const char* P = "CC1mu1piXp_";

  // ---- muon momentum, relative, BOTH MCS and range estimators ----
  one_panel(t, SEL, Form("%scandidate_muon_mom_true",P),
            Form("(%scandidate_muon_mom_mcs-%scandidate_muon_mom_true)/%scandidate_muon_mom_true",P,P,P),
            "p_{#mu} resolution (MCS)", "true p_{#mu} [GeV/c]", "(reco#minustrue)/true",
            24, 0.15, 2.0, -1.0, 1.0, "resolution_pmu_mcs");
  one_panel(t, SEL, Form("%scandidate_muon_mom_true",P),
            Form("(%scandidate_muon_mom_range-%scandidate_muon_mom_true)/%scandidate_muon_mom_true",P,P,P),
            "p_{#mu} resolution (range)", "true p_{#mu} [GeV/c]", "(reco#minustrue)/true",
            24, 0.15, 2.0, -1.0, 1.0, "resolution_pmu_range");

  // ---- pion momentum, relative ----
  one_panel(t, SEL, Form("%scandidate_pion_mom_true",P),
            Form("(%scandidate_pion_mom_reco-%scandidate_pion_mom_true)/%scandidate_pion_mom_true",P,P,P),
            "p_{#pi} resolution", "true p_{#pi} [GeV/c]", "(reco#minustrue)/true",
            18, 0.175, 1.0, -1.0, 1.0, "resolution_ppi");

  // ---- muon cos(theta), absolute ----
  one_panel(t, SEL, Form("%scandidate_muon_costh_true",P),
            Form("(%scandidate_muon_costh_reco-%scandidate_muon_costh_true)",P,P),
            "cos#theta_{#mu} resolution", "true cos#theta_{#mu}", "reco#minustrue",
            20, -1.0, 1.0, -0.5, 0.5, "resolution_costhmu");

  // ---- pion cos(theta), absolute ----
  one_panel(t, SEL, Form("%scandidate_pion_costh_true",P),
            Form("(%scandidate_pion_costh_reco-%scandidate_pion_costh_true)",P,P),
            "cos#theta_{#pi} resolution", "true cos#theta_{#pi}", "reco#minustrue",
            20, -1.0, 1.0, -0.5, 0.5, "resolution_costhpi");

  // ---- opening angle, absolute ----
  one_panel(t, SEL, Form("%strue_mu_pi_opening_angle",P),
            Form("(%smu_pi_opening_angle-%strue_mu_pi_opening_angle)",P,P),
            "#theta_{#mu#pi} resolution", "true #theta_{#mu#pi} [rad]", "reco#minustrue [rad]",
            13, 0.0, 2.6, -1.0, 1.0, "resolution_thmupi");

  f->Close();
  printf("done -> report/figures/resolution_*.png\n");
}
