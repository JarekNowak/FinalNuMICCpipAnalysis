// dsigma_current.C — differential cross-section result figures for the CURRENT
// Poisson-thrown central-value fake data. Reads the UnfolderNuMI closure dumps
// closure_hists_xsec_<CFG>_<obs>.root (h_unfolded_nuwro = unfolded fake data,
// h_fakedata_truth = A_C-smeared truth, h_genie_tune = MicroBooNE tune), overlays
// the four generator FTE predictions, and writes one multi-panel figure per config.
//   usage: root -l -b -q 'macros/dsigma_current.C("FHC5","newg4")'  (or "RHCFULL","rhc" / "COMB","comb")
// 2026-09-27 (review 6.1-6.3): physics notation and units from dsigma_style.h, larger text, the
// legend in the free sixth cell instead of in every panel, open-ended top bins hatched and labelled.
#include "dsigma_style.h"

// Remap the p_pi histogram (open top bin) onto an equal-width axis. The adopted bins are
// [0.175,0.205] and everything above, i.e. 0.030 against 0.795 GeV/c: drawn to scale the
// first is 4% of the axis and collapses into an invisible spike that also sets the y-range,
// leaving the panel looking empty. Equal width keeps both bins readable; the cost is that
// area is no longer proportional to cross section, which the slide text states.
TH1D* eq2bin( TH1D* h ) {
  if ( !h ) return nullptr;
  static int uid = 0;
  int n = h->GetNbinsX();
  TH1D* nh = new TH1D( Form("eq2_%d", uid++), h->GetTitle(), n, 0., n );
  nh->SetDirectory( 0 );
  for ( int i = 1; i <= n; ++i ) {
    nh->SetBinContent( i, h->GetBinContent(i) );
    nh->SetBinError  ( i, h->GetBinError(i)   );
  }
  nh->SetLineColor( h->GetLineColor() );   nh->SetLineStyle( h->GetLineStyle() );
  nh->SetLineWidth( h->GetLineWidth() );   nh->SetMarkerColor( h->GetMarkerColor() );
  nh->SetMarkerStyle( h->GetMarkerStyle() ); nh->SetMarkerSize( h->GetMarkerSize() );
  // labels from the bin edges; the last p_pi bin is open above (2026-09-26: three bins)
  for ( int i = 1; i <= n; ++i )
    nh->GetXaxis()->SetBinLabel( i, i < n ? Form( "%g-%g", h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i) )
                                          : Form( "> %g", h->GetXaxis()->GetBinLowEdge(i) ) );
  nh->GetXaxis()->SetLabelSize( 0.058 );
  return nh;
}

void dsigma_current(const char* cfg = "FHC5", const char* gtag = "newg4") {
  const char* PROC = "/data/uboone/processed/";
  const char* GP   = "../generator_predictions/newg4/";
  // 2026-09-26: binnings of the 0.50 migration criterion. p_pi is read from the three-bin
  // extraction (ppi3bin); theta_mu is no longer reported (cos theta_mu carries the muon angle).
  const int NOBS = 5;
  const char* obs[NOBS]    = {"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* src[NOBS]    = {"pmu","ppi3bin","costhmu","costhpi","thmupi"};
  // panel titles and axis labels (physics notation, units) come from ds_obs() in dsigma_style.h
  // Okabe-Ito colorblind-safe palette + distinct line styles (redundant encoding,
  // so the four generators are separable in grayscale and for all colour-vision types)
  int gcol[4]   = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"),
                    TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00") }; // blue,green,purple,vermillion
  int gstyle[4] = { 1, 2, 7, 9 };
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  std::vector<TObject*> keep;

  std::vector<DsPanel> P(NOBS);
  for (int o = 0; o < NOBS; ++o) {
    P[o].key = src[o];
    TFile* f = TFile::Open(Form("%sclosure_hists_xsec_%s_%s.root", PROC, cfg, src[o]));
    if (!f || f->IsZombie()) { printf("  missing closure %s %s\n", cfg, obs[o]); continue; }
    TH1D* hunf = (TH1D*)f->Get("h_unfolded_nuwro");   // unfolded fake data
    TH1D* htru = (TH1D*)f->Get("h_fakedata_truth");   // A_C-smeared truth
    TH1D* htun = (TH1D*)f->Get("h_genie_tune");       // uB tune
    if (!hunf) continue;
    hunf = (TH1D*)hunf->Clone(); hunf->SetDirectory(0);
    hunf->SetMarkerStyle(20); hunf->SetMarkerSize(1.1); hunf->SetLineColor(kBlack); hunf->SetMarkerColor(kBlack);
    // Generator overlays: use the A_C-SMEARED predictions the unfolder dumped into the
    // closure file (h_gen_<Label>), NOT the raw truth-level FTE files. The unfolded data
    // and the uB tune both live in the A_C-smeared measurement space (Wiener-SVD returns
    // A_C * truth), so every model must be smeared by the same additional-smearing matrix
    // A_C for a fair comparison; otherwise the raw generators sit artificially high
    // relative to the A_C-smeared tune. These h_gen_* are already in dsigma/dx units, so
    // they are drawn directly (no bin-width division). (void)GP; (void)gtag; kept for API.
    const char* genHist[4] = { "h_gen_GENIE", "h_gen_GiBUU", "h_gen_NEUT", "h_gen_NuWro" };
    (void)GP; (void)gtag;
    std::vector<TH1D*> gh; std::vector<int> gidx;
    for (int g = 0; g < 4; ++g) {
      TH1D* hgd = (TH1D*)f->Get(genHist[g]);
      if (hgd && hgd->GetNbinsX()==hunf->GetNbinsX()) {
        TH1D* hg = (TH1D*)hgd->Clone(Form("g_%s_%d_%d",cfg,o,g)); hg->SetDirectory(0);
        hg->SetLineColor(gcol[g]); hg->SetLineStyle(gstyle[g]); hg->SetLineWidth(2); hg->SetMarkerSize(0);
        gh.push_back(hg); gidx.push_back(g); keep.push_back(hg);
      } else {
        printf("  [warn] %s %s: %s missing or bin-mismatch -> DROPPED\n", cfg,obs[o],genHist[g]);
      }
    }
    if (htru) { htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); keep.push_back(htru); }
    if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
    if (gh.empty()) P[o].gen_note = "no generator predictions";   // (every inclusive sidecar carries them)
    // open-ended top bin: its physical lower edge, from the result file (before any remapping)
    if ( ds_obs(src[o])->open_top ) {
      P[o].open_lo = hunf->GetXaxis()->GetBinLowEdge(hunf->GetNbinsX());
      P[o].open_hi = hunf->GetXaxis()->GetXmax();
    }
    if ( std::string(src[o]) == "ppi3bin" ) {
      hunf = eq2bin(hunf); if (htru) htru = eq2bin(htru); if (htun) htun = eq2bin(htun);
      for (auto& hg : gh) hg = eq2bin(hg);
      keep.push_back(hunf); if (htru) keep.push_back(htru); if (htun) keep.push_back(htun);
      for (auto hg : gh) if (hg) keep.push_back(hg);
      P[o].index_axis = true;
    }

    // y-axis max over data(+error), truth, tune and every generator curve
    double ymax = 0.;
    for (int b=1;b<=hunf->GetNbinsX();++b) ymax = std::max(ymax, hunf->GetBinContent(b)+hunf->GetBinError(b));
    if (htru) for (int b=1;b<=htru->GetNbinsX();++b) ymax = std::max(ymax, htru->GetBinContent(b));
    if (htun) for (int b=1;b<=htun->GetNbinsX();++b) ymax = std::max(ymax, htun->GetBinContent(b));
    for (auto hg : gh) for (int b=1;b<=hg->GetNbinsX();++b) ymax = std::max(ymax, hg->GetBinContent(b));
    hunf->SetMinimum(0); hunf->SetMaximum(1.35*ymax);
    if (htru) { htru->SetLineColor(kGray+2); htru->SetLineWidth(2); htru->SetLineStyle(2); }
    if (htun) { htun->SetLineColor(kBlack); htun->SetLineWidth(2); }
    keep.push_back(hunf);
    P[o].data = hunf; P[o].truth = htru; P[o].tune = htun; P[o].gen = gh; P[o].gidx = gidx;
  }
  std::vector<DsPanel*> all;
  for (auto& p : P) if (p.data) all.push_back(&p);

  // text in pixels for a 1600 px wide figure printed at 0.98 of the text width (~8.5 pt titles);
  // 560 px tall panels leave room for the full vertical-axis title with its units
  const DsStyle S = { /*title*/30, /*label*/27, /*binlabel*/24, /*legend*/30, /*ptitle*/33, /*note*/26,
                      /*lm*/0.215, /*rm*/0.045, /*tm*/0.095, /*bm*/0.16, /*xoff*/1.3, /*yoff*/1.85 };
  TString out = Form("unfold_output/dsigma_%s.pdf", cfg);
  ds_figure(out.Data(), all, 3, 2, 1600, 560, 0, S, cfg, keep);
}
