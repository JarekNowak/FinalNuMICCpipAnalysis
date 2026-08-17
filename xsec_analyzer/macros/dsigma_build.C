// dsigma_build.C — progressive-reveal ("build") version of the inclusive differential
// cross-section montage, for talks. Produces one figure per build step, each adding one
// more curve on identical axes so the slides overlay cleanly:
//
//   step 1 : A_C-smeared fake-data truth          (what the measurement should recover)
//   step 2 : + unfolded fake data                 (the closure)
//   step 3 : + MicroBooNE tune
//   step 4 : + GENIE
//   step 5 : + GiBUU
//   step 6 : + NEUT
//   step 7 : + NuWro
//
// The y-axis range is fixed across all steps from the FULL set of curves, so nothing
// jumps between slides.
//   usage: root -l -b -q 'macros/dsigma_build.C("FHC5")'
#include <vector>

// Remap a two-bin p_pi histogram onto an equal-width axis. The adopted bins are
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
  nh->GetXaxis()->SetBinLabel( 1, "0.175-0.205" );
  nh->GetXaxis()->SetBinLabel( 2, "> 0.205" );
  nh->GetXaxis()->SetLabelSize( 0.058 );
  return nh;
}

void dsigma_build(const char* cfg = "FHC5") {
  const char* PROC = "/data/uboone/processed/";
  const char* obs[6]  = {"pmu","ppi","costhmu","costhpi","thmupi","thetamu"};
  // p_pi is read from the ADOPTED two-bin extraction; the five-bin version is
  // superseded (only its lowest bin cleared the 0.68 smearing-diagonal criterion).
  const char* src[6]    = {"pmu","ppi2bin","costhmu","costhpi","thmupi","thetamu"};
  const char* obsX[6] = {"p_{#mu} [GeV/c]","p_{#pi} [GeV/c]","cos#theta_{#mu}",
                         "cos#theta_{#pi}","#theta_{#mu#pi} [rad]","#theta_{#mu} [rad]"};
  const char* gens[4]   = {"h_gen_GENIE","h_gen_GiBUU","h_gen_NEUT","h_gen_NuWro"};
  const char* glab[4]   = {"GENIE","GiBUU","NEUT","NuWro"};
  int gcol[4] = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"),
                  TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00") };
  int gsty[4] = { 1, 2, 7, 9 };
  gStyle->SetOptStat(0);

  const int NSTEP = 7;
  for (int step = 1; step <= NSTEP; ++step) {
    TCanvas c(Form("build_%s_%d",cfg,step), "", 1600, 900);
    c.Divide(3, 2);
    std::vector<TObject*> keep;

  // If the x-axis does not start at zero, ROOT only labels the first "round" tick, so the
  // bottom-left corner reads as the origin (0,0) even though it is not. Draw the actual
  // lower limit under the left edge of the frame.
  auto mark_xmin = [](TH1* h){
    double xmin = h->GetXaxis()->GetXmin();
    if ( xmin <= 0. ) return;
    TLatex* t = new TLatex();
    t->SetNDC(); t->SetTextSize(0.045); t->SetTextAlign(23);
    t->DrawLatex( gPad->GetLeftMargin(), gPad->GetBottomMargin()-0.055,
                  Form("%g", xmin) );
  };

    for (int o = 0; o < 6; ++o) {
      c.cd(o+1);
      gPad->SetBottomMargin(0.15); gPad->SetLeftMargin(0.16); gPad->SetTopMargin(0.08);
      TFile* f = TFile::Open(Form("%sclosure_hists_xsec_%s_%s.root", PROC, cfg, src[o]));
      if (!f || f->IsZombie()) { printf("  missing %s %s\n",cfg,obs[o]); continue; }
      TH1D* hdat = (TH1D*)f->Get("h_unfolded_nuwro");
      TH1D* htru = (TH1D*)f->Get("h_fakedata_truth");
      TH1D* htun = (TH1D*)f->Get("h_genie_tune");
      if (!hdat || !htru) { f->Close(); continue; }
      hdat=(TH1D*)hdat->Clone(); hdat->SetDirectory(0); keep.push_back(hdat);
      htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); keep.push_back(htru);
      if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
      std::vector<TH1D*> gh;
      for (int g = 0; g < 4; ++g) {
        TH1D* h = (TH1D*)f->Get(gens[g]);
        if (!h) { gh.push_back(nullptr); continue; }
        h=(TH1D*)h->Clone(); h->SetDirectory(0);
        h->SetLineColor(gcol[g]); h->SetLineStyle(gsty[g]); h->SetLineWidth(2);
        h->SetMarkerSize(0); gh.push_back(h); keep.push_back(h);
      }
      f->Close();

      if ( std::string(src[o]) == "ppi2bin" ) {
        hdat = eq2bin(hdat); htru = eq2bin(htru); if (htun) htun = eq2bin(htun);
        for (auto& h : gh) h = eq2bin(h);
        keep.push_back(hdat); keep.push_back(htru); if (htun) keep.push_back(htun);
        for (auto h : gh) if (h) keep.push_back(h);
      }

      // Fixed y-range from ALL curves so the axes never move between build steps
      double ymax = 0.;
      for (int b=1;b<=hdat->GetNbinsX();++b) ymax=std::max(ymax,hdat->GetBinContent(b)+hdat->GetBinError(b));
      for (int b=1;b<=htru->GetNbinsX();++b) ymax=std::max(ymax,htru->GetBinContent(b));
      if (htun) for (int b=1;b<=htun->GetNbinsX();++b) ymax=std::max(ymax,htun->GetBinContent(b));
      for (auto h : gh) if (h) for (int b=1;b<=h->GetNbinsX();++b) ymax=std::max(ymax,h->GetBinContent(b));

      // frame: always draw the truth first (step 1 shows it alone)
      htru->SetTitle(Form("%s;%s;d#sigma/dx [10^{-38} cm^{2}/Ar]", obs[o], obsX[o]));
      htru->SetLineColor(kGray+2); htru->SetLineWidth(3); htru->SetLineStyle(2);
      htru->SetMarkerSize(0);
      htru->GetXaxis()->SetTitleSize(0.050); htru->GetXaxis()->SetLabelSize(0.042);
      htru->GetXaxis()->SetTitleOffset(1.25);
      htru->GetYaxis()->SetTitleSize(0.048); htru->GetYaxis()->SetLabelSize(0.042);
      htru->GetYaxis()->SetTitleOffset(1.55);
      htru->SetMinimum(0); htru->SetMaximum(1.45*ymax);
      htru->Draw("hist");

      bool fwd = (o==2 || o==3);   // forward-peaked cos(theta): legend on the left
      TLegend* lg = new TLegend(fwd?0.18:0.42, 0.62, fwd?0.62:0.90, 0.90);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030); lg->SetNColumns(2);
      lg->AddEntry(htru, "Truth (A_{C})", "l");

      if (step >= 3 && htun) { htun->SetLineColor(kBlack); htun->SetLineWidth(2);
        htun->Draw("hist same"); lg->AddEntry(htun, "uB tune", "l"); }
      for (int g = 0; g < 4; ++g)
        if (step >= 4+g && gh[g]) { gh[g]->Draw("hist same"); lg->AddEntry(gh[g], glab[g], "l"); }

      if (step >= 2) {
        hdat->SetMarkerStyle(20); hdat->SetMarkerSize(0.9);
        hdat->SetLineColor(kBlack); hdat->SetMarkerColor(kBlack);
        hdat->Draw("E1 same");
        lg->AddEntry(hdat, "Unfolded fake data", "lep");
      }
      mark_xmin(htru);
    lg->Draw(); keep.push_back(lg);
    }
    TString out = Form("unfold_output/dsigma_build_%s_%d.pdf", cfg, step);
    c.SaveAs(out);
    printf("wrote %s\n", out.Data());
  }
}
