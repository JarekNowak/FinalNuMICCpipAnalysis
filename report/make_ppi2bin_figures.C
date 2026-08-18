// make_ppi2bin_figures.C -- result figures for the two-bin p_pi cross section.
//
// Reads the closure files written by UnfolderNuMI for the three beam configurations and
// draws d(sigma)/dp_pi: unfolded fake data against the A_C-smeared truth, the uB tune and
// the four generators. Every curve is A_C-smeared, so all of them live in the same space
// as the measurement and the comparison is like-for-like.
//
// PRESENTATION NOTE. The two bins are [0.175,0.205] and [0.205,1.000] GeV/c, i.e. 0.030
// and 0.795 wide. Drawn to scale the first bin is 4% of the axis and unreadable, so the
// panels use EQUAL-WIDTH bins with the true edges written on the axis. The consequence is
// that area on these panels is not proportional to cross section; the integrated sigma is
// printed on each panel instead, and the caption says so.
//
//   usage: root -l -b -q report/make_ppi2bin_figures.C

#include <vector>
#include <string>

static const char* FIG = "report/figures/";

// generator curves, in the drawing order used by the note table
static const char* GEN[4]  = { "h_gen_GENIE", "h_gen_GiBUU", "h_gen_NEUT", "h_gen_NuWro" };
static const char* GLAB[4] = { "GENIE", "GiBUU", "NEUT", "NuWro" };
static const int   GCOL[4] = { kGreen+2, kMagenta+1, kOrange+7, kAzure+2 };
// distinct dash patterns as well as colours: the four generators sit within ~0.1 of each
// other in bin 1, so colour alone does not separate them in print or on a projector
static const int   GSTY[4] = { 9, 7, 5, 10 };

// integrated sigma quoted on each panel (10^-38 cm^2/Ar), from the note table
struct Cfg { const char* tag; const char* file; const char* label; double sig; double chi2; double pval; };
static const Cfg CFG[3] = {
  { "FHC",  "FHC5",    "FHC (#nu_{#mu})",       0.976, 0.15, 0.93 },
  { "RHC",  "RHCFULL", "RHC (#bar{#nu}_{#mu})", 0.887, 0.24, 0.89 },
  { "COMB", "COMB",    "Combined",              0.903, 0.13, 0.94 }
};

// Re-bin a 2-bin histogram onto an equal-width axis so the narrow first bin stays visible.
TH1D* equalise( TH1D* src, const char* name ) {
  TH1D* h = new TH1D( name, "", 2, 0., 2. );
  h->SetDirectory( 0 );
  for ( int i = 1; i <= 2; ++i ) {
    h->SetBinContent( i, src->GetBinContent(i) );
    h->SetBinError  ( i, src->GetBinError(i)   );
  }
  h->GetXaxis()->SetBinLabel( 1, "0.175 - 0.205" );
  h->GetXaxis()->SetBinLabel( 2, "> 0.205" );
  h->GetXaxis()->SetLabelSize( 0.062 );
  return h;
}

void make_ppi2bin_figures() {

  gStyle->SetOptStat( 0 );
  gStyle->SetOptTitle( 0 );
  gStyle->SetEndErrorSize( 6 );

  // Canvases are saved as EPS and converted to tight PDFs at the end of this macro; see
  // report/eps2pdf.sh for why the direct-to-PDF route gives an A4 page with the plot
  // shrunk into a corner.

  for ( int c = 0; c < 3; ++c ) {

    TFile* f = TFile::Open(
      Form("/data/uboone/processed/closure_hists_xsec_%s_ppi2bin.root", CFG[c].file) );
    if ( !f || f->IsZombie() ) { printf( "  [skip] %s\n", CFG[c].file ); continue; }

    // Clone everything off the input file straight away and detach it. Writing or drawing
    // while an input TFile is the active directory has silently sent histograms into the
    // wrong file before; detaching first makes that impossible.
    TH1D* hdat = equalise( (TH1D*) f->Get("h_unfolded_nuwro"),  Form("dat_%s",CFG[c].tag) );
    TH1D* htru = equalise( (TH1D*) f->Get("h_fakedata_truth"),  Form("tru_%s",CFG[c].tag) );
    TH1D* htun = equalise( (TH1D*) f->Get("h_genie_tune"),      Form("tun_%s",CFG[c].tag) );
    TH1D* hgen[4];
    for ( int g = 0; g < 4; ++g )
      hgen[g] = equalise( (TH1D*) f->Get(GEN[g]), Form("%s_%s",GLAB[g],CFG[c].tag) );
    f->Close();

    // ---- styles ------------------------------------------------------------------
    hdat->SetMarkerStyle( 20 ); hdat->SetMarkerSize( 1.3 );
    hdat->SetLineColor( kBlack ); hdat->SetLineWidth( 2 ); hdat->SetMarkerColor( kBlack );

    htru->SetLineColor( kRed+1 );  htru->SetLineWidth( 3 );
    htun->SetLineColor( kBlue+1 ); htun->SetLineWidth( 3 ); htun->SetLineStyle( 2 );
    for ( int g = 0; g < 4; ++g ) {
      hgen[g]->SetLineColor( GCOL[g] );
      hgen[g]->SetLineWidth( 2 );
      hgen[g]->SetLineStyle( GSTY[g] );
    }

    // ---- frame -------------------------------------------------------------------
    double ymax = hdat->GetBinContent(1) + hdat->GetBinError(1);
    for ( int g = 0; g < 4; ++g ) ymax = std::max( ymax, hgen[g]->GetBinContent(1) );
    ymax = std::max( ymax, htru->GetBinContent(1) );

    TCanvas cv( Form("c_%s",CFG[c].tag), "", 640, 560 );
    cv.SetLeftMargin( 0.16 ); cv.SetBottomMargin( 0.14 ); cv.SetTopMargin( 0.05 );

    hdat->SetTitle( ";p_{#pi} [GeV/c];"
                    "d#sigma/dp_{#pi} [10^{-38} cm^{2}/(GeV/c)/Ar]" );
    hdat->GetYaxis()->SetTitleSize( 0.052 ); hdat->GetYaxis()->SetTitleOffset( 1.30 );
    hdat->GetXaxis()->SetTitleSize( 0.052 );
    hdat->GetYaxis()->SetLabelSize( 0.048 );
    hdat->GetYaxis()->SetRangeUser( 0., 1.32*ymax );

    hdat->Draw( "E1" );
    htru->Draw( "hist same" );
    htun->Draw( "hist same" );
    for ( int g = 0; g < 4; ++g ) hgen[g]->Draw( "hist same" );
    hdat->Draw( "E1 same" );   // data on top

    // ---- legend ------------------------------------------------------------------
    // Legend on the right, panel stats on the left in short lines that cannot reach it.
    TLegend lg( 0.47, 0.53, 0.95, 0.94 );
    lg.SetBorderSize( 0 ); lg.SetFillStyle( 0 ); lg.SetTextSize( 0.042 );
    lg.AddEntry( hdat, "unfolded data", "lep" );
    lg.AddEntry( htru, "A_{C} truth", "l" );
    lg.AddEntry( htun, "uB tune", "l" );
    for ( int g = 0; g < 4; ++g ) lg.AddEntry( hgen[g], GLAB[g], "l" );
    lg.Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize( 0.050 );
    tx.DrawLatex( 0.18, 0.885, Form("#bf{%s}", CFG[c].label) );
    tx.SetTextSize( 0.042 );
    tx.DrawLatex( 0.18, 0.825, Form("#sigma_{int} = %.3f", CFG[c].sig) );
    tx.DrawLatex( 0.18, 0.768, Form("#chi^{2}/ndf = %.2f/2", CFG[c].chi2) );
    tx.DrawLatex( 0.18, 0.711, Form("p = %.2f", CFG[c].pval) );

    cv.SaveAs( Form("%sppi2bin_xsec_%s.eps", FIG, CFG[c].tag) );
    
  }

  // ---- the additional-smearing matrix, one panel per configuration ----------------
  for ( int c = 0; c < 3; ++c ) {
    TFile* f = TFile::Open(
      Form("/data/uboone/processed/closure_hists_xsec_%s_ppi2bin.root", CFG[c].file) );
    if ( !f || f->IsZombie() ) continue;
    TH2D* src = (TH2D*) f->Get( "h_A_C" );
    if ( !src ) { f->Close(); continue; }
    TH2D* h = (TH2D*) src->Clone( Form("ac_%s",CFG[c].tag) );
    h->SetDirectory( 0 );
    f->Close();

    TCanvas cv( Form("cac_%s",CFG[c].tag), "", 520, 460 );
    cv.SetLeftMargin( 0.26 ); cv.SetRightMargin( 0.16 ); cv.SetBottomMargin( 0.14 );
    gStyle->SetPaintTextFormat( "5.2f" );
    h->SetTitle( ";true bin;smeared bin" );
    h->SetMarkerSize( 2.4 );
    for ( int b = 1; b <= 2; ++b ) {
      h->GetXaxis()->SetBinLabel( b, b == 1 ? "0.175-0.205" : "> 0.205" );
      h->GetYaxis()->SetBinLabel( b, b == 1 ? "0.175-0.205" : "> 0.205" );
    }
    h->GetXaxis()->SetLabelSize( 0.052 ); h->GetYaxis()->SetLabelSize( 0.052 );
    h->GetYaxis()->SetTitleOffset( 2.4 );
    h->GetXaxis()->SetTitle( "true bin [GeV/c]" );
    h->GetYaxis()->SetTitle( "smeared bin [GeV/c]" );
    h->Draw( "colz" );
    h->Draw( "text same" );
    TLatex tx; tx.SetNDC(); tx.SetTextSize( 0.05 );
    tx.DrawLatex( 0.14, 0.945, Form("#bf{A_{C}, %s}", CFG[c].label) );
    cv.SaveAs( Form("%sppi2bin_AC_%s.eps", FIG, CFG[c].tag) );
    
  }


  // ---- proton-tagged (CC1mu1pi1p) two-bin p_pi ------------------------------------
  // Same scheme and same edges as the inclusive result, so the two are directly
  // comparable. These closure files carry no generator curves: the ccpi1p xsec configs
  // declare only the tune and the fake data. Configurations still running are skipped.
  {
    // tag, closure-file tag, panel label, sigma_int, chi2, p, diagonals -- the last four
    // are filled in from the unfold log as each configuration lands
    struct P { const char* tag; const char* file; const char* label;
               double sig; double chi2; double pval; const char* diag; };
    const P onep[3] = {
      { "FHC",  "ccpi1p_FHC5",    "Proton-tagged, FHC",      0.575, 0.15, 0.93, "89 / 74%" },
      { "RHC",  "ccpi1p_RHCFULL", "Proton-tagged, RHC",      0.495, 2.37, 0.31, "89 / 75%" },
      { "COMB", "ccpi1p_COMB",    "Proton-tagged, combined", 0.,    0.,   0.,   nullptr    }
    };

    for ( int c = 0; c < 3; ++c ) {
      TFile* f = TFile::Open(
        Form("/data/uboone/processed/closure_hists_xsec_%s_ppi2bin.root", onep[c].file) );
      if ( !f || f->IsZombie() ) {
        printf( "  [skip] proton-tagged %s (not yet unfolded)\n", onep[c].tag );
        continue;
      }
      TH1D* hdat = equalise( (TH1D*) f->Get("h_unfolded_nuwro"), Form("dat1p_%s",onep[c].tag) );
      TH1D* htru = equalise( (TH1D*) f->Get("h_fakedata_truth"), Form("tru1p_%s",onep[c].tag) );
      TH1D* htun = equalise( (TH1D*) f->Get("h_genie_tune"),     Form("tun1p_%s",onep[c].tag) );
      f->Close();
      if ( !hdat || !htru ) continue;

      hdat->SetMarkerStyle( 20 ); hdat->SetMarkerSize( 1.3 );
      hdat->SetLineColor( kBlack ); hdat->SetLineWidth( 2 ); hdat->SetMarkerColor( kBlack );
      htru->SetLineColor( kRed+1 );  htru->SetLineWidth( 3 );
      if ( htun ) { htun->SetLineColor( kBlue+1 ); htun->SetLineWidth( 3 ); htun->SetLineStyle( 2 ); }

      double ymax = hdat->GetBinContent(1) + hdat->GetBinError(1);
      ymax = std::max( ymax, htru->GetBinContent(1) );

      TCanvas cv( Form("c1p_%s",onep[c].tag), "", 640, 560 );
      cv.SetLeftMargin( 0.16 ); cv.SetBottomMargin( 0.14 ); cv.SetTopMargin( 0.05 );
      hdat->SetTitle( ";p_{#pi} [GeV/c];"
                      "d#sigma/dp_{#pi} [10^{-38} cm^{2}/(GeV/c)/Ar]" );
      hdat->GetYaxis()->SetTitleSize( 0.052 ); hdat->GetYaxis()->SetTitleOffset( 1.30 );
      hdat->GetXaxis()->SetTitleSize( 0.052 ); hdat->GetYaxis()->SetLabelSize( 0.048 );
      hdat->GetYaxis()->SetRangeUser( 0., 1.32*ymax );
      hdat->Draw( "E1" );
      htru->Draw( "hist same" );
      if ( htun ) htun->Draw( "hist same" );
      hdat->Draw( "E1 same" );

      TLegend lg( 0.58, 0.71, 0.96, 0.93 );
      lg.SetBorderSize( 0 ); lg.SetFillStyle( 0 ); lg.SetTextSize( 0.040 );
      lg.AddEntry( hdat, "unfolded data", "lep" );
      lg.AddEntry( htru, "A_{C} truth", "l" );
      if ( htun ) lg.AddEntry( htun, "uB tune", "l" );
      lg.Draw();

      TLatex tx; tx.SetNDC(); tx.SetTextSize( 0.044 );
      tx.DrawLatex( 0.19, 0.888, Form("#bf{%s}", onep[c].label) );
      tx.SetTextSize( 0.042 );
      if ( onep[c].sig > 0. ) {
        tx.DrawLatex( 0.19, 0.825, Form("#sigma_{int} = %.3f", onep[c].sig) );
        tx.DrawLatex( 0.19, 0.768, Form("#chi^{2}/ndf = %.2f/2", onep[c].chi2) );
        tx.DrawLatex( 0.19, 0.711, Form("p = %.2f", onep[c].pval) );
      }
      if ( onep[c].diag ) tx.DrawLatex( 0.19, 0.654, Form("diagonals %s", onep[c].diag) );

      cv.SaveAs( Form("%sppi2bin_xsec_1p_%s.eps", FIG, onep[c].tag) );
    }
  }

  // Convert to tight PDFs so that an \includegraphics width in the note sizes the plot
  // rather than the surrounding whitespace.
  gSystem->Exec( Form("report/eps2pdf.sh %sppi2bin_*.eps", FIG) );
}
