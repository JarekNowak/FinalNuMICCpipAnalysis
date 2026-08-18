// make_thetamu_vs_costhmu.C -- the muon-angle comparison figure: the SAME events and the
// SAME binning (theta edges are the arccos of the cos-theta edges), plotted once in
// cos(theta_mu) and once in theta_mu, so the only difference is the variable.
//
// There was no macro for this figure before: it existed only as a hand-made PNG carrying
// an FHC row and an RHC row. This version drops RHC and keeps the two FHC panels, which
// is the comparison the slide is actually making.
//
// Every curve is A_C-smeared (the unfolder dumps them that way into the closure files), so
// all of them live in the same space as the measurement.
//
//   usage: root -l -b -q report/make_thetamu_vs_costhmu.C

void make_thetamu_vs_costhmu() {

  gStyle->SetOptStat( 0 );
  gStyle->SetOptTitle( 1 );
  // Without this the pad frame picks up a dashed line style from the last histogram
  // drawn and a dashed black rectangle appears around the plot.
  gStyle->SetFrameLineStyle( 1 );

  const char* PROC = "/data/uboone/processed/";
  const char* obs[2]   = { "thetamu", "costhmu" };
  const char* title[2] = { "FHC #theta_{#mu}", "FHC cos#theta_{#mu}" };
  const char* xlab[2]  = { "#theta_{#mu} [rad]", "cos#theta_{#mu}" };

  // same Okabe-Ito palette and dash pattern as the dsigma montages
  const char* gens[4] = { "h_gen_GENIE", "h_gen_GiBUU", "h_gen_NEUT", "h_gen_NuWro" };
  const char* glab[4] = { "GENIE", "GiBUU", "NEUT", "NuWro" };
  int gcol[4] = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"),
                  TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00") };
  int gsty[4] = { 1, 2, 7, 9 };

  TCanvas c( "c_tvc", "", 1300, 470 );
  c.Divide( 2, 1 );
  std::vector<TObject*> keep;

  for ( int o = 0; o < 2; ++o ) {
    c.cd( o + 1 );
    gPad->SetBottomMargin( 0.14 ); gPad->SetLeftMargin( 0.13 ); gPad->SetTopMargin( 0.09 );

    TFile* f = TFile::Open( Form("%sclosure_hists_xsec_FHC5_%s.root", PROC, obs[o]) );
    if ( !f || f->IsZombie() ) { printf( "  missing closure for %s\n", obs[o] ); continue; }
    TH1D* hdat = (TH1D*) f->Get( "h_unfolded_nuwro" );
    TH1D* htun = (TH1D*) f->Get( "h_genie_tune" );
    if ( !hdat ) { f->Close(); continue; }
    hdat = (TH1D*) hdat->Clone(); hdat->SetDirectory(0); keep.push_back(hdat);
    if ( htun ) { htun = (TH1D*) htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
    std::vector<TH1D*> gh;
    for ( int g = 0; g < 4; ++g ) {
      TH1D* h = (TH1D*) f->Get( gens[g] );
      if ( h ) { h = (TH1D*) h->Clone(); h->SetDirectory(0); keep.push_back(h); }
      gh.push_back( h );
    }
    f->Close();

    double ymax = 0.;
    for ( int b = 1; b <= hdat->GetNbinsX(); ++b )
      ymax = std::max( ymax, hdat->GetBinContent(b) + hdat->GetBinError(b) );
    if ( htun ) for ( int b = 1; b <= htun->GetNbinsX(); ++b )
      ymax = std::max( ymax, htun->GetBinContent(b) );
    for ( auto h : gh ) if ( h ) for ( int b = 1; b <= h->GetNbinsX(); ++b )
      ymax = std::max( ymax, h->GetBinContent(b) );

    hdat->SetTitle( Form("%s;%s;d#sigma/dx [10^{-38} cm^{2}/Ar]", title[o], xlab[o]) );
    hdat->SetMarkerStyle( 20 ); hdat->SetMarkerSize( 1.0 );
    hdat->SetLineColor( kBlack ); hdat->SetLineWidth( 2 ); hdat->SetMarkerColor( kBlack );
    hdat->GetXaxis()->SetTitleSize( 0.050 ); hdat->GetXaxis()->SetLabelSize( 0.042 );
    hdat->GetYaxis()->SetTitleSize( 0.048 ); hdat->GetYaxis()->SetLabelSize( 0.042 );
    hdat->GetYaxis()->SetTitleOffset( 1.25 );
    hdat->SetMinimum( 0. ); hdat->SetMaximum( 1.45*ymax );
    hdat->Draw( "E1" );

    if ( htun ) {
      htun->SetLineColor( kBlack ); htun->SetLineWidth( 2 ); htun->SetMarkerSize( 0 );
      htun->Draw( "hist same" );
    }
    for ( int g = 0; g < 4; ++g ) {
      if ( !gh[g] ) continue;
      gh[g]->SetLineColor( gcol[g] ); gh[g]->SetLineStyle( gsty[g] );
      gh[g]->SetLineWidth( 2 ); gh[g]->SetMarkerSize( 0 );
      gh[g]->Draw( "hist same" );
    }
    hdat->Draw( "E1 same" );

    TLegend* lg = new TLegend( 0.40, 0.66, 0.93, 0.90 );
    lg->SetBorderSize( 0 ); lg->SetFillStyle( 0 ); lg->SetTextSize( 0.038 );
    lg->SetNColumns( 2 );
    lg->AddEntry( hdat, "unfolded data", "lep" );
    if ( htun ) lg->AddEntry( htun, "uB tune", "l" );
    for ( int g = 0; g < 4; ++g ) if ( gh[g] ) lg->AddEntry( gh[g], glab[g], "l" );
    lg->Draw();
    gPad->RedrawAxis();          // frame and ticks on top of the filled histograms
  }

  c.SaveAs( "report/figures/thetamu_vs_costhmu.eps" );
  gSystem->Exec( "report/eps2pdf.sh report/figures/thetamu_vs_costhmu.eps" );
  printf( "  wrote report/figures/thetamu_vs_costhmu.pdf\n" );
}
