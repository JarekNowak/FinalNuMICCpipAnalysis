// syst_breakdown.C — fractional systematic-uncertainty breakdown by source.
// Reads the per-source covariance matrices the framework builds
// (xsec.unfolded_cov_matrix_map_) and the unfolded signal, and plots the
// per-bin fractional uncertainty for the total and for each source group
// (cross section, flux/beamline, detector, reinteraction, stats, normalisation).
// Fractional uncertainty is unit-invariant, so event-count units are used.
//
//   syst_breakdown XSEC_Config OBSNAME
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

int main( int argc, char* argv[] ) {
  if ( argc != 3 ) { std::cout << "Usage: syst_breakdown XSEC_Config OBSNAME\n"; return 1; }
  CrossSectionExtractor extr( argv[1] );
  auto xsec = extr.get_unfolded_events();
  const TMatrixD& u = *xsec.result_.unfolded_signal_;
  int N = u.GetNrows();
  auto& m = xsec.unfolded_cov_matrix_map_;

  struct Src { std::string label; std::vector<std::string> keys; int color; int width; };
  std::vector<Src> srcs = {
    { "Total",          { "total" },                       kBlack,     3 },
    { "Cross section",  { "xsec_total" },                  kRed+1,     2 },
    { "Flux (PPFX)",    { "flux" },                        kBlue+1,    2 },
    { "Beamline geom.", { "flux_beamline" },               kCyan+2,    2 },
    { "Detector",       { "detVar_total" },                kGreen+2,   2 },
    { "Reinteraction",  { "reint" },                       kMagenta+1, 2 },
    { "MC + data stats",{ "SimulationStats", "DataStats" },kOrange+7,  2 },
    { "POT + targets",  { "POT", "numTargets" },           kGray+2,    2 },
  };

  std::vector<TH1D*> H;
  for ( size_t s = 0; s < srcs.size(); ++s ) {
    TH1D* h = new TH1D( Form("h_%zu",s), "", N, 0.5, N+0.5 );
    for ( int i = 0; i < N; ++i ) {
      double var = 0.;
      for ( const auto& k : srcs[s].keys ) if ( m.count(k) ) var += (*m[k])(i,i);
      double frac = ( u(i,0) > 0. ) ? 100.0*std::sqrt(var)/u(i,0) : 0.;
      h->SetBinContent( i+1, frac );
    }
    h->SetLineColor( srcs[s].color ); h->SetLineWidth( srcs[s].width );
    H.push_back( h );
  }

  gStyle->SetOptStat(0);
  TCanvas c( "c", "", 900, 600 ); c.SetGridy(); c.SetLeftMargin(0.11);
  double ymax = 0.; for ( auto h : H ) ymax = std::max( ymax, h->GetMaximum() );
  H[0]->SetTitle( Form("Fractional systematic uncertainty: %s;True bin;Fractional uncertainty [%%]", argv[2]) );
  H[0]->GetYaxis()->SetRangeUser( 0., ymax*1.30 );
  H[0]->Draw( "hist" );
  for ( size_t i = 1; i < H.size(); ++i ) H[i]->Draw( "hist same" );
  TLegend leg( 0.50, 0.60, 0.88, 0.88 ); leg.SetBorderSize(0); leg.SetFillStyle(0);
  for ( size_t i = 0; i < H.size(); ++i ) leg.AddEntry( H[i], srcs[i].label.c_str(), "l" );
  leg.Draw();
  c.SaveAs( Form("unfold_output/syst_breakdown_%s.pdf", argv[2]) );

  std::cout << "=== " << argv[2] << " bin-averaged fractional uncertainty ===\n";
  for ( size_t i = 0; i < H.size(); ++i ) {
    double avg = 0.; for ( int b = 1; b <= N; ++b ) avg += H[i]->GetBinContent(b); avg /= N;
    printf( "  %-16s %5.1f%%\n", srcs[i].label.c_str(), avg );
  }
  return 0;
}
