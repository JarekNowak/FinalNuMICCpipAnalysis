// dump_xsec.C — dump the framework differential cross section dsigma/dx per bin
// (value + total uncertainty) and the integrated sigma, for building the note's
// results tables on the optimised binning.
//   dump_xsec XSEC_Config OBSNAME
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include "XSecAnalyzer/CrossSectionExtractor.hh"

int main( int argc, char* argv[] ) {
  if ( argc != 3 ) { std::cout << "Usage: dump_xsec XSEC_Config OBSNAME\n"; return 1; }
  std::string obs = argv[2];
  std::map<std::string, std::vector<double>> E = {
    {"pmu",     {0.150,0.350,0.550,0.750,0.950,1.250,1.750,3.000}},
    {"ppi",     {0.175,0.250,0.320,0.420,0.550,1.000}},
    {"costhmu", {-1.000,0.450,0.650,0.800,0.900,1.000}},
    {"costhpi", {-1.000,-0.100,0.350,0.550,0.750,1.000}},
    {"thmupi",  {0.000,0.600,0.850,1.100,1.300,1.520,1.850,2.600}},
  };
  auto edges = E.at(obs);

  CrossSectionExtractor extr( argv[1] );
  double conv = extr.conversion_factor();
  auto xsec = extr.get_unfolded_events();
  const TMatrixD& u   = *xsec.result_.unfolded_signal_;
  const TMatrixD& cov = *xsec.result_.cov_matrix_;
  int n = u.GetNrows();
  if ( (int)edges.size() != n+1 ) {
    std::cerr << "edge/bin mismatch: " << edges.size()-1 << " bins vs " << n << "\n"; return 1;
  }
  double sig_int = 0.;
  std::cout << "OBS " << obs << "\n";
  for ( int i = 0; i < n; ++i ) {
    double w    = edges[i+1] - edges[i];
    double dsdx = u(i,0) / ( conv * w );
    double unc  = std::sqrt( std::max(0.0, cov(i,i)) ) / ( conv * w );
    sig_int += u(i,0) / conv;
    printf( "BIN %d %.3f %.3f %.4f %.4f\n", i+1, edges[i], edges[i+1], dsdx, unc );
  }
  printf( "INT %.4f\n", sig_int );
  return 0;
}
