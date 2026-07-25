// ac_diagnostic.C — check the Wiener-SVD additional smearing matrix A_C and
// whether the unfolded estimator satisfies x_hat = A_C * x_true (as Wiener-SVD
// requires). Reports A_C's column sums / total gain, and compares x_hat to
// A_C * x_true (the fake-data truth) bin by bin.
//
//   ac_diagnostic XSEC_Config   (config must declare a "Fakedata" prediction)
#include <iostream>
#include <cmath>
#include "XSecAnalyzer/CrossSectionExtractor.hh"

int main( int argc, char* argv[] ) {
  if ( argc != 2 ) { std::cout << "Usage: ac_diagnostic XSEC_Config\n"; return 1; }
  CrossSectionExtractor extr( argv[1] );
  auto xsec = extr.get_unfolded_events();
  const TMatrixD& A_C  = *xsec.result_.add_smear_matrix_;
  const TMatrixD& xhat = *xsec.result_.unfolded_signal_;
  const auto& pred_map = extr.get_prediction_map();
  auto it = pred_map.find( "Fakedata" );
  if ( it == pred_map.end() ) { std::cout << "no \"Fakedata\" prediction in config\n"; return 1; }
  TMatrixD truth = it->second->get_prediction();

  int n = A_C.GetNrows();
  std::cout << "dims: A_C " << n << "x" << A_C.GetNcols()
            << "  xhat " << xhat.GetNrows() << "  truth " << truth.GetNrows() << "\n";

  double csmin=1e30, csmax=-1e30, cssum=0, diagsum=0;
  for ( int j=0; j<n; ++j ) { double cs=0; for ( int i=0; i<n; ++i ) cs += A_C(i,j);
    csmin=std::min(csmin,cs); csmax=std::max(csmax,cs); cssum+=cs; }
  for ( int i=0; i<n; ++i ) diagsum += A_C(i,i);
  std::cout << "A_C column sums: min=" << csmin << " max=" << csmax << " mean=" << cssum/n
            << "   (=1 would conserve the total)\n";
  std::cout << "A_C diagonal mean = " << diagsum/n << "   (=1 would be identity)\n";

  TMatrixD Ac_truth( A_C, TMatrixD::kMult, truth );
  double sx=0, st=0, sat=0;
  for ( int i=0; i<n; ++i ) { sx+=xhat(i,0); st+=truth(i,0); sat+=Ac_truth(i,0); }
  std::cout << "\nsum xhat=" << sx << "  sum truth=" << st << "  sum A_C*truth=" << sat << "\n";
  std::cout << "  xhat / truth        = " << sx/st  << "   <- the recurring ratio\n";
  std::cout << "  (A_C*truth) / truth = " << sat/st << "   <- what A_C alone does\n";
  std::cout << "  xhat / (A_C*truth)  = " << sx/sat << "   <- should be ~1 if x_hat = A_C*x_true\n";
  std::cout << "\nbin |   xhat   |  truth  | A_C*truth | xhat/(A_C*truth)\n";
  for ( int i=0; i<n; ++i )
    std::cout << "  " << i << " | " << xhat(i,0) << " | " << truth(i,0) << " | "
              << Ac_truth(i,0) << " | " << ( Ac_truth(i,0)!=0 ? xhat(i,0)/Ac_truth(i,0) : 0 ) << "\n";
  return 0;
}
