// raw_genie_tune.C — print the framework's RAW (pre-A_C) GENIE tune (and NuWro
// prediction if present) integrated cross section per observable, in
// 10^-38 cm^2/Ar, for absolute-scale comparison against standalone generators.
//
//   raw_genie_tune XSEC_Config SLICE_Config

#include <iostream>
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

int main( int argc, char* argv[] ) {
  if ( argc != 3 ) {
    std::cout << "Usage: raw_genie_tune XSEC_Config SLICE_Config\n";
    return 1;
  }
  CrossSectionExtractor extr( argv[1] );
  auto xsec = extr.get_unfolded_events();
  double conv_factor = extr.conversion_factor();
  const auto& pred_map = extr.get_prediction_map();

  SliceBinning sb( argv[2] );

  for ( size_t sl = 0u; sl < sb.slices_.size(); ++sl ) {
    const auto& slice = sb.slices_.at( sl );

    // width of any "other" (integrated-over) variables
    double other_w = 1.;
    for ( const auto& ov : slice.other_vars_ ) {
      double hi = ov.high_bin_edge_, lo = ov.low_bin_edge_;
      if ( hi != lo && std::abs(hi-lo) < 1e300 ) other_w *= (hi-lo);
    }

    for ( const auto& pp : pred_map ) {
      const std::string& name = pp.first;
      // raw prediction event counts -> slice histogram
      SliceHistogram* sh = SliceHistogram::make_slice_histogram(
        pp.second->get_prediction(), slice, nullptr );
      int nb = sh->hist_->GetNbinsX();
      // integrate dsigma/dx over the slice: sum( content/(conv*width*other_w) )
      double integ = 0.;
      for ( int b = 1; b <= nb; ++b ) {
        double w = sh->hist_->GetBinWidth( b ) * other_w;
        integ += sh->hist_->GetBinContent( b ) / ( conv_factor * w ) * w; // = content/conv/other_w summed
      }
      // simpler: total sigma = sum(content)/conv_factor/other_w  [in 1e-38]
      double total_sigma = 0.;
      for ( int b = 1; b <= nb; ++b ) total_sigma += sh->hist_->GetBinContent( b );
      total_sigma = total_sigma / conv_factor / other_w;
      std::cout << "slice " << sl << "  pred \"" << name << "\"  RAW total sigma = "
        << total_sigma << " [1e-38 cm^2/Ar]\n";
      // per-bin dsigma/dx (content / conv_factor / binwidth / other_w)
      std::cout << "  perbin[" << name << "]:";
      for ( int b = 1; b <= nb; ++b ) {
        double dsdx = sh->hist_->GetBinContent( b )
          / conv_factor / ( sh->hist_->GetBinWidth( b ) * other_w );
        std::cout << " " << dsdx;
      }
      std::cout << "\n";
      delete sh;
    }
  }
  return 0;
}
