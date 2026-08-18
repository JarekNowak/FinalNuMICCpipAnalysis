#pragma once

#include <vector>

// ROOT includes
#include "TDecompChol.h"
#include "TDecompSVD.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Unfolder.hh"

// Implementation of the Wiener-SVD unfolding method
// W. Tang et al., J. Instrum. 12, P10002 (2017)
// https://arxiv.org/abs/1705.03568
class WienerSVDUnfolder : public Unfolder {

  public:

    enum RegularizationMatrixType { kIdentity, kFirstDeriv, kSecondDeriv };

    inline WienerSVDUnfolder( bool use_wiener_filter = true,
      RegularizationMatrixType type = kIdentity ) : Unfolder(),
      use_filter_( use_wiener_filter ), reg_type_( type ) {}

    // Trick taken from https://stackoverflow.com/a/18100999
    using Unfolder::unfold;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const override;

    inline bool use_filter() const { return use_filter_; }
    inline void set_use_filter( bool use_filter ) { use_filter_ = use_filter; }

    inline RegularizationMatrixType get_regularization_type() const
      { return reg_type_; }

    inline void set_regularization_type( const RegularizationMatrixType& type )
      { reg_type_ = type; }

    // Physical widths of the true bins, in the order of the true-bin index.
    // OPTIONAL: if left empty the derivative regularisation matrices fall back to
    // assuming uniform spacing, which is what the code did unconditionally before.
    // Supplying them matters: the (1,-2,1) stencil is only a second derivative on a
    // uniform grid, and on the strongly non-uniform binnings used here it is
    // near-singular, which lets the Wiener filter retain a meaningless near-null
    // direction and collapses A_C towards rank one.
    inline void set_bin_widths( const std::vector< double >& w )
      { bin_widths_ = w; }

    inline const std::vector< double >& get_bin_widths() const
      { return bin_widths_; }

  protected:

    // Helper function that sets the contents of the regularization matrix
    // based on the current value of reg_type_
    void set_reg_matrix( TMatrixD& C ) const;

    // Flag indicating whether the Wiener filter should be used. If it is
    // false, then the usual expression will be replaced with an identity
    // matrix
    bool use_filter_ = true;

    // Enum that determines the form to use for the regularization matrix C
    RegularizationMatrixType reg_type_ = kIdentity;

    // True-bin widths; empty means "assume uniform" (see set_bin_widths)
    std::vector< double > bin_widths_;
};
