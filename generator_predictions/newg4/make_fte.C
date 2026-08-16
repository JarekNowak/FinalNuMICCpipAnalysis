// Convert a *_newg4_final.root (per-observable dsigma/dx histos) to the
// FileTrueEvents (FTE) format the xsec configs consume: a flat-index histogram
// <obs>_fte (nbins = observable bins, x-range 0..N) whose bin content is the
// per-bin cross section dsigma/dx * binwidth [1e-38 cm^2/Ar]. Matches how the
// NuWro/GENIE FTE files were built (verified against nuwro_newg4_fte.root).
//
// THE THETA OBSERVABLES. theta_mu / theta_pi are not separate inputs: their bin
// edges are exactly acos() of the corresponding cos(theta) edges, so the same
// phase space is partitioned identically and the per-bin cross section is the
// cos(theta) content in REVERSED bin order (acos is monotonically decreasing, so
// cos-bin 1 <-> theta-bin N). Sigma is conserved exactly. This reversal used to
// be done by hand and was the only step in the chain not reproducible from
// committed code, so regenerating the FTE files silently dropped thetamu_fte /
// thetapi_fte, which the theta xsec configs require.
//
// PREFER add_theta_fte() OVER make_fte() FOR THE THETA HISTOGRAMS. The
// *_newg4_final.root inputs are STALE with respect to the coarsened analysis
// binning (they still carry ppi 6 / costhpi 5 / thmupi 7 bins, against the
// 5 / 4 / 5 the analysis and the on-disk FTE files use). Running make_fte() on
// them would therefore REGRESS those three observables. add_theta_fte() derives
// the theta histograms from the cos(theta) histograms inside an existing FTE
// file, so it is correct regardless of the state of *_final.root -- it was
// verified to reproduce all twelve on-disk FTE files bit-for-bit.

#include <cmath>

namespace {
  // Write <dst> as the exact bin reversal of <src>, preserving per-bin cross
  // section and errors. Returns false if src is missing.
  bool reverse_into( TFile& fo, TH1* src, const char* dst_name ) {
    if ( !src ) return false;
    int nb = src->GetNbinsX();
    TH1D* h = new TH1D( dst_name, "", nb, 0, nb );
    double sum_src = 0., sum_dst = 0.;
    for ( int b = 1; b <= nb; ++b ) {
      h->SetBinContent( nb + 1 - b, src->GetBinContent(b) );
      h->SetBinError(   nb + 1 - b, src->GetBinError(b) );
      sum_src += src->GetBinContent(b);
    }
    for ( int b = 1; b <= nb; ++b ) sum_dst += h->GetBinContent(b);
    if ( std::abs(sum_src - sum_dst) > 1e-9 * std::abs(sum_src) ) {
      printf( "  [WARN] %s: sigma not conserved in reversal (%.6e -> %.6e)\n",
        dst_name, sum_src, sum_dst );
    }
    h->SetDirectory( &fo );
    h->Write( dst_name, TObject::kOverwrite );
    return true;
  }
}

// Add (or refresh) thetamu_fte / thetapi_fte inside an EXISTING FTE file by
// reversing its own costhmu_fte / costhpi_fte. This is the reproducible,
// binning-safe way to produce the theta predictions.
//   usage: root -l -b -q 'make_fte.C+' then add_theta_fte("genie_newg4_fte.root")
//      or: root -l -b -q 'make_fte.C("","")' is NOT needed -- call directly:
//          root -l -b -q 'add_theta_fte.C' style wrappers are unnecessary since
//          this file defines both functions.
void add_theta_fte( const char* ftefile ) {
  TFile fo( ftefile, "UPDATE" );
  if ( fo.IsZombie() ) { printf( "cannot open %s\n", ftefile ); return; }
  const char* src[2] = { "costhmu_fte", "costhpi_fte" };
  const char* dst[2] = { "thetamu_fte", "thetapi_fte" };
  int n_ok = 0;
  for ( int r = 0; r < 2; ++r ) {
    TH1* h = (TH1*) fo.Get( src[r] );
    if ( !h ) { printf( "  missing %s in %s -> cannot build %s\n",
                        src[r], ftefile, dst[r] ); continue; }
    if ( reverse_into( fo, h, dst[r] ) ) ++n_ok;
  }
  fo.Close();
  printf( "add_theta_fte: %s (%d/2 theta observables written)\n", ftefile, n_ok );
}

void make_fte( const char* infile, const char* outfile ) {
  const char* obs[5] = { "pmu", "ppi", "costhmu", "costhpi", "thmupi" };
  TFile fi( infile ), fo( outfile, "recreate" );
  if ( fi.IsZombie() ) { printf( "cannot open %s\n", infile ); return; }

  for ( int o = 0; o < 5; ++o ) {
    TH1D* h = (TH1D*) fi.Get( obs[o] );
    if ( !h ) { printf( "missing %s in %s\n", obs[o], infile ); continue; }
    int nb = h->GetNbinsX();
    TH1D* fte = new TH1D( Form("%s_fte", obs[o]), "", nb, 0, nb );
    for ( int b = 1; b <= nb; ++b ) {
      fte->SetBinContent( b, h->GetBinContent(b) * h->GetBinWidth(b) );
      fte->SetBinError(   b, h->GetBinError(b)   * h->GetBinWidth(b) );
    }
    fte->SetDirectory( &fo );
    fte->Write();
  }

  // Derive the theta observables from the cos(theta) histograms just written.
  for ( int r = 0; r < 2; ++r ) {
    const char* s = r == 0 ? "costhmu_fte" : "costhpi_fte";
    const char* d = r == 0 ? "thetamu_fte" : "thetapi_fte";
    reverse_into( fo, (TH1*) fo.Get(s), d );
  }
  fo.Close();

  // Staleness guard: the coarsened analysis binning is authoritative. If an FTE
  // file already exists alongside, compare bin counts and shout on a mismatch --
  // silently regressing ppi/costhpi/thmupi to the pre-coarsening binning is the
  // failure mode this chain is prone to.
  TString ref( outfile );
  if ( !ref.EndsWith("_fte.root") ) ref += "";
  TFile fchk( outfile, "READ" );
  const char* names[7] = { "pmu_fte", "ppi_fte", "costhmu_fte", "costhpi_fte",
                           "thmupi_fte", "thetamu_fte", "thetapi_fte" };
  const int expect[7] = { 7, 5, 5, 4, 5, 5, 4 };   // current coarsened analysis binning
  bool warned = false;
  for ( int i = 0; i < 7; ++i ) {
    TH1* h = (TH1*) fchk.Get( names[i] );
    if ( !h ) continue;
    if ( h->GetNbinsX() != expect[i] ) {
      if ( !warned ) {
        printf( "\n  [STALE INPUT] bin counts disagree with the current analysis"
                " binning -- is %s pre-coarsening?\n", infile );
        warned = true;
      }
      printf( "      %-12s produced %d bins, analysis uses %d\n",
              names[i], h->GetNbinsX(), expect[i] );
    }
  }
  fchk.Close();
  if ( warned ) printf( "      -> regenerate the *_final.root inputs with the"
                        " coarsened binning before trusting this file.\n" );
  printf( "wrote %s (5 direct + 2 reversed theta observables)\n", outfile );
}
