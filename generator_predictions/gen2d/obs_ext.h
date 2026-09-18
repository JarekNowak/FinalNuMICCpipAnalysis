// obs_ext.h -- generator-side definitions of the observables added 2026-09-16:
//   proton-tagged 1D : theta_p (4 bins) and the pion-proton opening angle theta_pipr (3 bins)
//   2D (flattened)   : 4 inclusive and 4 proton-tagged pairs
// THE EDGES MUST MATCH the bin configs written by xsec_analyzer/scripts/make_bin_configs.py from
// xsec_analyzer/scripts/newobs_candidates.py; check_ext_bins.py compares the two.
// 2D analysis bins are flattened Y slice by Y slice (j*nx + i), stored as a TH1D over the bin
// index (width 1), so each reader's usual "/width" normalisation leaves the per-bin cross section
// and the generic combiner's "*width" turns it into the FTE content unchanged.
// Axis semantics follow the TRUE bins of the configs: lowOpen -> values below the first edge go
// to the first bin, topOpen -> values above the last edge go to the last bin, otherwise out.
// Truth frame: generator +z is the neutrino direction.
#pragma once
#include <vector>
#include <string>
#include "TH1D.h"
#include "TVector3.h"
#include "../wtki_gen.h"

namespace ext {
  struct Axis { std::vector<double> e; bool lowOpen; bool topOpen; };
  inline int bin( double v, const Axis& a ) {
    int n = (int)a.e.size() - 1;
    if ( v < a.e[0] ) return a.lowOpen ? 0 : -1;
    for ( int b = 0; b < n; ++b ) if ( v < a.e[b+1] ) return b;
    return a.topOpen ? n - 1 : -1;
  }
  struct Pair { const char* name; Axis x; Axis y; };
  inline int flat( const Pair& p, double x, double y ) {
    int i = bin( x, p.x ), j = bin( y, p.y );
    if ( i < 0 || j < 0 ) return -1;
    return j * ( (int)p.x.e.size() - 1 ) + i;
  }
  inline int nflat( const Pair& p ) { return ( (int)p.x.e.size() - 1 ) * ( (int)p.y.e.size() - 1 ); }

  const double PI = 3.1416;
  // ---- inclusive 2D (X | Y) --------------------------------------------------------------
  inline std::vector<Pair> incl_pairs() {
    std::vector<Pair> v;
    Pair a = { "costhpi_costhmu", { {-1, 0.35, 1}, false, false }, { {-1, 0.65, 0.85, 1}, false, false } }; v.push_back( a );
    Pair b = { "costhmu_pmu",     { {-1, 0.8, 1},  false, false }, { {0.15, 0.55, 3.0}, false, true } };    v.push_back( b );
    Pair c = { "costhpi_thmupi",  { {-1, 0.35, 1}, false, false }, { {0, 0.85, 1.5, 2.6}, false, false } }; v.push_back( c );
    Pair d = { "costhpi_ppi",     { {-1, 0.35, 1}, false, false }, { {0.175, 0.205, 1.0}, false, true } };  v.push_back( d );
    return v;
  }
  // ---- proton-tagged 2D (X | Y) ----------------------------------------------------------
  inline std::vector<Pair> p1_pairs() {
    std::vector<Pair> v;
    Pair a = { "thetap_dpt",   { {0, 0.8, PI}, false, false }, { {0, 0.3, 2.5}, true, true } };      v.push_back( a );
    Pair b = { "thpipr_dpt",   { {0, 1.2, PI}, false, false }, { {0, 0.3, 2.5}, true, true } };      v.push_back( b );
    Pair c = { "thetap_Wpipr", { {0, 0.8, PI}, false, false }, { {1.08, 1.19, 2.9}, true, true } };  v.push_back( c );
    Pair d = { "thpipr_Wpipr", { {0, 1.2, PI}, false, false }, { {1.08, 1.19, 2.9}, true, true } };  v.push_back( d );
    return v;
  }
  // ---- proton-tagged 1D ------------------------------------------------------------------
  inline Axis thetap_axis() { Axis a = { {0, 0.659, 1.068, 1.539, PI}, false, false }; return a; }  // RHC maximin, adopted 2026-09-18
  inline Axis thpipr_axis() { Axis a = { {0, 1.193, 2.010, PI}, false, false }; return a; }

  struct Hists {
    std::vector<TH1D*> all;
    std::vector<Pair> pairs; std::vector<TH1D*> h2;
    TH1D* hthetap = 0; TH1D* hthpipr = 0;
    explicit Hists( bool proton_tagged ) {
      pairs = proton_tagged ? p1_pairs() : incl_pairs();
      for ( size_t k = 0; k < pairs.size(); ++k ) {
        TH1D* h = new TH1D( pairs[k].name, "", nflat( pairs[k] ), 0., (double)nflat( pairs[k] ) );
        h2.push_back( h ); all.push_back( h );
      }
      if ( proton_tagged ) {
        Axis t = thetap_axis(), o = thpipr_axis();
        hthetap = new TH1D( "thetap", "", (int)t.e.size() - 1, t.e.data() ); all.push_back( hthetap );
        hthpipr = new TH1D( "thpipr", "", (int)o.e.size() - 1, o.e.data() ); all.push_back( hthpipr );
      }
    }
    void fill1d( TH1D* h, const Axis& a, double v, double w ) { int b = bin( v, a ); if ( b >= 0 ) h->Fill( h->GetBinCenter( b + 1 ), w ); }
    void fill2d( size_t k, double x, double y, double w ) { int i = flat( pairs[k], x, y ); if ( i >= 0 ) h2[k]->Fill( i + 0.5, w ); }
    // inclusive signal event (after the 0.175 GeV/c pion phase space)
    void fill_incl( const TVector3& mu, const TVector3& pi, double w = 1. ) {
      double cmu = mu.CosTheta(), cpi = pi.CosTheta(), th = mu.Angle( pi );
      fill2d( 0, cpi, cmu, w ); fill2d( 1, cmu, mu.Mag(), w ); fill2d( 2, cpi, th, w ); fill2d( 3, cpi, pi.Mag(), w );
    }
    // proton-tagged signal event
    void fill_1p( const TVector3& mu, const TVector3& pi, const TVector3& pr, double w = 1. ) {
      wtki::Obs o = wtki::compute( mu, pi, pr );
      double thp = pr.Theta(), tpp = pi.Angle( pr );
      fill1d( hthetap, thetap_axis(), thp, w ); fill1d( hthpipr, thpipr_axis(), tpp, w );
      fill2d( 0, thp, o.dpt, w ); fill2d( 1, tpp, o.dpt, w ); fill2d( 2, thp, o.Wpipr, w ); fill2d( 3, tpp, o.Wpipr, w );
    }
  };
}
