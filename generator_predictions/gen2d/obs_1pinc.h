// obs_1pinc.h -- generator-side definitions of the muon/pion observables measured on the PROTON-TAGGED
// (CC1mu1pi1p) sample, at the proton-tagged binnings (2026-09-27). Filled by the *_1p_inc readers, which
// are the 1p readers (analyze_gst_1p.C, gibuu_1p.C, neut_1p.C, nuwro_1p.cc) with these fills added, so the
// signal selection and normalisation are exactly those of the released 1p_ext histograms.
// THE EDGES MUST MATCH the true bins of
//   xsec_analyzer/configs/ccpi1p_pmu_bin_config.txt       p_mu        6 bins, last open
//   xsec_analyzer/configs/ccpi1p_ppi_bin_config_2bin.txt  p_pi        2 bins, last open
//   xsec_analyzer/configs/ccpi1p_costhmu_bin_config.txt   cos th_mu   8 bins, closed at 1
//   xsec_analyzer/configs/ccpi1p_costhpi_bin_config.txt   cos th_pi   4 bins, last open (>= 0.75)
//   xsec_analyzer/configs/ccpi1p_thmupi_bin_config.txt    th_mu_pi    5 bins, last open
// (check_1pinc_bins.py compares them). The outer edge of an open bin is the slice-config plotting edge; it
// cancels in the FTE (the readers divide by the width, the combiner multiplies by it again).
// Observables as in the inclusive readers: |p_mu|, |p_pi|, cos theta about +z (the generator +z is the
// neutrino direction) and the 3D mu-pi opening angle. Axis semantics: those of ext::bin (obs_ext.h).
#pragma once
#include "obs_ext.h"

namespace p1inc {
  inline ext::Axis pmu_axis()     { ext::Axis a = { {0.15, 0.30, 0.425, 0.575, 0.825, 1.25, 3.0}, false, true }; return a; }
  inline ext::Axis ppi_axis()     { ext::Axis a = { {0.175, 0.205, 1.0}, false, true }; return a; }
  inline ext::Axis costhmu_axis() { ext::Axis a = { {-1, 0.4, 0.625, 0.75, 0.825, 0.9, 0.95, 0.975, 1}, false, false }; return a; }
  inline ext::Axis costhpi_axis() { ext::Axis a = { {-1, -0.1, 0.35, 0.75, 1}, false, true }; return a; }
  inline ext::Axis thmupi_axis()  { ext::Axis a = { {0, 0.6, 0.85, 1.3, 1.85, 2.6}, false, true }; return a; }

  struct Hists {
    std::vector<TH1D*> all;
    TH1D *hpmu, *hppi, *hcmu, *hcpi, *hth;
    TH1D* book( const char* n, const ext::Axis& a ) { TH1D* h = new TH1D( n, "", (int)a.e.size() - 1, a.e.data() ); all.push_back( h ); return h; }
    Hists() {
      hpmu = book( "pmu1p", pmu_axis() ); hppi = book( "ppi1p", ppi_axis() ); hcmu = book( "costhmu1p", costhmu_axis() );
      hcpi = book( "costhpi1p", costhpi_axis() ); hth = book( "thmupi1p", thmupi_axis() );
    }
    void fill1d( TH1D* h, const ext::Axis& a, double v, double w ) { int b = ext::bin( v, a ); if ( b >= 0 ) h->Fill( h->GetBinCenter( b + 1 ), w ); }
    // proton-tagged signal event (after the 1p reader's full signal selection)
    void fill( const TVector3& mu, const TVector3& pi, double w = 1. ) {
      fill1d( hpmu, pmu_axis(), mu.Mag(), w ); fill1d( hppi, ppi_axis(), pi.Mag(), w );
      fill1d( hcmu, costhmu_axis(), mu.CosTheta(), w ); fill1d( hcpi, costhpi_axis(), pi.CosTheta(), w );
      fill1d( hth, thmupi_axis(), mu.Angle( pi ), w );
    }
  };
}
