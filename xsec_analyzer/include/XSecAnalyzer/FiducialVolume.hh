#pragma once
#include "TVector3.h"
#include "Constants.hh"
#include <iostream>
// Definition of the fiducial volume for the analysis, together with a utility
// function template to check if a given point lies within it or not

//DB
//Fiducial volume definitions are actually different between Stephen's and Samantha's analyses. See the definition for Stephen's above
//(now being used as the default). Samantha's defintion is here: https://github.com/ssfehlberg/CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples/helper_funcs.h#L18
//Discusions with Afro stated that as long as a consistent definition of FV being used, it should be ok:
//https://microboone.slack.com/archives/D04A8CUB1EW/p1691423102739129 . Thus using Stephen's defintion as the default

struct FiducialVolume{
  double X_Min, X_Max, Y_Min, Y_Max, Z_Min, Z_Max = BOGUS;
};

// Use a template here so that this function can take float or double values as
// input
template <typename Number> bool point_inside_FV( FiducialVolume FV, Number x, Number y, Number z ) {
  bool x_inside_FV = ( FV.X_Min < x ) && ( x < FV.X_Max );
  bool y_inside_FV = ( FV.Y_Min < y ) && ( y < FV.Y_Max );
  bool z_inside_FV = ( FV.Z_Min < z ) && ( z < FV.Z_Max );
  
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}

inline bool point_inside_FV( FiducialVolume FV, const TVector3& pos ) {
  return point_inside_FV( FV, pos.X(), pos.Y(), pos.Z() );
}

// Returns the number of Ar nuclei inside the fiducial volume.
// Optionally subtracts an excluded z dead-region slab (spanning the full x-y
// face) when dead_z_max > dead_z_min; defaults to no subtraction so existing
// callers are unaffected.
inline double num_Ar_targets_in_FV(FiducialVolume FV,
  double dead_z_min = 0., double dead_z_max = 0.) {
  double volume = ( FV.X_Max - FV.X_Min ) * ( FV.Y_Max - FV.Y_Min )
    * ( FV.Z_Max - FV.Z_Min ); // cm^3
  if ( dead_z_max > dead_z_min ) {
    volume -= ( FV.X_Max - FV.X_Min ) * ( FV.Y_Max - FV.Y_Min )
      * ( dead_z_max - dead_z_min );
  }
  constexpr double m_mol_Ar = 39.948; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_LAr = 1.3836; // g/cm^3

  double num_Ar = volume * mass_density_LAr * N_Avogadro / m_mol_Ar;
  return num_Ar;
}

// Returns the total BNB muon neutrino flux (numu / cm^2) in the fiducial
// volume as a function of a given beam exposure (measured in
// protons-on-target)
// NOTE: This is currently approximated using the flux in the *active volume*.
// TODO: Revisit this approximation
inline double integrated_numu_flux_in_FV( double pot ) {
  // Obtained using the histogram hEnumu_cv (the central-value numu flux
  // in the MicroBooNE active volume as a function of neutrino energy)
  // stored in /pnfs/uboone/persistent/uboonebeam/bnb_gsimple
  // /bnb_gsimple_fluxes_01.09.2019_463_hist/. The ROOT commmands executed
  // were
  // root [3] hEnumu_cv->Scale( 1/(4997.*5e8)/(256.35*233.) )
  // root [4] hEnumu_cv->Integral()
  // (double) 7.3762291e-10
  // See the README file in that same folder for details.
  double flux;
  if (useNuMI) {
    // NuMI FHC numu + numubar integrated flux (E > 60 MeV), per POT per cm^2.
    //
    // Authoritative value from the flux-file author (flux_newg4.root, the new
    // GEANT4 NuMI flux): numu = 4.43515e-10, numubar = 2.37644e-10 -> combined
    // 6.81159e-10 /POT/cm^2, for E > 60 MeV (the 60 MeV cut removes an
    // unphysical ~25-30 MeV artifact spike).
    //
    // The OLD value (2.372512e-9, from the superseded uboone_numi_flux_
    // histograms.root) was ~3.48x too large, which made every extracted cross
    // section ~3.5x too SMALL. Confirmed four independent ways: standalone
    // NuWro (2.7x low) and GENIE (3.2x low) on the old flux; the author's
    // numbers (3.48x); and a GENIE closure -- the framework GENIE tune scaled
    // by 3.48x matches a standalone GENIE on the corrected flux to 3%.
    // See generator_predictions/ for the full cross-check.
    //
    // NOTE: still the active-volume-averaged flux used as an approximation for
    // the fiducial-volume flux (a documented ~few-% residual; see the header
    // comment above). num_Ar counts targets in the FV.
    double numu_per_cm2_per_POT = 6.81159e-10;
    flux = pot * numu_per_cm2_per_POT;
  }
  else {
    constexpr double numu_per_cm2_per_POT_in_AV = 7.3762291e-10;
    flux = pot * numu_per_cm2_per_POT_in_AV; // numu / cm^2
  }
  return flux;
}
