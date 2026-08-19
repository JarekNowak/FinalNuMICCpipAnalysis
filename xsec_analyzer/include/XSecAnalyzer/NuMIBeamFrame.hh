#pragma once

// NuMI beam geometry in MicroBooNE detector coordinates.
//
// The NuMI neutrinos do NOT travel along the detector z axis: the mean direction
// sits 28.6 degrees away from it. Every transverse-kinematic-imbalance variable and
// every polar angle is defined with respect to the NEUTRINO direction, so anything
// built on the detector z axis is measuring about the wrong axis. STVTools takes z
// as longitudinal throughout (correct for BNB, where the beam IS along z), so NuMI
// selections must rotate their momenta into the beam frame before calling it.
//
// Two axes must not be confused:
//   * the NuMI BEAMLINE axis, which sits 9.1 degrees from the neutrinos seen here
//     (MicroBooNE is an off-axis detector), and
//   * the TARGET -> DETECTOR direction, which is what the neutrinos actually follow.
// The second is the correct one. Measured against the true per-event neutrino
// direction over 307k in-detector events, the RMS opening angle is
//   detector z .................. 532 mrad   <- the uncorrected convention
//   NuMI beamline axis .......... 242 mrad
//   target -> detector .......... 178 mrad   <- used here
//   per-event target -> vertex .. 178 mrad   (no gain: see below)
//   flux-weighted mean .......... 175 mrad   (the best any fixed axis can do)
//
// The residual ~178 mrad is the physical divergence of the beam - neutrinos reaching
// the detector come from decays spread along the pipe - and is NOT recoverable from
// the interaction vertex: a 684 m baseline across a 2.5 m detector subtends only
// ~4 mrad, which is why the per-event vertex direction is no better than a constant.
// Truth-level observables therefore use the exact per-event neutrino direction (which
// matches the generator predictions, where the neutrino defines +z by construction)
// and reco-level observables use the fixed axis below; the difference between them is
// a resolution effect and is absorbed by the response matrix like any other.

#include <TRotation.h>
#include <TVector3.h>

namespace NuMIBeam {

  // Detector -> beam rotation, as used by the NuMI beamline geometry weights
  // (src/app/NuMI/addBeamlineGeometryWeightsToMap.cpp).
  inline const TRotation& det_to_beam() {
    static const TRotation R = []{
      TRotation r;
      r.RotateAxes(
        TVector3( 0.92103853804025681562,   0.022713504803924120662, 0.38880857519374290021  ),
        TVector3( 4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359 ),
        TVector3(-0.38947144863934973769,   0.053832413938664107345, 0.91946400794392302291  ) );
      return r;
    }();
    return R;
  }

  // Position of the MicroBooNE detector in beam coordinates (cm).
  inline TVector3 detector_in_beam_coords() { return TVector3( 5502., 7259., 67270. ); }

  // Mean neutrino direction in DETECTOR coordinates: (0.46237, 0.04885, 0.88534),
  // i.e. 27.7 degrees from the detector z axis.
  inline const TVector3& axis() {
    static const TVector3 a =
      ( det_to_beam().Inverse() * detector_in_beam_coords() ).Unit();
    return a;
  }

  // NuMI target position in DETECTOR coordinates (cm), ~684 m upstream.
  inline const TVector3& target() {
    static const TVector3 t =
      -( det_to_beam().Inverse() * detector_in_beam_coords() );
    return t;
  }

  // The angle beta of the NuMI nue CC1pi internal note (v3.2 Sec. 5.3) and of
  // K. Mistry's thesis: approximate the neutrino direction by the vector from the
  // NuMI target to the reconstructed interaction vertex. This is the reco-level
  // convention; at truth level the exact neutrino direction is used instead, exactly
  // as that note does.
  //
  // It differs from the fixed axis() above by a median 2.4 mrad (max 7.3 mrad), which
  // moves 0.42% of muons across a cos(theta_mu) analysis bin edge - so the two are
  // interchangeable in practice. beta is used because it is the established NuMI
  // convention, not because it measurably improves anything: against the true
  // per-event direction it gives a mean opening angle of 1.95 degrees where the fixed
  // axis gives 2.00.
  inline TVector3 nu_dir_from_vertex( double vx, double vy, double vz ) {
    TVector3 d = TVector3( vx, vy, vz ) - target();
    if ( d.Mag() <= 0. ) return axis();
    return d.Unit();
  }

  // Rotate a detector-frame vector into a frame whose +z is `nu_dir`. Any such
  // rotation will do: delta_pT and p_n are magnitudes, and delta_alpha_T / delta_phi_T
  // are angles measured within the transverse plane relative to the muon's own
  // transverse momentum, so all four are invariant under rotations about the beam.
  inline TVector3 to_beam_frame( const TVector3& v, const TVector3& nu_dir ) {
    TVector3 out = v;
    out.RotateZ( -nu_dir.Phi() );
    out.RotateY( -nu_dir.Theta() );
    return out;
  }

  // Same, about the fixed axis (the only option for reconstructed quantities).
  inline TVector3 to_beam_frame( const TVector3& v ) { return to_beam_frame( v, axis() ); }

}
