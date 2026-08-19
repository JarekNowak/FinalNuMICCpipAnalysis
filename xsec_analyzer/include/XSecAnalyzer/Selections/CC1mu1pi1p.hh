#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"

// Proton-tagged subsample of the CC1mu1piXp (nu_mu CC 1mu 1pi^{+-} X p) selection:
// requires (in addition to the parent signal/selection) at least one identified
// proton in the final state, so the hadronic system can be reconstructed. Adds:
//   - W_pipr : invariant mass of the (charged pion + leading proton) system
//              (directly reconstructed, beam-direction independent; the Delta/N*
//              resonance mass)
//   - W_had  : calorimetric hadronic invariant mass from STVTools (E_cal, Q^2)
//   - transverse kinematic imbalance (delta alpha_T, delta phi_T, delta p_T, p_n)
//     of the muon vs. the visible (pion+proton) hadronic system.
// It INHERITS from CC1mu1piXp and reuses the muon/pion candidate finding, cuts and
// signal definition; the overrides only add the proton requirement and the new
// observables.
class CC1mu1pi1p : public CC1mu1piXp {

public:

  CC1mu1pi1p();

  bool define_signal( AnalysisEvent* event ) override;
  bool selection( AnalysisEvent* event ) override;
  void compute_reco_observables( AnalysisEvent* event ) override;
  void compute_true_observables( AnalysisEvent* event ) override;
  void define_output_branches() override;
  void reset() override;

private:

  // Fill W_pipr, W_had and the TKI variables from the muon, pion and proton
  // 3-momenta and energies (used identically for reco and truth).
  // `nu_dir` is the neutrino direction in detector coordinates, about which the
  // TKI variables are defined: the true per-event direction on the truth side, the
  // fixed NuMIBeam::axis() on the reco side. NOT the detector z axis - see
  // NuMIBeamFrame.hh.
  void compute_had_observables( const TVector3& mu, double Emu,
    const TVector3& pi, double Epi, const TVector3& pr, double Epr,
    const TVector3& nu_dir,
    double& W_pipr, double& W_had, double& dalphaT, double& dphiT,
    double& dpT, double& pn );

  // index of the identified reco proton candidate (leading, most proton-like)
  int CandidateProtonIdx_ = -1;

  // ---- reco proton-tagged observables ----
  double reco_W_pipr_;
  double reco_W_had_;
  double reco_deltaAlphaT_;
  double reco_deltaPhiT_;
  double reco_deltaPt_;
  double reco_pn_;
  double reco_proton_mom_;
  double reco_proton_costh_;
  int    reco_n_proton_;
  // diagnostic: LLR PID score of the identified proton candidate (lower = more
  // proton-like); stored so a tighter proton-PID working point can be emulated
  // offline without reprocessing.
  double reco_proton_llr_;

  // ---- true proton-tagged observables ----
  double true_W_pipr_;
  double true_W_had_;
  double true_deltaAlphaT_;
  double true_deltaPhiT_;
  double true_deltaPt_;
  double true_pn_;
  double true_proton_mom_;
  double true_proton_costh_;
  int    true_lead_proton_idx_;
  double lead_true_proton_mom_;
};
