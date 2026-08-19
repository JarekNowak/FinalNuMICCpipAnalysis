// ROOT includes
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu1pi1p.hh"
#include "XSecAnalyzer/STVTools.hh"
#include "XSecAnalyzer/NuMIBeamFrame.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/AnalysisEvent.hh"

namespace {
  // True neutrino direction in detector coordinates (see NuMIBeamFrame.hh). Truth-level
  // angles are referred to it so they match the generator predictions exactly.
  TVector3 true_nu_dir( const AnalysisEvent* ev ) {
    TVector3 d( ev->mc_nu_px_, ev->mc_nu_py_, ev->mc_nu_pz_ );
    if ( d.Mag() > 0. ) return d.Unit();
    return NuMIBeam::axis();
  }

  // Reco-level neutrino direction: the beta convention (target -> reco vertex).
  TVector3 reco_nu_dir( const AnalysisEvent* ev ) {
    return NuMIBeam::nu_dir_from_vertex( ev->nu_vx_, ev->nu_vy_, ev->nu_vz_ );
  }

  // Proton identification LLR-PID working point. Tightened from the framework
  // DEFAULT_PROTON_PID_CUT (0.2) to 0.05 after a cut scan on FHC: this maximises
  // S/sqrt(B) (+2%) and lifts purity (+1.7 pts) at ~0.3% efficiency cost, with the
  // W_pipr resolution unchanged (~9%). Lower LLR = more proton-like.
  constexpr double PROTON_LLR_CUT = 0.05;

  // safe sqrt: returns 0 for a (small) negative argument from reco fluctuations
  inline double ssqrt( double x ) { return x > 0. ? std::sqrt(x) : 0.; }
  // proton momentum (GeV/c) from proton-hypothesis kinetic energy (GeV)
  inline double proton_mom_from_KE( double KE ) {
    return ssqrt( KE*KE + 2.*PROTON_MASS*KE );
  }
}

CC1mu1pi1p::CC1mu1pi1p() : CC1mu1piXp() {
  // Rename so all output branches use the "CC1mu1pi1p_" prefix (the parent
  // constructor set it to "CC1mu1piXp").
  this->set_selection_name( "CC1mu1pi1p" );
}

// ---------------------------------------------------------------------------
// Signal: parent CC1mu1piXp signal PLUS at least one true proton above the
// tracking threshold (0.3 GeV/c).
bool CC1mu1pi1p::define_signal( AnalysisEvent* Event ) {
  bool base = CC1mu1piXp::define_signal( Event );

  double lead_p = 0.;
  size_t n = Event->mc_nu_daughter_pdg_->size();
  for ( size_t i = 0; i < n; ++i ) {
    if ( Event->mc_nu_daughter_pdg_->at(i) != PROTON ) continue;
    TVector3 p( Event->mc_nu_daughter_px_->at(i),
      Event->mc_nu_daughter_py_->at(i), Event->mc_nu_daughter_pz_->at(i) );
    if ( p.Mag() > lead_p ) lead_p = p.Mag();
  }
  lead_true_proton_mom_ = lead_p;

  return base && ( lead_p > 0.3 ); // GeV/c
}

// ---------------------------------------------------------------------------
// Selection: parent CC1mu1piXp selection PLUS an identified proton candidate
// (leading, most proton-like track that is contained and above threshold).
bool CC1mu1pi1p::selection( AnalysisEvent* Event ) {
  bool base = CC1mu1piXp::selection( Event ); // sets CandidateMuonIndex/PionIndex

  CandidateProtonIdx_ = BOGUS_INDEX;
  reco_n_proton_ = 0;
  double best_llr = 1.e9; // lowest LLR PID score = most proton-like

  size_t n = Event->track_llr_pid_score_->size();
  for ( size_t i = 0; i < n; ++i ) {
    if ( (int)i == CandidateMuonIndex || (int)i == CandidatePionIndex ) continue;
    // track-like (reject showers)
    if ( Event->pfp_track_score_->at(i) < TRACK_SCORE_CUT ) continue;
    // proton-like PID
    double llr = Event->track_llr_pid_score_->at(i);
    if ( llr >= PROTON_LLR_CUT ) continue;   // tightened proton PID (was DEFAULT_PROTON_PID_CUT=0.2)
    // above the proton tracking threshold
    double pmom = proton_mom_from_KE( Event->track_kinetic_energy_p_->at(i) );
    if ( pmom < 0.3 ) continue;
    // contained (track end inside the reco FV)
    if ( !point_inside_FV( this->reco_FV(), Event->track_endx_->at(i),
      Event->track_endy_->at(i), Event->track_endz_->at(i) ) ) continue;

    ++reco_n_proton_;
    if ( llr < best_llr ) { best_llr = llr; CandidateProtonIdx_ = (int)i; }
  }

  return base && ( CandidateProtonIdx_ != BOGUS_INDEX );
}

// ---------------------------------------------------------------------------
// Shared W / TKI kinematics from muon, pion, proton 3-momenta and energies.
void CC1mu1pi1p::compute_had_observables( const TVector3& mu, double Emu,
  const TVector3& pi, double Epi, const TVector3& pr, double Epr,
  const TVector3& nu_dir,
  double& W_pipr, double& W_had, double& dalphaT, double& dphiT,
  double& dpT, double& pn ) {

  // W_pipr : invariant mass of the (pion + proton) system. Directly measured,
  // no neutrino-direction assumption; the Delta/N* resonance mass.
  TLorentzVector p4pi( pi, Epi );
  TLorentzVector p4pr( pr, Epr );
  W_pipr = ( p4pi + p4pr ).M();

  // TKI + calorimetric W_had from the muon vs. the visible (pion+proton)
  // hadronic system, via the framework STV tool.
  //
  // STVTools takes the DETECTOR z axis as longitudinal throughout - correct for BNB,
  // where the beam runs along z, but wrong for NuMI, whose neutrinos arrive 28 degrees
  // away from it. Rotating the momenta into a frame whose +z is the neutrino direction
  // makes that assumption true again, so STVTools itself is left untouched and the BNB
  // selections that share it are unaffected. Without this the axis error swamps the
  // observable: mean delta_pT comes out at 0.690 GeV/c about z against 0.267 GeV/c
  // about the beam. See NuMIBeamFrame.hh.
  TVector3 mu_b  = NuMIBeam::to_beam_frame( mu, nu_dir );
  TVector3 had_b = NuMIBeam::to_beam_frame( pi + pr, nu_dir );
  double Ehad = Epi + Epr;
  STVTools stv;
  stv.CalculateSTVs( mu_b, had_b, Emu, Ehad, CalcType );
  dpT     = stv.ReturnPt();
  pn      = stv.ReturnPn();
  dalphaT = stv.ReturnDeltaAlphaT();
  dphiT   = stv.ReturnDeltaPhiT();

  double Ecal = stv.ReturnECal();
  double Q2   = stv.ReturnQ2();
  W_had = ssqrt( PROTON_MASS*PROTON_MASS + 2.*PROTON_MASS*( Ecal - Emu ) - Q2 );
}

// ---------------------------------------------------------------------------
void CC1mu1pi1p::compute_reco_observables( AnalysisEvent* Event ) {
  // Parent sets CandidateMuonIndex/PionIndex, candidate_muon/pion_mom_reco.
  CC1mu1piXp::compute_reco_observables( Event );

  if ( CandidateProtonIdx_ == BOGUS_INDEX || CandidateMuonIndex < 0
    || CandidatePionIndex < 0 ) return;

  // muon
  TVector3 mu( Event->track_dirx_->at(CandidateMuonIndex),
    Event->track_diry_->at(CandidateMuonIndex),
    Event->track_dirz_->at(CandidateMuonIndex) );
  mu.SetMag( candidate_muon_mom_reco );
  double Emu = ssqrt( candidate_muon_mom_reco*candidate_muon_mom_reco
    + MUON_MASS*MUON_MASS );

  // pion
  TVector3 pi( Event->track_dirx_->at(CandidatePionIndex),
    Event->track_diry_->at(CandidatePionIndex),
    Event->track_dirz_->at(CandidatePionIndex) );
  pi.SetMag( candidate_pion_mom_reco );
  double Epi = ssqrt( candidate_pion_mom_reco*candidate_pion_mom_reco
    + PI_PLUS_MASS*PI_PLUS_MASS );

  // proton
  double KEp = Event->track_kinetic_energy_p_->at( CandidateProtonIdx_ );
  double pmom = proton_mom_from_KE( KEp );
  TVector3 pr( Event->track_dirx_->at(CandidateProtonIdx_),
    Event->track_diry_->at(CandidateProtonIdx_),
    Event->track_dirz_->at(CandidateProtonIdx_) );
  pr.SetMag( pmom );
  double Epr = KEp + PROTON_MASS;
  reco_proton_mom_ = pmom;
  reco_proton_costh_ = pr.Unit().Dot( reco_nu_dir(Event) );
  reco_proton_llr_ = Event->track_llr_pid_score_->at( CandidateProtonIdx_ );

  // Reconstructed: the beta convention of the NuMI nue CC1pi internal note - the
  // neutrino direction approximated by the NuMI target to the reconstructed vertex.
  // The true per-event direction is not recoverable in data; what beta leaves behind
  // is the physical divergence of the beam, a resolution effect the response absorbs.
  compute_had_observables( mu, Emu, pi, Epi, pr, Epr, reco_nu_dir(Event),
    reco_W_pipr_, reco_W_had_,
    reco_deltaAlphaT_, reco_deltaPhiT_, reco_deltaPt_, reco_pn_ );
}

// ---------------------------------------------------------------------------
void CC1mu1pi1p::compute_true_observables( AnalysisEvent* Event ) {
  // Parent sets TrueCandidateMuonP / TrueCandidatePionP (in define_signal).
  CC1mu1piXp::compute_true_observables( Event );

  // leading true proton
  double lead_p = 0.; int idx = BOGUS_INDEX;
  size_t n = Event->mc_nu_daughter_pdg_->size();
  for ( size_t i = 0; i < n; ++i ) {
    if ( Event->mc_nu_daughter_pdg_->at(i) != PROTON ) continue;
    TVector3 p( Event->mc_nu_daughter_px_->at(i),
      Event->mc_nu_daughter_py_->at(i), Event->mc_nu_daughter_pz_->at(i) );
    if ( p.Mag() > lead_p ) { lead_p = p.Mag(); idx = (int)i; }
  }
  true_lead_proton_idx_ = idx;
  lead_true_proton_mom_ = lead_p;

  if ( idx == BOGUS_INDEX || TrueCandidateMuonP.Mag() <= 0.
    || TrueCandidatePionP.Mag() <= 0. ) return;

  TVector3 pr( Event->mc_nu_daughter_px_->at(idx),
    Event->mc_nu_daughter_py_->at(idx), Event->mc_nu_daughter_pz_->at(idx) );
  double Epr = ssqrt( pr.Mag2() + PROTON_MASS*PROTON_MASS );
  double Emu = ssqrt( TrueCandidateMuonP.Mag2() + MUON_MASS*MUON_MASS );
  double Epi = ssqrt( TrueCandidatePionP.Mag2() + PI_PLUS_MASS*PI_PLUS_MASS );
  true_proton_mom_ = pr.Mag();
  true_proton_costh_ = pr.Unit().Dot( true_nu_dir(Event) );

  // Truth: use the exact per-event neutrino direction, which is how every generator
  // prediction defines these variables (the neutrino is +z by construction there).
  // Fall back to the fixed axis only if the ntuple lacks the truth momentum branches.
  compute_had_observables( TrueCandidateMuonP, Emu, TrueCandidatePionP, Epi,
    pr, Epr, true_nu_dir(Event), true_W_pipr_, true_W_had_, true_deltaAlphaT_,
    true_deltaPhiT_, true_deltaPt_, true_pn_ );
}

// ---------------------------------------------------------------------------
void CC1mu1pi1p::define_output_branches() {
  // Parent registers the shared branches (signal flags, muon/pion observables)
  // under the (now CC1mu1pi1p_) prefix.
  CC1mu1piXp::define_output_branches();

  set_branch( &reco_W_pipr_,      "W_pipr_reco" );
  set_branch( &reco_W_had_,       "W_had_reco" );
  set_branch( &reco_deltaAlphaT_, "deltaAlphaT_reco" );
  set_branch( &reco_deltaPhiT_,   "deltaPhiT_reco" );
  set_branch( &reco_deltaPt_,     "deltaPt_reco" );
  set_branch( &reco_pn_,          "pn_reco" );
  set_branch( &reco_proton_mom_,  "proton_mom_reco" );
  set_branch( &reco_proton_costh_,"proton_costh_reco" );
  set_branch( &reco_n_proton_,    "n_proton_reco" );
  set_branch( &reco_proton_llr_,  "proton_llr_reco" );

  set_branch( &true_W_pipr_,      "W_pipr_true" );
  set_branch( &true_W_had_,       "W_had_true" );
  set_branch( &true_deltaAlphaT_, "deltaAlphaT_true" );
  set_branch( &true_deltaPhiT_,   "deltaPhiT_true" );
  set_branch( &true_deltaPt_,     "deltaPt_true" );
  set_branch( &true_pn_,          "pn_true" );
  set_branch( &true_proton_mom_,  "proton_mom_true" );
  set_branch( &true_proton_costh_,"proton_costh_true" );
  set_branch( &lead_true_proton_mom_, "lead_proton_mom_true" );
}

// ---------------------------------------------------------------------------
void CC1mu1pi1p::reset() {
  CC1mu1piXp::reset();

  CandidateProtonIdx_ = BOGUS_INDEX;
  reco_n_proton_ = 0;

  reco_W_pipr_ = reco_W_had_ = BOGUS;
  reco_deltaAlphaT_ = reco_deltaPhiT_ = reco_deltaPt_ = reco_pn_ = BOGUS;
  reco_proton_mom_ = reco_proton_costh_ = BOGUS;
  reco_proton_llr_ = BOGUS;

  true_W_pipr_ = true_W_had_ = BOGUS;
  true_deltaAlphaT_ = true_deltaPhiT_ = true_deltaPt_ = true_pn_ = BOGUS;
  true_proton_mom_ = true_proton_costh_ = BOGUS;
  true_lead_proton_idx_ = BOGUS_INDEX;
  lead_true_proton_mom_ = BOGUS;
}
