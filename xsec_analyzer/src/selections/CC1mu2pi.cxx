#include "XSecAnalyzer/Selections/CC1mu2pi.hh"

CC1mu2pi::CC1mu2pi() : CC1mu1piXp() {
  // Rename so all output branches use the "CC1mu2pi_" prefix.
  this->set_selection_name( "CC1mu2pi" );
}

// Reco selection: the inherited CC1mu1piXp selection PLUS exactly 2 reconstructed
// charged-pion candidates (the base gates topology/containment but not the pion count).
bool CC1mu2pi::selection( AnalysisEvent* Event ) {
  bool base = CC1mu1piXp::selection( Event );
  return base && ( pion_number == required_charged_pions() );
}
