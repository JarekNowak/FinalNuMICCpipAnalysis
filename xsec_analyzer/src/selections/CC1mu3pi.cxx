#include "XSecAnalyzer/Selections/CC1mu3pi.hh"

CC1mu3pi::CC1mu3pi() : CC1mu1piXp() {
  // Rename so all output branches use the "CC1mu3pi_" prefix.
  this->set_selection_name( "CC1mu3pi" );
}

// Reco selection: the inherited CC1mu1piXp selection PLUS exactly 3 reconstructed
// charged-pion candidates (the base gates topology/containment but not the pion count).
bool CC1mu3pi::selection( AnalysisEvent* Event ) {
  bool base = CC1mu1piXp::selection( Event );
  return base && ( pion_number == required_charged_pions() );
}
