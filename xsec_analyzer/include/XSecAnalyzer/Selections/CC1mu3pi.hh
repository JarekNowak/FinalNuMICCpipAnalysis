#pragma once
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"

// nu_mu / nu_mu-bar CC with exactly 3 charged pions (and any number of protons)
// in the final state. A thin multi-pion subselection of CC1mu1piXp: it inherits the
// full selection, signal definition and observables and overrides only the charged-
// pion multiplicity -- signal requires 3 true charged pions, and the reco selection
// additionally requires pion_number == 3. Low statistics: intended first for the
// cut-flow / feasibility study (differential cross sections may not be viable).
class CC1mu3pi : public CC1mu1piXp {
 public:
  CC1mu3pi();
  bool selection( AnalysisEvent* event ) override;
 protected:
  int required_charged_pions() const override { return 3; }
};
