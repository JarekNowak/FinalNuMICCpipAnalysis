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
  // Multi-pion tuning matching the dedicated CC-Npi study (I. Pophale MSc):
  // wider pion vertex window, looser pion ID (LLR only), tolerate one uncontained
  // pion, and drop the 1-pion opening-angle / shower cuts. See CC1mu1piXp.hh.
  double pion_vtx_distance_cut()  const override { return 9.5; }
  bool   loose_pion_id()          const override { return true; }
  int    max_uncontained_pions()  const override { return 1; }
  bool   apply_opening_angle_cut() const override { return false; }
  bool   apply_shower_veto()       const override { return false; }
};
