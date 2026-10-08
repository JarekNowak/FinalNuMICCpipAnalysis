#pragma once
// Study selections with the multi-pion particle classifier as the pion identification
// (report/planning/multipion/PHASE1_SUMMARY.md). Each is its parent selection with use_new_pid() on; every
// other cut is inherited. The thresholds give the pion efficiency of the parent's pion
// identification on the held-out runs 3 and 5 (scripts/mp_pid_train.py --unmatched --common, the
// model on the 43 inputs every sample carries): for the single-pion
// selection on its own pool (contained, vertex distance <= 4 cm, length > 20 cm), for the
// multi-pion selections that of the split mp_pion_bdt. Not part of any released result.
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi1p.hh"
#include "XSecAnalyzer/Selections/CC1mu2pi.hh"
#include "XSecAnalyzer/Selections/CC1mu3pi.hh"

// The threshold can be set at construction: the factory builds "<name>_tNN" as the selection with
// P(pi) > NN/100 (threshold scans; scripts/mp_pid_eval.py).
class CC1mu1piXpNewPID : public CC1mu1piXp {
 public:
  explicit CC1mu1piXpNewPID( double cut = 0.3087, const std::string& name = "CC1mu1piXpNewPID" )
    : CC1mu1piXp(), cut_( cut ) { this->set_selection_name( name ); }
 protected:
  bool   use_new_pid() const override { return true; }
  double new_pid_cut() const override { return cut_; }
  // keeps the 20 cm length requirement of the single-pion pion identification, which
  // matches the measured phase space p_pi > 0.175 GeV/c
  double new_pid_min_length() const override { return 20.; }
 private:
  double cut_;
};

class CC1mu2piNewPID : public CC1mu2pi {
 public:
  explicit CC1mu2piNewPID( double cut = 0.2276, const std::string& name = "CC1mu2piNewPID" )
    : CC1mu2pi(), cut_( cut ) { this->set_selection_name( name ); }
 protected:
  bool   use_new_pid() const override { return true; }
  double new_pid_cut() const override { return cut_; }
 private:
  double cut_;
};

class CC1mu3piNewPID : public CC1mu3pi {
 public:
  explicit CC1mu3piNewPID( double cut = 0.2276, const std::string& name = "CC1mu3piNewPID" )
    : CC1mu3pi(), cut_( cut ) { this->set_selection_name( name ); }
 protected:
  bool   use_new_pid() const override { return true; }
  double new_pid_cut() const override { return cut_; }
 private:
  double cut_;
};

// The released single-pion selection with the PID diagnostic switched on: every decision and
// output of CC1mu1piXp, plus the per-candidate BDT inputs and outputs of the pion pool and the
// muon candidate (pdiag_*, mudiag_*), for the data/prediction comparison of the PID variables in
// the opened control regions (scripts/cr_pid_diag_check.py).
class CC1mu1piXpPIDDiag : public CC1mu1piXp {
 public:
  CC1mu1piXpPIDDiag() : CC1mu1piXp() { this->set_selection_name( "CC1mu1piXpPIDDiag" ); }
 protected:
  bool store_pid_diag() const override { return true; }
};

// The released single-pion selections without the Bragg-pion >= 0.08 cut, in every sample, written
// to size the fail-open effect of the older-production ntuples at event level
// (scripts/bragg_cut_eval.py). Since the cut was dropped from the released selections (2026-10-02)
// they are identical to their parents.
class CC1mu1piXpNoBragg : public CC1mu1piXp {
 public:
  CC1mu1piXpNoBragg() : CC1mu1piXp() { this->set_selection_name( "CC1mu1piXpNoBragg" ); }
 protected:
  bool apply_bragg_pion_cut() const override { return false; }
};

class CC1mu1pi1pNoBragg : public CC1mu1pi1p {
 public:
  CC1mu1pi1pNoBragg() : CC1mu1pi1p() { this->set_selection_name( "CC1mu1pi1pNoBragg" ); }
 protected:
  bool apply_bragg_pion_cut() const override { return false; }
};

// The current two-pion selection with the event-level track dump switched on (etrk_*, evt_*): every
// decision and output of CC1mu2pi, plus the class probabilities and kinematics of every track of the
// particle classifier's domain, for the assignment and the event classifier of Phase 2
// (scripts/mp_evt_*.py).
class CC1mu2piEvt : public CC1mu2pi {
 public:
  CC1mu2piEvt() : CC1mu2pi() { this->set_selection_name( "CC1mu2piEvt" ); }
 protected:
  bool store_event_tracks() const override { return true; }
};

// The two-pion selection of Phases 2 and 4 (report/planning/multipion/PHASE2_SUMMARY.md, PHASE4_SUMMARY.md): the event
// classifier on the particle classifier with the LLR score and geometry only, signal region at score > 0.92, the four
// sidebands in <name>_evt_region. Its selection is the classifier's, not the pion count of CC1mu2pi.
class CC1mu2piBDT : public CC1mu2pi {
 public:
  CC1mu2piBDT() : CC1mu2pi() { this->set_selection_name( "CC1mu2piBDT" ); }
  bool selection( AnalysisEvent* event ) override { return CC1mu1piXp::selection( event ); }
 protected:
  bool use_event_classifier() const override { return true; }
  std::string mp_pid_model_file() const override { return "mp_pid_rbdt_unmatched_common_llronly.root"; }
  std::string evt_clf_model_file() const override { return "evt_clf_n2_llronly_pid_rbdt.root"; }
};
