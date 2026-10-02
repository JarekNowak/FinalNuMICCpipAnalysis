#pragma once
// Study selections with the multi-pion particle classifier as the pion identification
// (report/multipion/PHASE1_SUMMARY.md). Each is its parent selection with use_new_pid() on; every
// other cut is inherited. The thresholds give the pion efficiency of the parent's pion
// identification on the held-out runs 3 and 5 (scripts/mp_pid_train.py --unmatched --common, the
// model on the 43 inputs every sample carries): for the single-pion
// selection on its own pool (contained, vertex distance <= 4 cm, length > 20 cm), for the
// multi-pion selections that of the split mp_pion_bdt. Not part of any released result.
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"
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
