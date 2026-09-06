// sb_guard.h -- blinding guard for the sideband macros (2026-09-06).
// A beam-on file may be read by a control-region macro only if it is the
// signal-region-stripped skim written by ProcessNTuples with XSEC_CR_SKIM=1 (TNamed
// marker XSEC_CR_SKIM), or if it is fake data (is_mc true: a throw of the neutrino MC).
// A real beam-on file WITHOUT the marker is refused unless XSEC_UNBLIND=1, mirroring
// the guard in SystematicsCalculator. Call sb_guard_data(path) on every entry of a
// macro's beam-on list before opening it.
#pragma once
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include <stdexcept>
#include <string>
inline void sb_guard_data(const std::string& path){
  TFile f(path.c_str(),"read");
  if(f.IsZombie()) throw std::runtime_error("sb_guard: cannot open "+path);
  if(f.Get("XSEC_CR_SKIM")) return;                 // the control-region skim: allowed
  TTree* t=(TTree*)f.Get("stv_tree");
  if(!t) throw std::runtime_error("sb_guard: no stv_tree in "+path);
  bool is_mc=false; t->SetBranchStatus("*",0); t->SetBranchStatus("is_mc",1);
  t->SetBranchAddress("is_mc",&is_mc); if(t->GetEntries()>0) t->GetEntry(0);
  if(is_mc) return;                                  // fake data (thrown from the MC): allowed
  const char* u=std::getenv("XSEC_UNBLIND");
  if(u && std::string(u)=="1") return;
  throw std::runtime_error("BLINDING GUARD (sb_guard): "+path+" is a beam-on file without MC truth"
    " and without the XSEC_CR_SKIM marker, i.e. an unstripped real beam-on file. The sideband"
    " macros read the control-region skim only (ProcessNTuples with XSEC_CR_SKIM=1).");
}
