// throw_guard.h -- one throw file per ensemble member, created once and never replaced (2026-09-26).
// Arrays of different observables share a member's throw (same seed, identical content). Writing a
// temporary and renaming it over the final name (2026-09-20) still failed on NFS: the rename unlinks the
// old file while a job on another node is reading it, and that reader gets short reads (FHC p_mu members
// 8 and 9, 2026-09-26). The final name is therefore created only by link(2) from a closed temporary, which
// fails when the name exists, and an existing throw is reused when its seed and POT scale match and it is
// newer than every input. Only a stale throw is replaced.
#pragma once
#include "TFile.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <vector>

inline Long_t throw_mtime(const TString& f){
  Long_t id, fl, mt; Long64_t sz;
  return gSystem->GetPathInfo(f, &id, &sz, &fl, &mt) == 0 ? mt : -1;
}

// a complete throw with this seed and POT scale, newer than every input
inline bool throw_reusable(const TString& fin, const std::vector<TString>& inputs, int seed, double potscale){
  Long_t t = throw_mtime(fin); if (t < 0) return false;
  for (auto& i : inputs) if (throw_mtime(i) >= t) return false;
  TFile f(fin, "read"); if (f.IsZombie()) return false;
  auto* s = dynamic_cast<TParameter<int>*>(f.Get("throw_seed"));
  auto* p = dynamic_cast<TParameter<double>*>(f.Get("throw_potscale"));
  return s && p && s->GetVal() == seed && std::fabs(p->GetVal() / potscale - 1.) < 1e-9;
}

// the throw's identity, written into the output file before it is closed
inline void throw_stamp(int seed, double potscale){
  TParameter<int>("throw_seed", seed).Write("throw_seed", TObject::kOverwrite);
  TParameter<double>("throw_potscale", potscale).Write("throw_potscale", TObject::kOverwrite);
}

// put a closed temporary under the final name
inline void throw_install(const TString& ftmp, const TString& fin, const std::vector<TString>& inputs, int seed, double potscale){
  if (gSystem->Link(ftmp, fin) == 0) { gSystem->Unlink(ftmp); return; }
  // the name exists: created meanwhile by another job (identical content), or stale
  if (throw_reusable(fin, inputs, seed, potscale)) { gSystem->Unlink(ftmp); return; }
  if (gSystem->Rename(ftmp, fin) != 0) { printf("FAILED rename %s\n", ftmp.Data()); gSystem->Unlink(ftmp); }
}
