// nuwro_softpion.cc -- softest charged-pion momentum of the two-pion signal in a NuWro output file (post-FSI), as
// softpion_gst_gibuu.py: CC numu/numubar sample, one muon above 0.15 GeV/c, exactly two charged pions, no pi0, kaon or
// heavier meson; counts per bin of GEDGES written as JSON. Build and run with run_nuwro_neut.sh (SL7 container).
#include <cmath>
#include <cstdio>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "event1.h"

static bool is_meson(int pdg) {
  int a = std::abs(pdg);
  if (a >= 9900000 || (a / 1000) % 10 != 0 || (a / 100) % 10 == 0) return false;
  if ((a >= 901 && a <= 930) || a == 110 || a == 990 || a == 998 || a == 999 || a == 100) return false;
  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s in.root out.json\n", argv[0]); return 1; }
  const std::vector<double> E = {0.0, 0.05, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.60, 1e9};
  std::vector<double> h(E.size() - 1, 0.);
  TFile fin(argv[1]); TTree* t = (TTree*)fin.Get("treeout");
  event* e = new event(); t->SetBranchAddress("e", &e);
  long N = t->GetEntries();
  for (long i = 0; i < N; ++i) {
    t->GetEntry(i);
    int nmu = 0, npi = 0, npi0 = 0, nk = 0, nh = 0; double pmu = 0., psoft = 1e9;
    for (size_t j = 0; j < e->post.size(); ++j) {
      int pdg = e->post[j].pdg, a = std::abs(pdg);
      double p = TVector3(e->post[j].x, e->post[j].y, e->post[j].z).Mag() / 1000.;
      if (a == 13 && p > 0.) { ++nmu; pmu = p; }
      else if (a == 211 && p > 0.) { ++npi; if (p < psoft) psoft = p; }
      else if (pdg == 111) ++npi0;
      else if (a == 321 || a == 311 || pdg == 310 || pdg == 130) ++nk;
      else if (is_meson(pdg)) ++nh;
    }
    if (nmu != 1 || pmu <= 0.15 || npi != 2 || npi0 || nk || nh) continue;
    for (size_t b = 0; b + 1 < E.size(); ++b) if (psoft >= E[b] && psoft < E[b + 1]) { h[b] += 1.; break; }
  }
  FILE* f = fopen(argv[2], "w");
  fprintf(f, "{\"edges\": [0.0, 0.05, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.60, 1e9], \"counts\": [");
  for (size_t b = 0; b < h.size(); ++b) fprintf(f, "%s%.1f", b ? ", " : "", h[b]);
  fprintf(f, "], \"n_events\": %ld}\n", N); fclose(f);
  printf("%s: %ld events, two-pion %.0f\n", argv[2], N, [&]{ double s = 0; for (double x : h) s += x; return s; }());
  return 0;
}
