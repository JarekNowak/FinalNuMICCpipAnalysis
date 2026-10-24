// nuwro_2pi.cc -- two-pion signal events of a NuWro output file (post-FSI) for scripts/mp2pi_gen/mp2pi_gen.py:
// CC, exactly one muon, exactly two charged pions (any momentum), no pi0, kaon or heavier meson. One line per such
// event: muon momentum vector, the two pion momentum vectors (GeV/c), the number of protons above 0.3 GeV/c, the
// neutrino energy (GeV); the
// header carries the number of events and the flux-averaged cross section per nucleon (cm2, the event weight).
#include <cmath>
#include <cstdio>
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
  if (argc < 3) { fprintf(stderr, "usage: %s in.root out.txt\n", argv[0]); return 1; }
  TFile fin(argv[1]); TTree* t = (TTree*)fin.Get("treeout");
  event* e = new event(); t->SetBranchAddress("e", &e);
  long N = t->GetEntries(); double sigma = 0.; long nsig = 0;
  FILE* f = fopen(argv[2], "w");
  fprintf(f, "# N %ld\n", N);
  for (long i = 0; i < N; ++i) {
    t->GetEntry(i); sigma = e->weight;
    int nmu = 0, npi = 0, npi0 = 0, nk = 0, nh = 0, np = 0; TVector3 mu, pi[2];
    for (size_t j = 0; j < e->post.size(); ++j) {
      int pdg = e->post[j].pdg, a = std::abs(pdg);
      TVector3 p(e->post[j].x / 1000., e->post[j].y / 1000., e->post[j].z / 1000.);
      if (a == 13 && p.Mag() > 0.) { ++nmu; mu = p; }
      else if (a == 211 && p.Mag() > 0.) { if (npi < 2) pi[npi] = p; ++npi; }
      else if (pdg == 111) ++npi0;
      else if (a == 321 || a == 311 || pdg == 310 || pdg == 130) ++nk;
      else if (is_meson(pdg)) ++nh;
      else if (pdg == 2212 && p.Mag() > 0.3) ++np;
    }
    if (nmu != 1 || npi != 2 || npi0 || nk || nh) continue;
    ++nsig;
    fprintf(f, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %d %.6g\n", mu.X(), mu.Y(), mu.Z(), pi[0].X(), pi[0].Y(), pi[0].Z(),
            pi[1].X(), pi[1].Y(), pi[1].Z(), np, e->in[0].t / 1000.);
  }
  fprintf(f, "# sigma %.6e\n", sigma); fclose(f);
  printf("%s: %ld events, %ld two-pion, sigma %.4e\n", argv[2], N, nsig, sigma);
  return 0;
}
