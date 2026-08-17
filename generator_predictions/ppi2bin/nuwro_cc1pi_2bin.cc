// nuwro_cc1pi.cc
// Apply the CC1mu1piXp signal definition to a NuWro output file at the post-FSI
// truth level and produce differential cross-section predictions in the exact
// analysis binning, in the framework's units (10^-38 cm^2 / Ar), for direct
// use as FileTrueEvents predictions.
//
// Signal (matches CC1mu1piXp::define_signal, kinematic part only -- the
// detector fiducial-volume cut does not apply to a per-nucleus generator
// prediction):
//   - numu CC                       (guaranteed by the run card: numu, CC-only)
//   - exactly 1 muon  (|pdg|==13, p>0)
//   - exactly 1 charged pion (|pdg|==211, p>0)
//   - 0 neutral pions (pdg==111)
//   - 0 kaons (|pdg| in {321,311}, or 310, 130)
//   - 0 heavier mesons (is_meson, excluding pi0/pi+-)
//   - any number of protons
//   - phase space: p_mu > 0.15 GeV/c, p_pi > 0.175 GeV/c, theta_mupi < 2.6 rad
//
// Observables (match compute_true_observables): muon & pion momentum magnitude
// (GeV/c), cos(theta) w.r.t. the beam (+z), and the 3D mu-pi opening angle.
//
// Normalisation: NuWro sets every event's weight to the flux-averaged total
// cross section sigma_tot, PER NUCLEON (confirmed empirically: the total CC
// sigma is the same for carbon A=12 and argon A=40, ratio 1.04). This analyzer
// therefore outputs the per-nucleon differential cross section
//   dsigma/dx|bin = sigma_tot * (count_in_bin / N) / dx    [cm^2 / nucleon / x]
// The conversion to the framework's per-nucleus, combined-flux, 10^-38 units is
// done in combine_nuwro.C, which weights the numu and numubar samples by their
// flux fractions and multiplies by A = 40 (nucleons per Ar).
//
// Build (inside the SL7 container, after sourcing setupNuWro.sh):
//   g++ nuwro_cc1pi.cc -o nuwro_cc1pi \
//       -I$NUWRO/src $(root-config --cflags --libs) -lEG
// Run:
//   ./nuwro_cc1pi <nuwro_out.root> <out_pred.root>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector3.h"

#include "event1.h"

// --- replicate is_meson_or_antimeson from Functions.hh -----------------------
static bool is_meson(int pdg) {
  int a = std::abs(pdg);
  if (a >= 9900000) return false;
  if ((a / 1000) % 10 != 0) return false;      // thousands digit (n_q1) must be 0
  if ((a / 100) % 10 == 0) return false;        // hundreds digit (n_q2) nonzero
  if (a >= 901 && a <= 930) return false;
  if (a == 110 || a == 990) return false;
  if (a == 998 || a == 999) return false;
  if (a == 100) return false;
  return true;
}
static bool is_kaon(int pdg) {
  int a = std::abs(pdg);
  return a == 321 || a == 311 || pdg == 310 || pdg == 130;
}

// Analysis binning (true-signal edges from configs/ccpi_*_bin_config.txt).
static const std::vector<double> EDGES_PMU = {0.150,0.350,0.550,0.750,0.950,1.250,1.750,3.000};
static const std::vector<double> EDGES_PPI = {0.175,0.205,99.0};
static const std::vector<double> EDGES_COSTHMU = {-1.0,0.45,0.65,0.80,0.90,1.0};
static const std::vector<double> EDGES_COSTHPI = {-1.0,-0.10,0.35,0.55,0.75,1.0};
static const std::vector<double> EDGES_THMUPI = {0.0,0.60,0.85,1.10,1.30,1.52,1.85,2.60};

static TH1D* make_hist(const char* name, const std::vector<double>& e) {
  return new TH1D(name, name, (int)e.size() - 1, e.data());
}

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s in.root out.root\n", argv[0]); return 1; }
  TFile fin(argv[1]);
  if (fin.IsZombie()) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  TTree* t = (TTree*)fin.Get("treeout");
  if (!t) { fprintf(stderr, "no treeout in %s\n", argv[1]); return 1; }
  event* e = new event();
  t->SetBranchAddress("e", &e);

  TH1D* h_pmu     = make_hist("pmu",     EDGES_PMU);
  TH1D* h_ppi     = make_hist("ppi",     EDGES_PPI);
  TH1D* h_costhmu = make_hist("costhmu", EDGES_COSTHMU);
  TH1D* h_costhpi = make_hist("costhpi", EDGES_COSTHPI);
  TH1D* h_thmupi  = make_hist("thmupi",  EDGES_THMUPI);

  long N = t->GetEntries();
  double sigma_tot = 0.;   // cm^2 per Ar; NuWro stores it as the event weight
  long n_signal = 0;

  for (long i = 0; i < N; ++i) {
    t->GetEntry(i);
    sigma_tot = e->weight;   // constant across events after saving

    int n_mu = 0, n_pipm = 0, n_pi0 = 0, n_kaon = 0, n_heavy = 0;
    TVector3 p_mu, p_pi;
    for (size_t j = 0; j < e->post.size(); ++j) {
      int pdg = e->post[j].pdg;
      int a = std::abs(pdg);
      // NuWro momenta are in MeV; convert to GeV.
      TVector3 p(e->post[j].x / 1000., e->post[j].y / 1000., e->post[j].z / 1000.);
      double mom = p.Mag();
      if (a == 13 && mom > 0.) { ++n_mu; p_mu = p; }
      else if (a == 211 && mom > 0.) { ++n_pipm; p_pi = p; }
      else if (pdg == 111) { ++n_pi0; }
      else if (is_kaon(pdg)) { ++n_kaon; }
      else if (a != 111 && a != 211 && is_meson(pdg)) { ++n_heavy; }
      // protons and everything else: unconstrained
    }

    // Topological signal
    if (n_mu != 1 || n_pipm != 1 || n_pi0 != 0 || n_kaon != 0 || n_heavy != 0)
      continue;

    double pmu = p_mu.Mag();
    double ppi = p_pi.Mag();
    double theta = p_mu.Angle(p_pi);

    // Phase-space cuts (ppi relaxed to 0.113 for the low bin; others standard 0.175)
    if (pmu <= 0.15 || ppi <= 0.175 || theta >= 2.6) continue;

    ++n_signal;
    h_ppi->Fill(ppi);                    // ppi extends to 0.113
    if (ppi > 0.175) {                   // other observables: standard 0.175 phase space
      h_pmu->Fill(pmu);
      h_costhmu->Fill(p_mu.CosTheta());
      h_costhpi->Fill(p_pi.CosTheta());
      h_thmupi->Fill(theta);
    }
  }

  // Normalise each histogram to the per-nucleon dsigma/dx [cm^2 / nucleon / x].
  //   dsigma/dx|bin = sigma_tot * (count/N) / binwidth
  auto normalize = [&](TH1D* h) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double count = h->GetBinContent(b);
      double w = h->GetBinWidth(b);
      double dsdx = sigma_tot * (count / (double)N) / w;
      h->SetBinContent(b, dsdx);
      double derr = sigma_tot * (std::sqrt(count) / (double)N) / w;
      h->SetBinError(b, derr);
    }
  };
  normalize(h_pmu); normalize(h_ppi); normalize(h_costhmu);
  normalize(h_costhpi); normalize(h_thmupi);

  TFile fout(argv[2], "recreate");
  h_pmu->Write(); h_ppi->Write(); h_costhmu->Write();
  h_costhpi->Write(); h_thmupi->Write();
  fout.Close();

  printf("nuwro_cc1pi: N=%ld  signal=%ld (%.3f%%)  sigma_tot=%.4e cm^2/Ar\n",
    N, n_signal, 100. * n_signal / (double)N, sigma_tot);
  printf("  wrote %s\n", argv[2]);
  return 0;
}
