// mares_reweight.C — apply an M_A^RES +/-1sigma reweight to the fake data in a
// COPY of the costhmu univmake, so unfolding produces the M_A^RES-shifted result.
// The per-true-bin factor f(t) (signal costhmu bins) comes from the standalone
// GENIE MaCCRES excursion (genie_mares_{nom,<var>}.root). The measured reco
// (onBNB_reco) is reweighted by the migration-propagated factor
// g(r) = sum_t mig(t,r) f(t) / sum_t mig(t,r); the closure truth (FakeDataMC
// true / 2d) by f(t). The catch-all (background) true bin is left unchanged.
//
//   root -l -b -q 'mares_reweight.C("copy.root","plus")'   // or "minus"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

void mares_reweight(const char* univ_copy, const char* variation) {
  const char* GENDIR = "/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4";
  TFile fnf(Form("%s/genie_mares_nom.root", GENDIR));
  TFile fvf(Form("%s/genie_mares_%s.root", GENDIR, variation));
  TH1D* hn = (TH1D*)fnf.Get("costhmu");
  TH1D* hv = (TH1D*)fvf.Get("costhmu");
  double f[14];
  for (int t = 1; t <= 12; t++)
    f[t] = (hn->GetBinContent(t) > 0.) ? hv->GetBinContent(t)/hn->GetBinContent(t) : 1.0;
  f[13] = 1.0;                       // catch-all (background) true bin unchanged
  fnf.Close(); fvf.Close();

  TFile fu(univ_copy, "UPDATE");
  TDirectoryFile* tot = (TDirectoryFile*)fu.Get(
    "ccpi_CC1mu1piXp_costhmu_1D/total_configs+file_properties_numi.txt");
  TH2D* mig = (TH2D*)tot->Get("FakeDataMC_0_2d");   // X=true(13), Y=reco(12)
  int NT = mig->GetNbinsX(), NR = mig->GetNbinsY();
  double g[64];
  for (int r = 1; r <= NR; r++) {
    double num = 0, den = 0;
    for (int t = 1; t <= NT; t++) { double m = mig->GetBinContent(t, r); num += m*f[t]; den += m; }
    g[r] = (den > 0.) ? num/den : 1.0;
  }
  // measured reco (drives the unfold) + fake-data-universe reco: reweight by g(r)
  for (const char* nm : {"onBNB_reco", "FakeDataMC_0_reco"}) {
    TH1D* h = (TH1D*)tot->Get(nm);
    if (!h) continue;
    for (int r = 1; r <= h->GetNbinsX(); r++) h->SetBinContent(r, h->GetBinContent(r)*g[r]);
    tot->cd(); h->Write(nm, TObject::kOverwrite);
  }
  // closure truth: reweight by f(t)
  TH1D* ht = (TH1D*)tot->Get("FakeDataMC_0_true");
  for (int t = 1; t <= ht->GetNbinsX(); t++) ht->SetBinContent(t, ht->GetBinContent(t)*f[t]);
  tot->cd(); ht->Write("FakeDataMC_0_true", TObject::kOverwrite);
  // migration: reweight the true axis by f(t)
  for (int t = 1; t <= NT; t++) for (int r = 1; r <= NR; r++)
    mig->SetBinContent(t, r, mig->GetBinContent(t, r)*f[t]);
  tot->cd(); mig->Write("FakeDataMC_0_2d", TObject::kOverwrite);
  fu.Close();
  printf("reweighted %s with M_A^RES %s: g(reco)[1,6,12]=%.3f,%.3f,%.3f\n",
    univ_copy, variation, g[1], g[6], g[12]);
}
