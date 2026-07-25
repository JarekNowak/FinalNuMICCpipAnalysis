// threshcheck.C — validate the momentum thresholds in the signal definition.
// Builds a THRESHOLD-FREE signal (sig_ component flags, no pmu/ppi/angle cut) and
// measures (a) the selection efficiency turn-on vs true momentum, and (b) the
// reco-vs-true bias, near the pmu>0.15 and ppi>0.175 GeV/c thresholds.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <cmath>

void threshcheck() {
  TFile f("/data/uboone/processed/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t = (TTree*)f.Get("stv_tree");
  t->SetBranchStatus("*",0);
  Bool_t sel, fv, ccnc, isnumu, one_mu, one_pi, no_pi0, no_hm;
  Double_t mut, mur, pit, pir;
  auto B=[&](const char* n, void* p){ t->SetBranchStatus(n,1); t->SetBranchAddress(n,p); };
  B("CC1mu1piXp_Selected",&sel);
  B("CC1mu1piXp_sig_truevertex_in_fv",&fv); B("CC1mu1piXp_sig_ccnc",&ccnc);
  B("CC1mu1piXp_sig_is_numu",&isnumu); B("CC1mu1piXp_sig_one_muon_above_thresh",&one_mu);
  B("CC1mu1piXp_sig_one_charged_pion",&one_pi); B("CC1mu1piXp_sig_no_pions",&no_pi0);
  B("CC1mu1piXp_sig_no_heavy_mesons",&no_hm);
  B("CC1mu1piXp_candidate_muon_mom_true",&mut); B("CC1mu1piXp_candidate_muon_mom_reco",&mur);
  B("CC1mu1piXp_candidate_pion_mom_true",&pit); B("CC1mu1piXp_candidate_pion_mom_reco",&pir);

  int NB=60; double lo=0.0, hi=0.6;
  TH1D hmd("hmd","",NB,lo,hi), hmn("hmn","",NB,lo,hi); // muon den/num vs true pmu
  TH1D hpd("hpd","",NB,lo,hi), hpn("hpn","",NB,lo,hi); // pion den/num vs true ppi
  TH1D bmu_s("bmu_s","",NB,lo,hi), bmu_n("bmu_n","",NB,lo,hi); // muon bias sum/count
  TH1D bpi_s("bpi_s","",NB,lo,hi), bpi_n("bpi_n","",NB,lo,hi); // pion bias sum/count

  Long64_t N=t->GetEntries();
  for (Long64_t i=0;i<N;i++){ t->GetEntry(i);
    bool tf = fv && ccnc && isnumu && one_mu && one_pi && no_pi0 && no_hm; // threshold-free signal
    if(!tf) continue;
    if(mut>lo&&mut<hi){ hmd.Fill(mut); if(sel){ hmn.Fill(mut); bmu_s.Fill(mut,mur-mut); bmu_n.Fill(mut,1);} }
    if(pit>lo&&pit<hi){ hpd.Fill(pit); if(sel){ hpn.Fill(pit); bpi_s.Fill(pit,pir-pit); bpi_n.Fill(pit,1);} }
  }
  printf("entries=%lld\n",N);
  auto report=[&](const char* nm, TH1D& d, TH1D& n, TH1D& bs, TH1D& bc, double thr){
    printf("\n===== %s  (phase-space threshold = %.3f GeV/c) =====\n", nm, thr);
    printf("  true mom | eff(sel/threshfree) | mean(reco-true) | Ndenom\n");
    for(int b=1;b<=NB;b++){ double c=d.GetBinCenter(b); if(c<0.05||c>0.40) continue;
      double D=d.GetBinContent(b), Nn=n.GetBinContent(b);
      double eff = D>0? Nn/D : 0;
      double bias = bc.GetBinContent(b)>0? bs.GetBinContent(b)/bc.GetBinContent(b):0;
      printf("   %.3f    |   %5.1f%%           |   %+.3f        | %.0f%s\n",
        c, 100*eff, bias, D, (fabs(c-thr)<d.GetBinWidth(1)?"  <-- threshold":""));
    }
  };
  report("MUON pmu", hmd,hmn,bmu_s,bmu_n, 0.15);
  report("PION ppi", hpd,hpn,bpi_s,bpi_n, 0.175);
}
