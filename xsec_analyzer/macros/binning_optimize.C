// binopt.C — measure per-observable resolution & occupancy, then greedily design
// variable-width bins that are >= k*resolution wide and carry >= Nmin selected
// signal events, to reduce the unfolding-amplified and statistical systematics.
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

struct Obs { std::string n, tb, rb; double lo, hi; bool frac; };

void binopt() {
  const char* FN = "/data/uboone/processed/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root";
  TFile f(FN); TTree* t = (TTree*)f.Get("stv_tree");
  std::vector<Obs> O = {
    {"pmu",   "CC1mu1piXp_candidate_muon_mom_true","CC1mu1piXp_candidate_muon_mom_reco", 0.15,3.0, true},
    {"ppi",   "CC1mu1piXp_candidate_pion_mom_true","CC1mu1piXp_candidate_pion_mom_reco", 0.175,1.0,true},
    {"costhmu","CC1mu1piXp_candidate_muon_costh_true","CC1mu1piXp_candidate_muon_costh_reco",-1,1,false},
    {"costhpi","CC1mu1piXp_candidate_pion_costh_true","CC1mu1piXp_candidate_pion_costh_reco",-1,1,false},
    {"thmupi","CC1mu1piXp_true_mu_pi_opening_angle","CC1mu1piXp_mu_pi_opening_angle",0,2.6,false},
  };
  // enable only needed branches (fast read of the 15.7 GB file)
  t->SetBranchStatus("*",0);
  Bool_t sel, sig; double tv[5], rv[5];
  t->SetBranchStatus("CC1mu1piXp_Selected",1); t->SetBranchAddress("CC1mu1piXp_Selected",&sel);
  t->SetBranchStatus("CC1mu1piXp_MC_Signal",1); t->SetBranchAddress("CC1mu1piXp_MC_Signal",&sig);
  for (int o=0;o<5;o++){ t->SetBranchStatus(O[o].tb.c_str(),1); t->SetBranchAddress(O[o].tb.c_str(),&tv[o]);
                         t->SetBranchStatus(O[o].rb.c_str(),1); t->SetBranchAddress(O[o].rb.c_str(),&rv[o]); }
  const int NF=200;
  std::vector<TH1D*> occ, sw, sw2; // occupancy, sum(reco-true), sum((reco-true)^2) per fine true bin
  for (auto&ob:O){ occ.push_back(new TH1D(("occ_"+ob.n).c_str(),"",NF,ob.lo,ob.hi));
                   sw.push_back(new TH1D(("sw_"+ob.n).c_str(),"",NF,ob.lo,ob.hi));
                   sw2.push_back(new TH1D(("sw2_"+ob.n).c_str(),"",NF,ob.lo,ob.hi)); }
  Long64_t N=t->GetEntries();
  double POT_SCALE = 3.283e20/2.328199e21; // data/MC POT (approx, for occupancy scale)
  for (Long64_t i=0;i<N;i++){ t->GetEntry(i); if(!(sel&&sig)) continue;
    for(int o=0;o<5;o++){ double tr=tv[o], re=rv[o]; if(tr<O[o].lo||tr>=O[o].hi) continue;
      double d=re-tr; occ[o]->Fill(tr); sw[o]->Fill(tr,d); sw2[o]->Fill(tr,d*d); } }
  printf("entries=%lld\n",N);
  // greedy binning per observable
  for(int o=0;o<5;o++){ auto&ob=O[o];
    double totsel = occ[o]->Integral()*POT_SCALE;
    printf("\n===== %s  [%.3f,%.3f]  selected-signal(dataPOT)=%.0f  frac=%d =====\n",
           ob.n.c_str(),ob.lo,ob.hi,totsel,ob.frac);
    // resolution profile (coarse print): 8 slices
    printf("  resolution by slice: ");
    for(int s=0;s<8;s++){ double a=ob.lo+(ob.hi-ob.lo)*s/8, b=ob.lo+(ob.hi-ob.lo)*(s+1)/8;
      int ba=occ[o]->FindBin(a+1e-6), bb=occ[o]->FindBin(b-1e-6); double n=0,S=0,S2=0;
      for(int k=ba;k<=bb;k++){n+=occ[o]->GetBinContent(k);S+=sw[o]->GetBinContent(k);S2+=sw2[o]->GetBinContent(k);}
      double rms = n>1? sqrt(std::max(0.0,S2/n-(S/n)*(S/n))):0;
      printf("[%.2f]%.3f ",(a+b)/2,rms); }
    printf("\n");
    // precompute a smoothed LOCAL resolution profile res_s[fine bin] using a
    // sliding window (occupancy-weighted), so variable-resolution observables
    // (e.g. costhmu: 0.86 backward, 0.10 forward) get local, not global, widths.
    double step=(ob.hi-ob.lo)/NF; std::vector<double> res_s(NF,0.);
    int W=12; // window half-width in fine bins
    for(int b=0;b<NF;b++){ double nn=0,SS=0,SS2=0;
      for(int k=std::max(1,b+1-W);k<=std::min(NF,b+1+W);k++){
        nn+=occ[o]->GetBinContent(k);SS+=sw[o]->GetBinContent(k);SS2+=sw2[o]->GetBinContent(k);}
      res_s[b]= nn>1? sqrt(std::max(0.0,SS2/nn-(SS/nn)*(SS/nn))):0; }
    // greedy: grow bin until width>=k*res_LOCAL(low edge) AND occ>=Nmin
    double kres = 1.0;              // bin width >= 1x local RMS resolution
    double Nmin = 45.0;             // min selected-signal events per bin (~15% stat target)
    std::vector<double> edges={ob.lo}; double cur=ob.lo;
    while(cur<ob.hi-1e-6){ int lb=std::min(NF-1,std::max(0,int((cur-ob.lo)/step)));
      double reslocal=res_s[lb]; double e=cur; double occacc=0;
      while(e<ob.hi-1e-6){ double c=e+step/2; int fb=occ[o]->FindBin(c);
        occacc+=occ[o]->GetBinContent(fb)*POT_SCALE; e+=step;
        if((e-cur)>=kres*reslocal && occacc>=Nmin) break; }
      if(e>ob.hi) e=ob.hi;
      edges.push_back(e); cur=e;
    }
    // merge a too-small FINAL bin into its neighbour (occurred at the tail only)
    int ne=edges.size();
    if(ne>=3){ double last=edges[ne-1]-edges[ne-2], prev=edges[ne-2]-edges[ne-3];
      if(last<0.5*prev) edges.erase(edges.end()-2); }
    // print proposed edges
    printf("  PROPOSED %d bins: ", (int)edges.size()-1);
    for(double x:edges) printf("%.3f ",x); printf("\n");
    printf("  (current had %s)\n", o==0?"22":o==1?"5":o==4?"9":"12");
  }
}
