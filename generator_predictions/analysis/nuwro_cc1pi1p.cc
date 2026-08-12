// nuwro_cc1pi1p.cc — proton-tagged (CC1mu1pi1p) W/TKI prediction from a NuWro output
// file (post-FSI). Inclusive signal (nuwro_cc1pi.cc) plus a leading true proton above
// 0.3 GeV/c; six W/TKI observables (wtki_gen.h), per-nucleon dsigma/dx in absolute cm^2
// (combine_wtki.C rescales to per-Ar 1e-38).
//   build (in the NuWro env): g++ nuwro_cc1pi1p.cc -o nuwro_cc1pi1p -I$NUWRO/src \
//     $(root-config --cflags --libs) -lEG $NUWRO_FQ_DIR/bin/event1.so
//   run:   ./nuwro_cc1pi1p out_numu.root pred_1p_numu.root
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
#include "../wtki_gen.h"

static bool is_meson(int pdg){ int a=std::abs(pdg);
  if(a>=9900000)return false; if((a/1000)%10!=0)return false; if((a/100)%10==0)return false;
  if(a>=901&&a<=930)return false; if(a==110||a==990||a==998||a==999||a==100)return false; return true; }
static bool is_kaon(int pdg){ int a=std::abs(pdg); return a==321||a==311||pdg==310||pdg==130; }

int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s in.root out.root\n",argv[0]); return 1; }
  TFile fin(argv[1]); if(fin.IsZombie()){ fprintf(stderr,"cannot open %s\n",argv[1]); return 1; }
  TTree* t=(TTree*)fin.Get("treeout"); if(!t){ fprintf(stderr,"no treeout\n"); return 1; }
  event* e=new event(); t->SetBranchAddress("e",&e);
  const char* OBS[6]={"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  TH1D* h[6]; for(int o=0;o<6;o++){ std::vector<double> ed=wtki::edges(OBS[o]); h[o]=new TH1D(OBS[o],OBS[o],(int)ed.size()-1,ed.data()); }
  long N=t->GetEntries(); double sigma_tot=0.; long n_signal=0;
  for(long i=0;i<N;i++){ t->GetEntry(i); sigma_tot=e->weight;
    int n_mu=0,n_pipm=0,n_pi0=0,n_kaon=0,n_heavy=0; TVector3 p_mu,p_pi,p_pr; double leadp=0.;
    for(size_t j=0;j<e->post.size();j++){ int pdg=e->post[j].pdg,a=std::abs(pdg);
      TVector3 p(e->post[j].x/1000.,e->post[j].y/1000.,e->post[j].z/1000.); double mom=p.Mag();
      if(a==13&&mom>0.){ n_mu++; p_mu=p; }
      else if(a==211&&mom>0.){ n_pipm++; p_pi=p; }
      else if(pdg==111){ n_pi0++; }
      else if(pdg==2212){ if(mom>leadp){ leadp=mom; p_pr=p; } }
      else if(is_kaon(pdg)){ n_kaon++; }
      else if(a!=111&&a!=211&&is_meson(pdg)){ n_heavy++; } }
    if(n_mu!=1||n_pipm!=1||n_pi0!=0||n_kaon!=0||n_heavy!=0) continue;
    if(leadp<=0.3) continue;
    double pmu=p_mu.Mag(),ppi=p_pi.Mag(),th=p_mu.Angle(p_pi);
    if(pmu<=0.15||ppi<=0.175||th>=2.6) continue;
    ++n_signal;
    wtki::Obs w=wtki::compute(p_mu,p_pi,p_pr);
    double v[6]={w.Wpipr,w.Whad,w.dpt,w.dalphat,w.dphit,w.pn};
    for(int o=0;o<6;o++){ std::vector<double> ed=wtki::edges(OBS[o]); h[o]->Fill(wtki::clamp(v[o],ed)); }
  }
  for(int o=0;o<6;o++){ for(int b=1;b<=h[o]->GetNbinsX();b++){ double c=h[o]->GetBinContent(b),w=h[o]->GetBinWidth(b);
    h[o]->SetBinContent(b,sigma_tot*(c/(double)N)/w); h[o]->SetBinError(b,sigma_tot*(std::sqrt(c)/(double)N)/w); } }
  TFile fo(argv[2],"recreate"); for(int o=0;o<6;o++) h[o]->Write(); fo.Close();
  printf("nuwro_cc1pi1p: N=%ld  proton-tagged W/TKI signal=%ld (%.3f%%)  wrote %s\n",N,n_signal,100.*n_signal/(double)N,argv[2]);
  return 0;
}
