// gibuu_cc1pi_1p.C — proton-tagged (CC1mu1pi1p) W/TKI prediction from a GiBUU
// FinalEvents.dat. Inclusive signal (gibuu_cc1pi.C) plus a leading true proton
// (GiBUU id==1, charge==1) above 0.3 GeV/c; six W/TKI observables (wtki_gen.h),
// perweight-weighted. Units: perweight is 10^-38 cm^2, multiplied by 1e-38 to give
// absolute cm^2 (combine_wtki.C divides by 1e-38 like the other generators).
//   root -l -b -q 'gibuu_cc1pi_1p.C("run_numu/FinalEvents.dat",20,"gibuu_1p_numu.root")'
#include <cstdio>
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TVector3.h"
#include "TFile.h"
#include "../../wtki_gen.h"
static bool is_heavy_meson_gibuu(int id){ return (id>=102&&id<=109)||(id>=112&&id<=122); }

void gibuu_cc1pi_1p(const char* infile, double num_runs, const char* outfile){
  const char* OBS[6]={"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  TH1D* h[6]; for(int o=0;o<6;o++){ auto e=wtki::edges(OBS[o]); h[o]=new TH1D(OBS[o],"",e.size()-1,e.data()); }
  FILE* f=fopen(infile,"r"); if(!f){ printf("cannot open %s\n",infile); return; }
  char line[1024];
  long cur_run=-1,cur_evt=-1; int n_mu=0,n_pipm=0,n_pi0=0,n_kaon=0,n_heavy=0;
  double w_evt=0,leadp=0; TVector3 p_mu,p_pi,p_pr;
  long n_events=0,n_signal=0;
  auto finish_event=[&](){
    if(cur_evt<0) return; ++n_events;
    if(n_mu==1&&n_pipm==1&&n_pi0==0&&n_kaon==0&&n_heavy==0 && leadp>0.3){
      double pmu=p_mu.Mag(),ppi=p_pi.Mag(),th=p_mu.Angle(p_pi);
      if(pmu>0.15&&ppi>0.175&&th<2.6){
        ++n_signal;
        wtki::Obs w=wtki::compute(p_mu,p_pi,p_pr);
        double v[6]={w.Wpipr,w.Whad,w.dpt,w.dalphat,w.dphit,w.pn};
        for(int o=0;o<6;o++){ auto e=wtki::edges(OBS[o]); h[o]->Fill(wtki::clamp(v[o],e),w_evt); }
      }
    }
  };
  while(fgets(line,sizeof(line),f)){
    if(line[0]=='#') continue;
    long run,evt; int id,charge; double pw,x,y,z,e,px,py,pz; long hist; int prod; double enu;
    int nr=sscanf(line,"%ld %ld %d %d %lf %lf %lf %lf %lf %lf %lf %lf %ld %d %lf",
      &run,&evt,&id,&charge,&pw,&x,&y,&z,&e,&px,&py,&pz,&hist,&prod,&enu);
    if(nr<12) continue;
    if(run!=cur_run||evt!=cur_evt){ finish_event(); cur_run=run; cur_evt=evt;
      n_mu=n_pipm=n_pi0=n_kaon=n_heavy=0; w_evt=0; leadp=0; }
    if(id==902){ n_mu++; p_mu.SetXYZ(px,py,pz); w_evt=pw; }
    else if(id==101){ if(charge!=0){ n_pipm++; p_pi.SetXYZ(px,py,pz);} else n_pi0++; }
    else if(id==1&&charge==1){ TVector3 p(px,py,pz); if(p.Mag()>leadp){ leadp=p.Mag(); p_pr=p; } }
    else if(id==110||id==111){ n_kaon++; }
    else if(is_heavy_meson_gibuu(id)){ n_heavy++; }
  }
  finish_event(); fclose(f);
  auto norm=[&](TH1D* hh){ for(int b=1;b<=hh->GetNbinsX();++b){
    double c=hh->GetBinContent(b),w=hh->GetBinWidth(b);
    hh->SetBinContent(b,c/num_runs/w*1e-38); hh->SetBinError(b,hh->GetBinError(b)/num_runs/w*1e-38); } };
  for(int o=0;o<6;o++) norm(h[o]);
  TFile fo(outfile,"recreate"); for(int o=0;o<6;o++) h[o]->Write(); fo.Close();
  printf("GiBUU 1p: events=%ld  proton-tagged W/TKI signal=%ld (%.3f%%)  wrote %s\n",
    n_events,n_signal,100.*n_signal/n_events,outfile);
}
