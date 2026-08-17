// gibuu_cc1pi.C — apply the CC1mu1piXp signal to a GiBUU FinalEvents.dat and
// produce per-nucleon dsigma/dx in the analysis binning (10^-38 cm^2/nucleon).
//
// FinalEvents.dat columns (1-based):
//   1 Run  2 Event  3 ID  4 Charge  5 perweight  6-8 position  9-12 momentum(0..3)
//   13 history  14 production_ID  15 enu
// GiBUU IDs: 902 = muon; 101 = pion (charge ±1 charged, 0 = pi0);
//   110 = K+/K0, 111 = K-/K0bar; 102-109,112-122 = heavier mesons.
// perweight is the per-nucleon cross section contribution in 10^-38 cm^2; the
// per-nucleon sigma for a selection = sum(perweight of the muon over selected
// events) / num_runs_SameEnergy. Differential: fill weighted by perweight,
// divide by num_runs and bin width.
//
// Run: root -l -b -q 'gibuu_cc1pi.C("FinalEvents.dat", NUM_RUNS, "out.root")'

#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TVector3.h"
#include "TFile.h"

static bool is_heavy_meson_gibuu(int id){
  // mesons 101..122; exclude pion(101) and kaons(110,111)
  return (id>=102 && id<=109) || (id>=112 && id<=122);
}

void gibuu_cc1pi_pmufix(const char* infile, double num_runs, const char* outfile){
  std::vector<double> EPMU={0.150,0.350,0.550,0.750,0.950,1.250,1.750,99.0};
  std::vector<double> EPPI={0.175,0.250,0.320,0.420,0.550,1.000};
  std::vector<double> ECMU={-1.0,0.45,0.65,0.80,0.90,1.0};
  std::vector<double> ECPI={-1.0,-0.10,0.35,0.55,0.75,1.0};
  std::vector<double> ETH={0.0,0.60,0.85,1.10,1.30,1.52,1.85,2.60};
  TH1D* hpmu=new TH1D("pmu","",EPMU.size()-1,EPMU.data());
  TH1D* hppi=new TH1D("ppi","",EPPI.size()-1,EPPI.data());
  TH1D* hcmu=new TH1D("costhmu","",ECMU.size()-1,ECMU.data());
  TH1D* hcpi=new TH1D("costhpi","",ECPI.size()-1,ECPI.data());
  TH1D* hth =new TH1D("thmupi","",ETH.size()-1,ETH.data());

  FILE* f=fopen(infile,"r");
  if(!f){ printf("cannot open %s\n",infile); return; }
  char line[1024];
  // group by (run,event): accumulate this event's particles
  long cur_run=-1, cur_evt=-1;
  int n_mu=0,n_pipm=0,n_pi0=0,n_kaon=0,n_heavy=0;
  double w_evt=0; TVector3 p_mu,p_pi;
  long n_events=0, n_signal=0;

  auto finish_event=[&](){
    if(cur_evt<0) return;
    ++n_events;
    if(n_mu==1 && n_pipm==1 && n_pi0==0 && n_kaon==0 && n_heavy==0){
      double pmu=p_mu.Mag(), ppi=p_pi.Mag(), th=p_mu.Angle(p_pi);
      if(pmu>0.15 && ppi>0.175 && th<2.6){
        ++n_signal;
        hpmu->Fill(pmu,w_evt); hppi->Fill(ppi,w_evt);
        hcmu->Fill(p_mu.CosTheta(),w_evt); hcpi->Fill(p_pi.CosTheta(),w_evt);
        hth->Fill(th,w_evt);
      }
    }
  };

  while(fgets(line,sizeof(line),f)){
    if(line[0]=='#') continue;
    long run,evt; int id,charge; double pw,x,y,z,e,px,py,pz;
    long hist; int prod; double enu;
    int nr=sscanf(line,"%ld %ld %d %d %lf %lf %lf %lf %lf %lf %lf %lf %ld %d %lf",
      &run,&evt,&id,&charge,&pw,&x,&y,&z,&e,&px,&py,&pz,&hist,&prod,&enu);
    if(nr<12) continue;
    if(run!=cur_run || evt!=cur_evt){
      finish_event();
      cur_run=run; cur_evt=evt;
      n_mu=n_pipm=n_pi0=n_kaon=n_heavy=0; w_evt=0;
    }
    if(id==902){ n_mu++; p_mu.SetXYZ(px,py,pz); w_evt=pw; }
    else if(id==101){ if(charge!=0){ n_pipm++; p_pi.SetXYZ(px,py,pz);} else n_pi0++; }
    else if(id==110||id==111){ n_kaon++; }
    else if(is_heavy_meson_gibuu(id)){ n_heavy++; }
  }
  finish_event();
  fclose(f);

  // Normalise: per-nucleon dsigma/dx = (sum perweight)/num_runs/binwidth.
  // perweight is in 10^-38 cm^2, but combine_newg4.C expects ABSOLUTE cm^2 (as the
  // NuWro/NEUT predictions provide) since it divides by 1e-38 itself -- so multiply
  // by 1e-38 here to match units. (Without this the combined result is 1e38x high.)
  auto norm=[&](TH1D* h){ for(int b=1;b<=h->GetNbinsX();++b){
    double c=h->GetBinContent(b), w=h->GetBinWidth(b);
    h->SetBinContent(b, c/num_runs/w*1e-38);
    h->SetBinError(b, h->GetBinError(b)/num_runs/w*1e-38); } };
  norm(hpmu);norm(hppi);norm(hcmu);norm(hcpi);norm(hth);

  TFile fo(outfile,"recreate");
  hpmu->Write();hppi->Write();hcmu->Write();hcpi->Write();hth->Write(); fo.Close();
  printf("GiBUU: events=%ld CC1pi-signal=%ld (%.2f%%)  num_runs=%.0f\n",
    n_events,n_signal,100.*n_signal/n_events,num_runs);
  printf("  wrote %s\n",outfile);
}
