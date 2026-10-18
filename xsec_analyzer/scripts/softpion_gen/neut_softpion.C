// neut_softpion.C -- softest charged-pion momentum of the two-pion signal in a NEUT vector file, as nuwro_softpion.cc
// (CC: |Mode| < 30; final-state particles fIsAlive && fStatus == 0). Run with run_nuwro_neut.sh (ROOT 5).
#include <cmath>
#include <cstdio>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "neutvect.h"

static bool is_meson(int pdg){ int a=std::abs(pdg);
  if(a>=9900000||(a/1000)%10!=0||(a/100)%10==0) return false;
  if((a>=901&&a<=930)||a==110||a==990||a==998||a==999||a==100) return false; return true; }

void neut_softpion(const char* infile, const char* outfile){
  const double E[12]={0.0,0.05,0.10,0.125,0.15,0.175,0.20,0.25,0.30,0.40,0.60,1e9};
  double h[11]={0};
  TFile f(infile); TTree* t=(TTree*)f.Get("neuttree");
  NeutVect* nv=new NeutVect(); t->SetBranchAddress("vectorbranch",&nv);
  long N=t->GetEntries(); double tot=0;
  for(long i=0;i<N;i++){ t->GetEntry(i);
    if(std::abs(nv->Mode)>=30) continue;
    int nmu=0,npi=0,npi0=0,nk=0,nh=0; double pmu=0,psoft=1e9;
    for(int j=0;j<nv->Npart();j++){ NeutPart* p=nv->PartInfo(j);
      if(!p->fIsAlive||p->fStatus!=0) continue;
      int pdg=p->fPID,a=std::abs(pdg); double m=TVector3(p->fP.X(),p->fP.Y(),p->fP.Z()).Mag()/1000.;
      if(a==13){nmu++;pmu=m;} else if(a==211){npi++; if(m<psoft) psoft=m;} else if(pdg==111) npi0++;
      else if(a==321||a==311||pdg==310||pdg==130) nk++; else if(is_meson(pdg)) nh++; }
    if(nmu!=1||pmu<=0.15||npi!=2||npi0||nk||nh) continue;
    for(int b=0;b<11;b++) if(psoft>=E[b]&&psoft<E[b+1]){ h[b]+=1; tot+=1; break; } }
  FILE* o=fopen(outfile,"w");
  fprintf(o,"{\"edges\": [0.0, 0.05, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.60, 1e9], \"counts\": [");
  for(int b=0;b<11;b++) fprintf(o,"%s%.1f",b?", ":"",h[b]);
  fprintf(o,"], \"n_events\": %ld}\n",N); fclose(o);
  printf("%s: %ld events, two-pion %.0f\n",outfile,N,tot);
}
