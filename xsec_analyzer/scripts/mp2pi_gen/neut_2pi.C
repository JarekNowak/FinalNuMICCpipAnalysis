// neut_2pi.C -- two-pion signal events of a NEUT vector file, as nuwro_2pi.cc (last column the neutrino energy) (CC: |Mode| < 30; final-state
// particles fIsAlive && fStatus == 0). The header carries N and the flux-averaged cross section per nucleon in cm2
// (harmonic mean of Totcrs, 1e-38 cm2, as generator_predictions/gen2d/neut_incl.C).
#include <cmath>
#include <cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "neutvect.h"
static bool is_meson(int pdg){ int a=std::abs(pdg);
  if(a>=9900000||(a/1000)%10!=0||(a/100)%10==0) return false;
  if((a>=901&&a<=930)||a==110||a==990||a==998||a==999||a==100) return false; return true; }
void neut_2pi(const char* infile, const char* outfile){
  TFile f(infile); TTree* t=(TTree*)f.Get("neuttree");
  NeutVect* nv=new NeutVect(); t->SetBranchAddress("vectorbranch",&nv);
  long N=t->GetEntries(), nsig=0; double inv=0;
  FILE* o=fopen(outfile,"w"); fprintf(o,"# N %ld\n",N);
  for(long i=0;i<N;i++){ t->GetEntry(i);
    if(nv->Totcrs>0.) inv+=1.0/nv->Totcrs;
    if(std::abs(nv->Mode)>=30) continue;
    int nmu=0,npi=0,npi0=0,nk=0,nh=0,np=0; TVector3 mu,pi[2];
    for(int j=0;j<nv->Npart();j++){ NeutPart* p=nv->PartInfo(j);
      if(!p->fIsAlive||p->fStatus!=0) continue;
      int pdg=p->fPID,a=std::abs(pdg); TVector3 m(p->fP.X()/1000.,p->fP.Y()/1000.,p->fP.Z()/1000.);
      if(a==13){nmu++;mu=m;} else if(a==211){ if(npi<2) pi[npi]=m; npi++; } else if(pdg==111) npi0++;
      else if(a==321||a==311||pdg==310||pdg==130) nk++; else if(is_meson(pdg)) nh++;
      else if(pdg==2212&&m.Mag()>0.3) np++; }
    if(nmu!=1||npi!=2||npi0||nk||nh) continue;
    nsig++;
    fprintf(o,"%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %d %.6g\n",mu.X(),mu.Y(),mu.Z(),pi[0].X(),pi[0].Y(),pi[0].Z(),
            pi[1].X(),pi[1].Y(),pi[1].Z(),np,nv->PartInfo(0)->fP.E()/1000.); }
  fprintf(o,"# sigma %.6e\n",(inv>0?N/inv:0.)*1e-38); fclose(o);
  printf("%s: %ld events, %ld two-pion\n",outfile,N,nsig);
}
