// neut_cc1pi1p.C — proton-tagged (CC1mu1pi1p) W/TKI prediction from a NEUT
// neutvect.root. Inclusive signal (neut_cc1pi.C) plus a leading final-state proton
// above 0.3 GeV/c; six W/TKI observables (wtki_gen.h), per-nucleon dsigma/dx (absolute
// cm^2; combine_wtki.C -> per-Ar 1e-38). Needs the NEUT class libs + ACLiC.
//   root -l -b -q 'neut_cc1pi1p.C+("neutvect_numu.root","pred_1p_numu.root")'
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector3.h"
#include "neutvect.h"
#include "../../wtki_gen.h"

static bool is_meson(int pdg){ int a=std::abs(pdg);
  if(a>=9900000)return false; if((a/1000)%10!=0)return false; if((a/100)%10==0)return false;
  if(a>=901&&a<=930)return false; if(a==110||a==990||a==998||a==999||a==100)return false; return true; }
static bool is_kaon(int pdg){ int a=std::abs(pdg); return a==321||a==311||pdg==310||pdg==130; }

void neut_cc1pi1p(const char* infile, const char* outfile){
  TFile f(infile); TTree* t=(TTree*)f.Get("neuttree");
  if(!t){ printf("no neuttree in %s\n",infile); return; }
  NeutVect* nv=new NeutVect(); t->SetBranchAddress("vectorbranch",&nv);
  const char* OBS[6]={"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  TH1D* h[6]; for(int o=0;o<6;o++){ std::vector<double> ed=wtki::edges(OBS[o]); h[o]=new TH1D(OBS[o],OBS[o],(int)ed.size()-1,ed.data()); }
  long N=t->GetEntries(); long sig=0; double sum_inv_totcrs=0.;
  for(long i=0;i<N;i++){ t->GetEntry(i);
    if(nv->Totcrs>0.) sum_inv_totcrs+=1.0/nv->Totcrs;
    if(std::abs(nv->Mode)>=30) continue;                       // CC only
    int nmu=0,npi=0,npi0=0,nk=0,nheavy=0; TVector3 pmu,ppi,ppr; double leadp=0.;
    for(int j=0;j<nv->Npart();j++){ NeutPart* p=nv->PartInfo(j);
      if(!p->fIsAlive||p->fStatus!=0) continue;                // final state
      int pdg=p->fPID,a=std::abs(pdg);
      TVector3 mom(p->fP.X()/1000.,p->fP.Y()/1000.,p->fP.Z()/1000.);
      if(a==13){ nmu++; pmu=mom; }
      else if(a==211){ npi++; ppi=mom; }
      else if(pdg==111){ npi0++; }
      else if(pdg==2212){ if(mom.Mag()>leadp){ leadp=mom.Mag(); ppr=mom; } }
      else if(is_kaon(pdg)){ nk++; }
      else if(a!=111&&a!=211&&is_meson(pdg)){ nheavy++; } }
    if(nmu!=1||npi!=1||npi0!=0||nk!=0||nheavy!=0) continue;
    if(leadp<=0.3) continue;
    double Pmu=pmu.Mag(),Ppi=ppi.Mag(),th=pmu.Angle(ppi);
    if(Pmu<=0.15||Ppi<=0.175||th>=2.6) continue;
    ++sig;
    wtki::Obs w=wtki::compute(pmu,ppi,ppr);
    double v[6]={w.Wpipr,w.Whad,w.dpt,w.dalphat,w.dphit,w.pn};
    for(int o=0;o<6;o++){ std::vector<double> ed=wtki::edges(OBS[o]); h[o]->Fill(wtki::clamp(v[o],ed)); }
  }
  double sperNuc=((sum_inv_totcrs>0.)?(double)N/sum_inv_totcrs:0.)*1e-38;  // per-nucleon, cm^2
  for(int o=0;o<6;o++){ for(int b=1;b<=h[o]->GetNbinsX();b++){ double c=h[o]->GetBinContent(b),w=h[o]->GetBinWidth(b);
    h[o]->SetBinContent(b,sperNuc*(c/(double)N)/w); h[o]->SetBinError(b,sperNuc*(std::sqrt(c)/(double)N)/w); } }
  TFile fo(outfile,"recreate"); for(int o=0;o<6;o++) h[o]->Write(); fo.Close();
  printf("NEUT 1p: N=%ld  proton-tagged W/TKI signal=%ld (%.3f%%)  wrote %s\n",N,sig,100.*sig/(double)N,outfile);
}
