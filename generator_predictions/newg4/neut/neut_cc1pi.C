// neut_cc1pi.C — apply the CC1mu1piXp signal to a NEUT neutvect.root and
// produce per-nucleon dsigma/dx in the analysis binning.
//
// NEUT is run with EVCT-MPV 3 (events ~ flux(E)*sigma(E)), so the signal
// fraction N_signal/N_total = sigma_signal / sigma_total. The flux-averaged
// total cross section sigma_tot_perNuc (per nucleon, 10^-38 cm^2) is passed in
// (read from the NEUT run output). Then dsigma/dx = sigma_tot_perNuc *
// (count_in_bin / N_total) / binwidth.
//
// Requires the NEUT class libraries (neutvect.so etc.) loaded first.
// Run: root -l -b -q 'neut_cc1pi.C("neutvect.root", SIGMA_TOT_PERNUC_1e38, "out.root")'

#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector3.h"
#include "neutvect.h"

static bool is_meson(int pdg){ int a=std::abs(pdg);
  if(a>=9900000) return false; if((a/1000)%10!=0) return false;
  if((a/100)%10==0) return false; if(a>=901&&a<=930) return false;
  if(a==110||a==990||a==998||a==999||a==100) return false; return true; }
static bool is_kaon(int pdg){ int a=std::abs(pdg); return a==321||a==311||pdg==310||pdg==130; }

void neut_cc1pi(const char* infile, double sigma_tot_perNuc, const char* outfile){
  TFile f(infile);
  TTree* t=(TTree*)f.Get("neuttree");
  if(!t){ printf("no neuttree in %s\n",infile); return; }
  NeutVect* nv=new NeutVect();
  t->SetBranchAddress("vectorbranch",&nv);

  std::vector<double> EPMU={0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.30,2.60,3.00};
  std::vector<double> EPPI={0.175,0.20,0.30,0.40,0.50,1.00};
  std::vector<double> ECMU={-1.0,-0.5,0.0,0.2,0.4,0.55,0.65,0.75,0.82,0.88,0.93,0.97,1.0};
  std::vector<double> ECPI={-1.0,-0.7,-0.4,-0.2,0.0,0.2,0.35,0.5,0.65,0.78,0.88,0.95,1.0};
  std::vector<double> ETH={0.0,0.314,0.628,0.942,1.257,1.571,1.885,2.199,2.513,2.6};
  TH1D* hpmu=new TH1D("pmu","",EPMU.size()-1,EPMU.data());
  TH1D* hppi=new TH1D("ppi","",EPPI.size()-1,EPPI.data());
  TH1D* hcmu=new TH1D("costhmu","",ECMU.size()-1,ECMU.data());
  TH1D* hcpi=new TH1D("costhpi","",ECPI.size()-1,ECPI.data());
  TH1D* hth =new TH1D("thmupi","",ETH.size()-1,ETH.data());

  long N=t->GetEntries(); long sig=0;
  for(long i=0;i<N;i++){ t->GetEntry(i);
    // CC only: NEUT Mode |mode|<30 is CC
    if(std::abs(nv->Mode)>=30) continue;
    int nmu=0,npi=0,npi0=0,nk=0,nheavy=0; TVector3 pmu,ppi;
    for(int j=0;j<nv->Npart();j++){
      NeutPart* p=nv->PartInfo(j);
      if(!p->fIsAlive || p->fStatus!=0) continue;   // final-state only
      int pdg=p->fPID, a=std::abs(pdg);
      TVector3 mom(p->fP.X()/1000.,p->fP.Y()/1000.,p->fP.Z()/1000.); // MeV->GeV
      if(a==13){ nmu++; pmu=mom; }
      else if(a==211){ npi++; ppi=mom; }
      else if(pdg==111){ npi0++; }
      else if(is_kaon(pdg)){ nk++; }
      else if(a!=111&&a!=211&&is_meson(pdg)){ nheavy++; }
    }
    if(nmu!=1||npi!=1||npi0!=0||nk!=0||nheavy!=0) continue;
    double Pmu=pmu.Mag(),Ppi=ppi.Mag(),th=pmu.Angle(ppi);
    if(Pmu<=0.15||Ppi<=0.175||th>=2.6) continue;
    ++sig; hpmu->Fill(Pmu);hppi->Fill(Ppi);hcmu->Fill(pmu.CosTheta());hcpi->Fill(ppi.CosTheta());hth->Fill(th);
  }
  auto norm=[&](TH1D* h){ for(int b=1;b<=h->GetNbinsX();++b){ double c=h->GetBinContent(b),w=h->GetBinWidth(b);
    h->SetBinContent(b, sigma_tot_perNuc*(c/(double)N)/w);
    h->SetBinError(b, sigma_tot_perNuc*(std::sqrt(c)/(double)N)/w); } };
  norm(hpmu);norm(hppi);norm(hcmu);norm(hcpi);norm(hth);
  TFile fo(outfile,"recreate"); hpmu->Write();hppi->Write();hcmu->Write();hcpi->Write();hth->Write(); fo.Close();
  printf("NEUT: N=%ld CC1pi-signal=%ld (%.2f%%)  sigma_tot_perNuc=%.3e\n",N,sig,100.*sig/N,sigma_tot_perNuc);
}
