// neut_cc1pi.C — apply the CC1mu1piXp signal to a NEUT neutvect.root and
// produce per-nucleon dsigma/dx in the analysis binning.
//
// NEUT is run with EVCT-MPV 3 (events ~ flux(E)*sigma_tot(E)), so:
//   - signal fraction N_signal/N_all = sigma_signal_fluxavg / sigma_tot_fluxavg
//   - sigma_tot_fluxavg = INT(sigma phi)/INT(phi) = 1/<1/Totcrs> (harmonic mean
//     of the per-event total cross section over ALL events; derived from the
//     flux*sigma sampling pdf).
// Totcrs is per NUCLEON (in 1e-38 cm^2): verified directly from neutvect_numu.root
// -- per-event Totcrs spans 0.001-9.67 (mean 1.71, harmonic mean 0.491), which is
// a per-nucleon magnitude at these energies (per-Ar would be ~20-40). So the
// flux-averaged harmonic mean IS already per-nucleon and must NOT be divided by A;
// combine_newg4.C multiplies by A=40 to reach the per-nucleus dsigma/dx.
//
// Requires the NEUT class libraries (neutvect.so etc.) loaded first, and ACLiC
// compilation (root 'neut_cc1pi.C+') since CINT (ROOT 5.34) can't handle this.
// Run: root -l -b -q 'neut_cc1pi.C+("neutvect.root","out.root")'

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

void neut_cc1pi(const char* infile, const char* outfile){
  TFile f(infile);
  TTree* t=(TTree*)f.Get("neuttree");
  if(!t){ printf("no neuttree in %s\n",infile); return; }
  NeutVect* nv=new NeutVect();
  t->SetBranchAddress("vectorbranch",&nv);

  std::vector<double> EPMU={0.150,0.350,0.550,0.750,0.950,1.250,1.750,3.000};
  std::vector<double> EPPI={0.113,0.175,0.250,0.320,0.420,0.550,1.000};
  std::vector<double> ECMU={-1.0,0.45,0.65,0.80,0.90,1.0};
  std::vector<double> ECPI={-1.0,-0.10,0.35,0.55,0.75,1.0};
  std::vector<double> ETH={0.0,0.60,0.85,1.10,1.30,1.52,1.85,2.60};
  TH1D* hpmu=new TH1D("pmu","",EPMU.size()-1,EPMU.data());
  TH1D* hppi=new TH1D("ppi","",EPPI.size()-1,EPPI.data());
  TH1D* hcmu=new TH1D("costhmu","",ECMU.size()-1,ECMU.data());
  TH1D* hcpi=new TH1D("costhpi","",ECPI.size()-1,ECPI.data());
  TH1D* hth =new TH1D("thmupi","",ETH.size()-1,ETH.data());

  long N=t->GetEntries(); long sig=0;
  double sum_inv_totcrs=0.;   // for the harmonic-mean flux-avg total xsec
  for(long i=0;i<N;i++){ t->GetEntry(i);
    if(nv->Totcrs>0.) sum_inv_totcrs += 1.0/nv->Totcrs;
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
    if(Pmu<=0.15||Ppi<=0.113||th>=2.6) continue;      // ppi relaxed to 0.113 (low bin)
    ++sig; hppi->Fill(Ppi);                            // ppi extends to 0.113
    if(Ppi>0.175){ hpmu->Fill(Pmu);hcmu->Fill(pmu.CosTheta());hcpi->Fill(ppi.CosTheta());hth->Fill(th); } // others: standard 0.175
  }
  // flux-averaged total cross section = harmonic mean of Totcrs (per Ar, 1e-38
  // cm^2). Per-nucleon in cm^2 for combine_newg4.C.
  double sigma_tot_perNuc_1e38 = (sum_inv_totcrs>0.) ? (double)N/sum_inv_totcrs : 0.;
  double sigma_tot_perNuc = sigma_tot_perNuc_1e38 * 1e-38; // Totcrs already per-nucleon
  auto norm=[&](TH1D* h){ for(int b=1;b<=h->GetNbinsX();++b){ double c=h->GetBinContent(b),w=h->GetBinWidth(b);
    h->SetBinContent(b, sigma_tot_perNuc*(c/(double)N)/w);
    h->SetBinError(b, sigma_tot_perNuc*(std::sqrt(c)/(double)N)/w); } };
  norm(hpmu);norm(hppi);norm(hcmu);norm(hcpi);norm(hth);
  TFile fo(outfile,"recreate"); hpmu->Write();hppi->Write();hcmu->Write();hcpi->Write();hth->Write(); fo.Close();
  printf("NEUT: N=%ld CC1pi-signal=%ld (%.2f%%)  sigma_tot_fluxavg=%.3e [1e-38/nucleon] = %.3e [1e-38/Ar]\n",
    N,sig,100.*sig/N,sigma_tot_perNuc_1e38,sigma_tot_perNuc_1e38*40.0);
}
