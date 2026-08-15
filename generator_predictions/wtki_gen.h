// wtki_gen.h — true-level W/TKI observables for the proton-tagged CC1mu1pi1p
// generator predictions, from the muon, charged-pion and (leading) proton
// three-momenta [GeV]. Inlines the same formulas as STVTools::CalculateSTVs
// (option 1, argon removal energy 24.78 MeV) so the generator predictions match the
// selection's true observables, with no framework link dependency. z = beam axis.
#pragma once
#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>
#include <string>
#include <vector>
namespace wtki {
  const double MMU=0.10565837, MPI=0.13957000, MP=0.93827208, MN=0.93956541, EB=0.02478;
  struct Obs { double Wpipr,Whad,dpt,dalphat,dphit,pn; };
  inline double ssq(double x){ return x>0.?std::sqrt(x):0.; }

  inline Obs compute(const TVector3& mu,const TVector3& pi,const TVector3& pr){
    double Emu=std::sqrt(mu.Mag2()+MMU*MMU), Epi=std::sqrt(pi.Mag2()+MPI*MPI), Epr=std::sqrt(pr.Mag2()+MP*MP);
    Obs o;
    TLorentzVector p4pi(pi,Epi), p4pr(pr,Epr);
    o.Wpipr=(p4pi+p4pr).M();
    TVector3 had=pi+pr; double Ehad=Epi+Epr;
    TVector3 muT(mu.X(),mu.Y(),0.), hadT(had.X(),had.Y(),0.);
    TVector3 ptv=muT+hadT; double pt=ptv.Mag();
    o.dpt=pt;
    double da=std::acos((-muT.Dot(ptv))/(muT.Mag()*pt))*180./M_PI; if(da>180)da-=180; if(da<0)da+=180; o.dalphat=da;
    double dp=std::acos((-muT.Dot(hadT))/(muT.Mag()*hadT.Mag()))*180./M_PI; if(dp>180)dp-=180; if(dp<0)dp+=180; o.dphit=dp;
    double ProtonKE=Ehad-MP;               // hadronic KE (STV tool treats the system with m_p)
    double Ecal=Emu+ProtonKE+EB;
    TLorentzVector nu(0.,0.,Ecal,Ecal), mu4(mu,Emu);
    double Q2=-(nu-mu4).Mag2();
    o.Whad=ssq(MP*MP+2.*MP*(Ecal-Emu)-Q2);
    double MA=22.*MN+18.*MP-0.34381, MAP=MA-MN+EB;
    double R=MA+mu.Z()+had.Z()-Emu-Ehad;
    double pL=0.5*R-(MAP*MAP+pt*pt)/(2.*R);
    o.pn=std::sqrt(pt*pt+pL*pL);
    return o;
  }
  // 3-bin analysis edges (match ccpi1p_*_bin_config.txt; coarsened for the low-stat
  // proton-tagged unfolding; first/last bins open -> clamp)
  // Six equal-population bins, matching the analysis bin configs
  // (xsec_analyzer/configs/ccpi1p_<obs>_bin_config.txt). The interior edges are copied
  // from those files; the first/last values are the open outer limits (the analysis
  // first and last bins are open, and clamp() below folds overflow into them).
  // Previously 3 bins -- the generator predictions were then silently dropped by the
  // bin-count match in dsigma_ccpi1p.C once the measurement moved to 6 bins.
  inline std::vector<double> edges(const std::string& k){
    if(k=="Wpipr")  return {1.08,1.19,1.23,1.27,1.34,1.47,2.00};
    if(k=="Whad")   return {0.0,0.78,1.08,1.21,1.31,1.44,1.90};
    if(k=="dpt")    return {0.0,0.43,0.59,0.72,0.88,1.17,1.60};
    if(k=="dalphat")return {0.,80.,120.,142.,158.,170.,180.};
    if(k=="dphit")  return {0.,37.,63.,84.,107.,139.,180.};
    if(k=="pn")     return {0.0,0.52,0.67,0.81,0.98,1.25,1.70};
    return {};
  }
  inline double clamp(double v,const std::vector<double>& e){
    double lo=e.front()+1e-6, hi=e.back()-1e-6; return v<lo?lo:(v>hi?hi:v);
  }
}
