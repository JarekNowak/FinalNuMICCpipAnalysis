// mcs_scale_fakedata_1p.C -- data-side MCS momentum-scale variation for the proton-tagged (w/) fake
// data (2026-09-27). Same prescription as mcs_scale_fakedata.C: every selected muon whose momentum is
// the MCS value (reco != range) gets reco -> reco*s in the fake DATA only. The proton-tagged tree also
// carries observables BUILT from the muon momentum (delta p_T, delta alpha_T, delta phi_T, p_n, W_had),
// precomputed by CC1mu1pi1p, so they are recomputed here from the scaled muon exactly as
// CC1mu1pi1p::compute_had_observables does (beam frame about target->vertex, STVTools kOpt1).
// delta phi_T depends on the muon direction only (recomputed values agree to 1e-11) and is copied.
// W_pipr, p_p, p_pi and the angles do not read the muon momentum and are left untouched.
// s=1 is the closure check: the recomputed values must reproduce the stored branches.
//   root -l -b -e 'gSystem->AddIncludePath("-I$PWD/include")' -q 'macros/mcs_scale_fakedata_1p.C+(1.00,"chk")'   // closure only
//   root -l -b -q 'macros/mcs_scale_fakedata_1p.C+(1.05,"up05")'         // FHC runs 1,2,4,5
//   root -l -b -q 'macros/mcs_scale_fakedata_1p.C+(1.05,"up05","rhc")'   // RHC runs 1-4
// Writes w/xsec-ana-fakedata_<beam>_run<i>_mcs<tag>.root next to the nominal files.
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "XSecAnalyzer/NuMIBeamFrame.hh"
#include "../src/utils/STVTools.cxx"   // needs -I include (see usage)

namespace {
  inline double ssqrt( double x ) { return x > 0. ? std::sqrt(x) : 0.; }
}

void mcs_scale_fakedata_1p(double s=1.05, const char* tag="up05", const char* beam="fhc"){
  const char* P="/data/uboone/processed/w/";
  const bool check = std::string(tag)=="chk";
  const std::vector<int> runs = std::string(beam)=="rhc" ? std::vector<int>{1,2,3,4} : std::vector<int>{1,2,4,5};
  const std::string S="CC1mu1pi1p_";
  for(int r : runs){
    TFile* fin=TFile::Open(Form("%sxsec-ana-fakedata_%s_run%d.root",P,beam,r));
    TTree* tin=(TTree*)fin->Get("stv_tree");
    int imu=-1, ipi=-1, ipr=-1; bool sel=false;
    double reco=0, range=0, ppi=0, Whad=0, dat=0, dpht=0, dpt=0, pn=0;
    float vx=0, vy=0, vz=0;
    std::vector<float> *dx=nullptr, *dy=nullptr, *dz=nullptr, *kep=nullptr;
    tin->SetBranchAddress((S+"Selected").c_str(),&sel);
    tin->SetBranchAddress((S+"CandidateMuonIndex").c_str(),&imu);
    tin->SetBranchAddress((S+"CandidatePionIndex").c_str(),&ipi);
    tin->SetBranchAddress((S+"CandidateProtonIndex").c_str(),&ipr);
    tin->SetBranchAddress((S+"candidate_muon_mom_reco").c_str(),&reco);
    tin->SetBranchAddress((S+"candidate_muon_mom_range").c_str(),&range);
    tin->SetBranchAddress((S+"candidate_pion_mom_reco").c_str(),&ppi);
    tin->SetBranchAddress((S+"W_had_reco").c_str(),&Whad);
    tin->SetBranchAddress((S+"deltaAlphaT_reco").c_str(),&dat);
    tin->SetBranchAddress((S+"deltaPhiT_reco").c_str(),&dpht);
    tin->SetBranchAddress((S+"deltaPt_reco").c_str(),&dpt);
    tin->SetBranchAddress((S+"pn_reco").c_str(),&pn);
    tin->SetBranchAddress("reco_nu_vtx_sce_x",&vx); tin->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
    tin->SetBranchAddress("reco_nu_vtx_sce_z",&vz);
    tin->SetBranchAddress("trk_dir_x_v",&dx); tin->SetBranchAddress("trk_dir_y_v",&dy);
    tin->SetBranchAddress("trk_dir_z_v",&dz); tin->SetBranchAddress("trk_energy_proton_v",&kep);

    TFile* out=nullptr; TTree* tout=nullptr;
    double o_reco=0, o_Whad=0, o_dat=0, o_dpht=0, o_dpt=0, o_pn=0;
    if(!check){
      out=new TFile(Form("%sxsec-ana-fakedata_%s_run%d_mcs%s.root",P,beam,r,tag),"recreate"); out->SetCompressionLevel(1);
      tout=tin->CloneTree(0);
      tout->SetBranchAddress((S+"candidate_muon_mom_reco").c_str(),&o_reco);
      tout->SetBranchAddress((S+"W_had_reco").c_str(),&o_Whad);
      tout->SetBranchAddress((S+"deltaAlphaT_reco").c_str(),&o_dat);
      tout->SetBranchAddress((S+"deltaPhiT_reco").c_str(),&o_dpht);
      tout->SetBranchAddress((S+"deltaPt_reco").c_str(),&o_dpt);
      tout->SetBranchAddress((S+"pn_reco").c_str(),&o_pn);
    }
    Long64_t N=tin->GetEntries(); long nmcs=0, ntki=0, nbad=0; double maxdev=0;
    for(Long64_t i=0;i<N;i++){ tin->GetEntry(i);
      o_reco=reco; o_Whad=Whad; o_dat=dat; o_dpht=dpht; o_dpt=dpt; o_pn=pn;
      const bool mcs = reco>0 && reco!=range;
      if(mcs){ o_reco=reco*s; nmcs++; }
      // the TKI block of CC1mu1pi1p::compute_reco_observables runs whenever all three candidates exist
      if(ipr>=0 && imu>=0 && ipi>=0 && (mcs || check)){
        TVector3 mu(dx->at(imu),dy->at(imu),dz->at(imu)); mu.SetMag(o_reco);
        double Emu=ssqrt(o_reco*o_reco+MUON_MASS*MUON_MASS);
        TVector3 pi(dx->at(ipi),dy->at(ipi),dz->at(ipi)); pi.SetMag(ppi);
        double Epi=ssqrt(ppi*ppi+PI_PLUS_MASS*PI_PLUS_MASS);
        double KEp=kep->at(ipr); double pmom=ssqrt(KEp*KEp+2.*PROTON_MASS*KEp);
        TVector3 pr(dx->at(ipr),dy->at(ipr),dz->at(ipr)); pr.SetMag(pmom);
        double Epr=KEp+PROTON_MASS;
        TVector3 nu=NuMIBeam::nu_dir_from_vertex(vx,vy,vz);
        TVector3 mu_b=NuMIBeam::to_beam_frame(mu,nu), had_b=NuMIBeam::to_beam_frame(pi+pr,nu);
        STVTools stv; stv.CalculateSTVs(mu_b,had_b,Emu,Epi+Epr,kOpt1);
        o_dpt=stv.ReturnPt(); o_pn=stv.ReturnPn(); o_dat=stv.ReturnDeltaAlphaT();
        // delta phi_T reads only the muon DIRECTION: keep the stored value so its term is exactly zero
        o_dpht=check ? stv.ReturnDeltaPhiT() : dpht;
        o_Whad=ssqrt(PROTON_MASS*PROTON_MASS+2.*PROTON_MASS*(stv.ReturnECal()-Emu)-stv.ReturnQ2());
        ntki++;
        if(check){
          double d=std::max({std::fabs(o_dpt-dpt),std::fabs(o_pn-pn),std::fabs(o_dat-dat),std::fabs(o_dpht-dpht),std::fabs(o_Whad-Whad)});
          if(d>maxdev) maxdev=d; if(d>1e-6) nbad++;
        }
      }
      if(tout) tout->Fill(); }
    if(check){
      printf("  run%d: %lld entries, %ld recomputed, %ld deviate >1e-6 (max %.3g)\n",r,N,ntki,nbad,maxdev);
    } else {
      tout->Write("",TObject::kOverwrite);
      TParameter<float>* sp=(TParameter<float>*)fin->Get("summed_pot"); if(sp){ out->cd(); sp->Write("summed_pot",TObject::kOverwrite); }
      printf("  run%d: %lld entries, %ld MCS-measured scaled by %.3f, %ld TKI/W_had recomputed -> %s\n",r,N,nmcs,s,ntki,
             Form("xsec-ana-fakedata_%s_run%d_mcs%s.root",beam,r,tag));
      out->Close(); }
    fin->Close(); }
}
