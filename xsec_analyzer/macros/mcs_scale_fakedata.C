// mcs_scale_fakedata.C -- data-side MCS momentum-scale variation (2026-09-19).
// The first MCS bound scaled the MCS momentum in the BIN CONFIGURATION, i.e. in the simulation
// and in the fake data alike; a common scale cancels in the extraction (response and data move
// together), so that test bounded only the regulariser's response to a remapped variable, not the
// systematic. A momentum-scale error is a data-vs-simulation mismatch, so here the fake DATA alone
// are scaled: every selected muon whose momentum is the MCS value (reco != range, 62.9% of the
// sample) gets reco -> reco*s; the simulation, the response and the selection flag are untouched.
// Writes xsec-ana-fakedata_fhc_run<i>_mcs<tag>.root next to the nominal files.
//   root -l -b -q 'macros/mcs_scale_fakedata.C(1.05,"up05")'
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include <cstdio>
void mcs_scale_fakedata(double s=1.05, const char* tag="up05"){
  const char* P="/data/uboone/processed/";
  for(int r : {1,2,4,5}){
    TFile* fin=TFile::Open(Form("%sxsec-ana-fakedata_fhc_run%d.root",P,r));
    TTree* tin=(TTree*)fin->Get("stv_tree");
    double reco=0, range=0; tin->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_reco",&reco);
    tin->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_range",&range);
    TFile* out=new TFile(Form("%sxsec-ana-fakedata_fhc_run%d_mcs%s.root",P,r,tag),"recreate"); out->SetCompressionLevel(1);
    TTree* tout=tin->CloneTree(0); double reco_out=0; tout->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_reco",&reco_out);
    Long64_t N=tin->GetEntries(); long nmcs=0, nsel=0;
    for(Long64_t i=0;i<N;i++){ tin->GetEntry(i); reco_out=reco;
      if(reco>0 && reco!=range){ reco_out=reco*s; nmcs++; }
      if(reco>0) nsel++;
      tout->Fill(); }
    tout->Write("",TObject::kOverwrite);
    TParameter<float>* sp=(TParameter<float>*)fin->Get("summed_pot"); if(sp) sp->Write("summed_pot",TObject::kOverwrite);
    printf("  run%d: %lld entries, %ld with a muon momentum, %ld MCS-measured scaled by %.3f -> %s\n",r,N,nsel,nmcs,s,Form("xsec-ana-fakedata_fhc_run%d_mcs%s.root",r,tag));
    out->Close(); fin->Close(); }
}
