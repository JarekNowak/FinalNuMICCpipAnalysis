// mcs_map_fakedata.C -- inclusive fake data with the MCS-measured muon momenta (reco != range) mapped by the
// control-region mapping of scripts/mcs_map_derive.py, reco -> f(reco), in the fake DATA only (impact study of the MCS
// mis-modelling of exiting muons, 2026-10-09). Copy of mcs_scale_fakedata.C (version 1.8) with reco*s replaced by
// f(reco); output in /data/uboone/processed/mcsmap_study/fakedata/, not next to the release files.
//   root -l -b -q 'macros/mcs_map_fakedata.C("fhc")'      root -l -b -q 'macros/mcs_map_fakedata.C("rhc")'
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
namespace {
  std::vector<double> gp, gf;
  void load_map( const char* f ) {
    std::ifstream in( f ); std::string l;
    while ( std::getline( in, l ) ) { if ( l.empty() || l[0] == '#' ) continue; std::istringstream ss( l ); double a, b; ss >> a >> b; gp.push_back( a ); gf.push_back( b ); }
  }
  double fmap( double p ) {
    if ( p <= gp.front() ) return gf.front();
    if ( p >= gp.back() ) return p * gf.back() / gp.back();
    size_t i = std::upper_bound( gp.begin(), gp.end(), p ) - gp.begin();
    const double t = ( p - gp[i-1] ) / ( gp[i] - gp[i-1] );
    return gf[i-1] + t * ( gf[i] - gf[i-1] );
  }
}
void mcs_map_fakedata(const char* beam="fhc", const char* mapfile="/data/uboone/processed/mcsmap_study/mcs_map.txt"){
  const char* P="/data/uboone/processed/";
  const char* O="/data/uboone/processed/mcsmap_study/fakedata/";
  load_map( mapfile );
  const std::vector<int> runs = std::string(beam)=="rhc" ? std::vector<int>{1,2,3,4} : std::vector<int>{1,2,4,5};
  for(int r : runs){
    TFile* fin=TFile::Open(Form("%sxsec-ana-fakedata_%s_run%d.root",P,beam,r));
    TTree* tin=(TTree*)fin->Get("stv_tree");
    double reco=0, range=0; tin->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_reco",&reco);
    tin->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_range",&range);
    TFile* out=new TFile(Form("%sxsec-ana-fakedata_%s_run%d.root",O,beam,r),"recreate"); out->SetCompressionLevel(1);
    TTree* tout=tin->CloneTree(0); double reco_out=0; tout->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_reco",&reco_out);
    Long64_t N=tin->GetEntries(); long nmcs=0, nsel=0;
    for(Long64_t i=0;i<N;i++){ tin->GetEntry(i); reco_out=reco;
      if(reco>0 && reco!=range){ reco_out=fmap(reco); nmcs++; }
      if(reco>0) nsel++;
      tout->Fill(); }
    tout->Write("",TObject::kOverwrite);
    TParameter<float>* sp=(TParameter<float>*)fin->Get("summed_pot"); if(sp) sp->Write("summed_pot",TObject::kOverwrite);
    printf("  run%d: %lld entries, %ld with a muon momentum, %ld MCS-measured mapped -> %s\n",r,N,nsel,nmcs,Form("%sxsec-ana-fakedata_%s_run%d.root",O,beam,r));
    out->Close(); fin->Close(); }
}
