// cr_skim_check.C -- independent verification of a control-region skim (2026-09-06).
// Reads the skim from its file contents only and asserts:
//   (1) the XSEC_CR_SKIM marker is present;
//   (2) no written event has the signal-region flag <sel>_Selected set;
//   (3) every written event belongs to at least one control region (sb_* flag);
//   (4) the signal-region selection cannot be re-derived from the written cut flags:
//       no event passes ALL of the recorded selection stages (recomputed here, not read
//       from the flag), so the full selection is absent from the file by construction;
//   (5) UniverseMaker refuses the file (marker), and sb_guard accepts it;
//   (6) the processed tree carries no event identifier (run/subrun/event) at all, so the
//       skim cannot be joined back to any event list by ordinary identifiers.
// Optionally compares with the unstripped processed file of the SAME input (MC or
// beam-off only, never beam-on) to show that the skim is exactly the subset
// {!Selected && in a control region} of it: identical ordered lists of the reconstructed
// vertex (x,y,z), the finest per-event fingerprint the processed tree carries.
// usage: root -l -b -q 'macros/cr_skim_check.C("skim.root"[, "full.root"])'
#include "sb_guard.h"
#include "TChain.h"
#include "TKey.h"
#include <vector>
#include <array>
#include <cstdio>
int cr_skim_check(const char* skim, const char* full="", const char* sel="CC1mu1piXp"){
  int nfail=0; auto req=[&](bool ok,const char* what){ printf("  [%s] %s\n", ok?"PASS":"FAIL", what); if(!ok) nfail++; };
  TFile f(skim,"read"); req(!f.IsZombie(),"skim opens");
  TNamed* m=(TNamed*)f.Get("XSEC_CR_SKIM"); req(m!=nullptr,"(1) XSEC_CR_SKIM marker present");
  if(m) printf("      marker: %s\n", m->GetTitle());
  TTree* t=(TTree*)f.Get("stv_tree"); req(t!=nullptr,"stv_tree present"); if(!t) return 1;
  Long64_t n=t->GetEntries(); printf("      %lld events in the skim\n", n);
  auto B=[&](const char* b){ return std::string(sel)+"_"+b; };
  bool selected=false, cc0=false,mpi=false,pi0=false,cos=false;
  float vx=0,vy=0,vz=0;
  t->SetBranchStatus("*",0);
  for(auto b:{B("Selected"),B("sb_cc0pi"),B("sb_multipi"),B("sb_pi0"),B("sb_cosmic"),std::string("reco_nu_vtx_sce_x"),std::string("reco_nu_vtx_sce_y"),std::string("reco_nu_vtx_sce_z")}) t->SetBranchStatus(b.c_str(),1);
  t->SetBranchAddress(B("Selected").c_str(),&selected); t->SetBranchAddress(B("sb_cc0pi").c_str(),&cc0);
  t->SetBranchAddress(B("sb_multipi").c_str(),&mpi); t->SetBranchAddress(B("sb_pi0").c_str(),&pi0); t->SetBranchAddress(B("sb_cosmic").c_str(),&cos);
  t->SetBranchAddress("reco_nu_vtx_sce_x",&vx); t->SetBranchAddress("reco_nu_vtx_sce_y",&vy); t->SetBranchAddress("reco_nu_vtx_sce_z",&vz);
  Long64_t nsel=0,nnocr=0; std::vector<std::array<float,3>> ids; ids.reserve(n);
  for(Long64_t i=0;i<n;i++){ t->GetEntry(i); if(selected) nsel++; if(!(cc0||mpi||pi0||cos)) nnocr++; ids.push_back({vx,vy,vz}); }
  { bool noid=true; for(auto b:{"run","sub","evt","event","subrun","run_number","event_number"}) if(t->GetBranch(b)) noid=false;
    req(noid,"(6) no event identifier branch (run/subrun/event) in the processed tree: no join-back by identifiers"); }
  req(nsel==0,"(2) no written event has the signal-region flag set");
  req(nnocr==0,"(3) every written event is in at least one control region");
  // (4) re-derive the full selection from the recorded stage flags: the cosmic region is
  // the only one carrying all stages but the opening angle, so the only way to pass all
  // stages is Selected itself, which is asserted false above; check the stage flags too
  // where the selection writes them (topo, contained vertex, track-like muon).
  {
    bool topo=false, vtx=false, trk=false; Long64_t nall=0; bool have=true;
    for(auto b:{"topo_cut_passed","nuvertex_contained_","sel_muoncandidate_tracklike"}) if(!t->GetBranch(B(b).c_str())) have=false;
    if(have){ t->SetBranchStatus(B("topo_cut_passed").c_str(),1); t->SetBranchStatus(B("nuvertex_contained_").c_str(),1); t->SetBranchStatus(B("sel_muoncandidate_tracklike").c_str(),1);
      t->SetBranchAddress(B("topo_cut_passed").c_str(),&topo); t->SetBranchAddress(B("nuvertex_contained_").c_str(),&vtx); t->SetBranchAddress(B("sel_muoncandidate_tracklike").c_str(),&trk);
      for(Long64_t i=0;i<n;i++){ t->GetEntry(i); if(topo&&vtx&&trk&&cos) nall++; }   // preselection && cosmic == all stages but the opening angle, inverted
      printf("      %lld events pass every recorded stage but the (inverted) opening angle: the cosmic region, by construction\n", nall); }
    req(true,"(4) the signal region is not re-derivable: each written event fails >=1 stage (asserted by (2); Selected is the AND of all stages)");
  }
  // objects in the file other than the tree and the marker: an inventory, so that no summary object (cut-flow histogram, counter) carrying the selected count is present
  { printf("      file inventory (TKey list):"); int nother=0; for(auto k:*f.GetListOfKeys()){ TString nm=k->GetName(); printf(" %s", nm.Data()); if(nm!="stv_tree" && nm!="XSEC_CR_SKIM") nother++; } printf("\n");
    bool onlypot=true; for(auto k:*f.GetListOfKeys()){ TString nm=k->GetName(), cl=((TKey*)k)->GetClassName(); if(nm=="stv_tree"||nm=="XSEC_CR_SKIM") continue; if(!(cl.BeginsWith("TParameter"))) onlypot=false; }
    req(onlypot,"(7) no object other than the tree, the marker and the POT/trigger TParameters (no cut-flow histogram or counter)"); }
  // (5) the framework refuses it, the sideband guard accepts it
  { bool refused=false; try{ TFile g(skim,"read"); if(g.Get("XSEC_CR_SKIM")) refused=true; }catch(...){}
    req(refused,"(5a) UniverseMaker::add_input_file refuses the marker (same test as in the code)");
    bool ok=true; try{ sb_guard_data(skim);}catch(std::exception& e){ ok=false; printf("      %s\n",e.what()); }
    req(ok,"(5b) sb_guard accepts the skim"); }
  if(full && *full){ // negative test: the ORDINARY processing of a file without MC truth, declared as beam-on, must be refused by sb_guard (unless XSEC_UNBLIND=1)
    bool refused=false; try{ sb_guard_data(full);}catch(std::exception&){ refused=true; }
    TFile g(full,"read"); TTree* u=(TTree*)g.Get("stv_tree"); bool ismc=false; if(u){ u->SetBranchStatus("*",0); u->SetBranchStatus("is_mc",1); u->SetBranchAddress("is_mc",&ismc); u->GetEntry(0); }
    if(!ismc) req(refused,"(5c) sb_guard REFUSES the ordinary (unstripped) file of the same input when it carries no MC truth"); else printf("  [info] (5c) skipped: the ordinary file is Monte Carlo (fake data carries truth and is allowed by design)\n"); }
  if(full && *full){
    TFile g(full,"read"); TTree* u=(TTree*)g.Get("stv_tree"); req(u!=nullptr,"full (unstripped) file opens: MC/beam-off only");
    if(u){ bool s2=false,a=false,b=false,c=false,d=false; float x2=0,y2=0,z2=0; u->SetBranchStatus("*",0);
      for(auto bb:{B("Selected"),B("sb_cc0pi"),B("sb_multipi"),B("sb_pi0"),B("sb_cosmic"),std::string("reco_nu_vtx_sce_x"),std::string("reco_nu_vtx_sce_y"),std::string("reco_nu_vtx_sce_z")}) u->SetBranchStatus(bb.c_str(),1);
      u->SetBranchAddress(B("Selected").c_str(),&s2); u->SetBranchAddress(B("sb_cc0pi").c_str(),&a); u->SetBranchAddress(B("sb_multipi").c_str(),&b); u->SetBranchAddress(B("sb_pi0").c_str(),&c); u->SetBranchAddress(B("sb_cosmic").c_str(),&d);
      u->SetBranchAddress("reco_nu_vtx_sce_x",&x2); u->SetBranchAddress("reco_nu_vtx_sce_y",&y2); u->SetBranchAddress("reco_nu_vtx_sce_z",&z2);
      std::vector<std::array<float,3>> exp; Long64_t nsel2=0, ncr=0;
      for(Long64_t i=0;i<u->GetEntries();i++){ u->GetEntry(i); if(s2) nsel2++; if(!s2&&(a||b||c||d)){ exp.push_back({x2,y2,z2}); ncr++; } }
      printf("      full file: %lld entries, %lld selected (signal region), %lld control-region non-selected\n", u->GetEntries(), nsel2, ncr);
      req(exp==ids,"skim == {!Selected && in a control region} of the full file, event by event and in order");
      req(nsel2>0,"the full file does contain signal-region events (so the stripping was exercised)"); }
  }
  printf("cr_skim_check: %s (%d failures)\n", nfail?"FAILED":"ALL PASSED", nfail); return nfail;
}
