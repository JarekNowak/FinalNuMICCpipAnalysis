// Per-run-period data/prediction in the four control regions, beam-on skim.
// Distinguishes a wrong per-run exposure (one run off, others fine) from a
// background-model discrepancy (all runs shift together).
#include "sb_guard.h"
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void sb_perrun_data(){
  const char* P="/data/uboone/processed/sb_pi0/"; /* current selection; sb/ predates the beam-frame fix (totals identical) */ const char* B="/data/uboone/processed/beamon_skim/";
  const char* E="/data/uboone/processed/ext_perrun/xsec-ana-";
  const char* R1="neutrinoselection_filt_run1_beamoff.root", *R3B="neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root","numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5="numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCC=0.98;
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  struct Grp { const char* tag; std::vector<Src> mc, ext; std::vector<std::string> dat; double pot; };
  std::vector<Grp> g;
  // ---- FHC
  g.push_back({"FHC run1",{{std::string(P)+"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",0.14101}},
    {{std::string(E)+R1,OCC*5748692./G1}},{std::string(B)+"xsec-ana-beamon_fhc_run1.root"},3.283});
  {std::vector<Src> e; for(auto f:R4) e.push_back({std::string(E)+f,OCC*4131149./G4});
   g.push_back({"FHC run4",{{std::string(P)+"xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",0.07323}},e,
    {std::string(B)+"xsec-ana-beamon_fhc_run4c.root",std::string(B)+"xsec-ana-beamon_fhc_run4d.root"},2.075});}
  g.push_back({"FHC run5",{{std::string(P)+"xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root",0.11560}},
    {{std::string(E)+R5,OCC*5154196./G5}},{std::string(B)+"xsec-ana-beamon_fhc_run5.root"},2.231});
  // ---- RHC
  g.push_back({"RHC run1",{{std::string(P)+"xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root",0.06728}},
    {{std::string(E)+R1,OCC*1458253./G1}},{std::string(B)+"xsec-ana-beamon_rhc_run1.root"},0.6053});
  {std::vector<Src> m; for(auto s:{"aa","ab","ac","ad","ae"}) m.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066});
   g.push_back({"RHC run3",m,{{std::string(E)+R3B,OCC*10349610./G3}},{std::string(B)+"xsec-ana-beamon_rhc_run3b.root"},5.003});}
  {std::vector<Src> m,e; for(auto s:{"Run4a_rhc","Run4b_rhc","Run4c_rhc"}) m.push_back({std::string(P)+"xsec-ana-"+std::string(s)+"_new_numi_flux_rhc_pandora_ntuple.root",0.08847});
   for(auto f:R4) e.push_back({std::string(E)+f,OCC*6304167./G4});
   g.push_back({"RHC run4",m,e,{std::string(B)+"xsec-ana-beamon_rhc_run4a.root",std::string(B)+"xsec-ana-beamon_rhc_run4b.root"},2.883});}
  const char* reg[4]={"sb_cc0pi","sb_multipi","sb_pi0","sb_cosmic"}; const char* rn[4]={"CC0pi","multipi","pi0","cosmic"};
  auto sum=[&](const std::string& f,const TString& cut,bool w)->double{ TChain c("stv_tree"); c.Add(f.c_str());
    TH1D h("h1","",1,-0.5,1.5); h.SetDirectory(gROOT); c.Draw("0.5>>h1",(w?TString(CVW):TString("1"))+"*("+cut+")","goff"); double v=h.Integral(0,2); return v; };
  printf("\n%-10s %8s | %s\n","period","POT/e20","    CC0pi          multipi           pi0            cosmic      (data / pred)");
  for(auto& G:g){
    printf("%-10s %8.3f |",G.tag,G.pot);
    for(int i=0;i<4;i++){ TString F=TString("CC1mu1piXp_")+reg[i];
      double pr=0,dd=0;
      for(auto&x:G.mc) pr+=x.scale*sum(x.file,F,true);
      for(auto&x:G.ext) pr+=x.scale*sum(x.file,F,false);
      { bool isf=TString(G.tag).BeginsWith("FHC"); double sc_full=isf?0.092402:0.071666, pot_full=isf?8.857:11.082;
        pr+=0.65*sc_full*(G.pot/pot_full)*sum(std::string(P)+"xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",F,true); }
      for(auto&d:G.dat){ sb_guard_data(d); dd+=sum(d,F,false); }
      printf(" %7.0f/%7.0f=%5.3f",dd,pr,pr>0?dd/pr:0.);
    }
    printf("\n");
  }
}
