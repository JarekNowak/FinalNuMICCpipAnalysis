// binning_screen.C -- cheap pre-screen of candidate binnings for every released inclusive observable:
// column-normalised migration diagonal per true bin, efficiency per true bin and the expected selected
// events (signal + background) per reco bin, from the processed ntuples with the per-run POT scales and
// the CV weight under the safe-weight rule. No universes are built; candidates that pass the diagonal
// criterion (>0.68, the BNB criterion the note quotes) go on to univmake + UnfolderNuMI.
//   root -l -b -q 'macros/binning_screen.C("fhc")'   (rhc)
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
struct Cand { std::string obs, label, tvar, rvar; std::vector<double> edges; bool openTop; };
void binning_screen(const char* mode="fhc"){
  const char* P="/data/uboone/processed/"; std::vector<Src> mc;
  auto fhcMC=[&](){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560}; for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); };
  auto rhcMC=[&](){ const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847}; for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]}); for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); };
  if(std::string(mode)=="fhc") fhcMC(); else rhcMC();
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* PPIR="sqrt(pow(sqrt(pow(CC1mu1piXp_candidate_pion_mom_reco,2)+0.011164)-0.10566+0.13957,2)-0.019480)";
  std::vector<Cand> C={
    {"pmu","released 7","CC1mu1piXp_candidate_muon_mom_true","CC1mu1piXp_candidate_muon_mom_reco",{0.15,0.35,0.55,0.75,0.95,1.25,1.75,3.0},true},
    {"pmu","cand 8 (split 1.25-1.75)","CC1mu1piXp_candidate_muon_mom_true","CC1mu1piXp_candidate_muon_mom_reco",{0.15,0.35,0.55,0.75,0.95,1.25,1.5,1.75,3.0},true},
    {"pmu","cand 9","CC1mu1piXp_candidate_muon_mom_true","CC1mu1piXp_candidate_muon_mom_reco",{0.15,0.30,0.45,0.60,0.75,0.95,1.2,1.5,2.0,3.0},true},
    {"costhmu","released 5","CC1mu1piXp_candidate_muon_costh_true","CC1mu1piXp_candidate_muon_costh_reco",{-1,0.45,0.65,0.8,0.9,1},false},
    {"costhmu","cand 6","CC1mu1piXp_candidate_muon_costh_true","CC1mu1piXp_candidate_muon_costh_reco",{-1,0.3,0.55,0.7,0.8,0.9,1},false},
    {"costhmu","cand 7","CC1mu1piXp_candidate_muon_costh_true","CC1mu1piXp_candidate_muon_costh_reco",{-1,0.2,0.45,0.6,0.7,0.8,0.9,1},false},
    {"thetamu","released 5","TMath::ACos(CC1mu1piXp_candidate_muon_costh_true)","TMath::ACos(CC1mu1piXp_candidate_muon_costh_reco)",{0,0.43,0.62,0.84,1.17,3.15},false},
    {"thetamu","cand 7","TMath::ACos(CC1mu1piXp_candidate_muon_costh_true)","TMath::ACos(CC1mu1piXp_candidate_muon_costh_reco)",{0,0.35,0.5,0.65,0.8,1.0,1.3,3.15},false},
    {"costhpi","released 4","CC1mu1piXp_candidate_pion_costh_true","CC1mu1piXp_candidate_pion_costh_reco",{-1,-0.1,0.35,0.75,1},false},
    {"costhpi","cand 5","CC1mu1piXp_candidate_pion_costh_true","CC1mu1piXp_candidate_pion_costh_reco",{-1,-0.3,0.2,0.5,0.75,1},false},
    {"costhpi","cand 6","CC1mu1piXp_candidate_pion_costh_true","CC1mu1piXp_candidate_pion_costh_reco",{-1,-0.4,0,0.3,0.55,0.75,1},false},
    {"thmupi","released 5","CC1mu1piXp_true_mu_pi_opening_angle","CC1mu1piXp_mu_pi_opening_angle",{0,0.6,0.85,1.3,1.85,2.6},false},
    {"thmupi","cand 6","CC1mu1piXp_true_mu_pi_opening_angle","CC1mu1piXp_mu_pi_opening_angle",{0,0.55,0.75,0.95,1.3,1.85,2.6},false},
    {"thmupi","cand 7","CC1mu1piXp_true_mu_pi_opening_angle","CC1mu1piXp_mu_pi_opening_angle",{0,0.5,0.7,0.9,1.15,1.45,1.9,2.6},false},
    {"ppi","released 2","CC1mu1piXp_candidate_pion_mom_true",PPIR,{0.175,0.205,1.0},true},
    {"ppi","cand 3","CC1mu1piXp_candidate_pion_mom_true",PPIR,{0.175,0.205,0.30,1.0},true},
    {"ppi","cand 4","CC1mu1piXp_candidate_pion_mom_true",PPIR,{0.175,0.205,0.26,0.34,1.0},true},
  };
  printf("\n==== %s: binning pre-screen (col-normalised diagonal | efficiency | selected sig+bkg per reco bin)\n",mode);
  for(auto& c:C){ int nb=c.edges.size()-1; std::vector<double> gen(nb,0),sel(nb,0),reco(nb,0),diagn(nb,0);
    for(auto& s:mc){ TFile f(s.file.c_str()); TTree* t=(TTree*)f.Get("stv_tree"); if(!t) continue;
      auto binexpr=[&](const std::string& v,int b,bool tr){ char lo[64],hi[64]; snprintf(lo,64,"%.4f",c.edges[b]); snprintf(hi,64,"%.4f",c.edges[b+1]);
        std::string e="("+v+">="+lo+")"; if(!(c.openTop && b==nb-1 && tr)) e+="&&("+v+"<"+hi+")"; return e; };
      for(int b=0;b<nb;b++){ std::string tb=binexpr(c.tvar,b,true);
        TH1D h1("h1","",1,-0.5,1.5); t->Draw("0.5>>h1",(std::string(CVW)+"*(CC1mu1piXp_MC_Signal&&"+tb+")").c_str(),"goff"); gen[b]+=h1.Integral(0,2)*s.scale;
        TH1D h2("h2","",1,-0.5,1.5); t->Draw("0.5>>h2",(std::string(CVW)+"*(CC1mu1piXp_MC_Signal&&CC1mu1piXp_Selected&&"+tb+")").c_str(),"goff"); sel[b]+=h2.Integral(0,2)*s.scale;
        std::string rb=binexpr(c.rvar,b,false); if(c.openTop && b==nb-1) rb="("+c.rvar+">="+std::to_string(c.edges[b])+")";
        TH1D h3("h3","",1,-0.5,1.5); t->Draw("0.5>>h3",(std::string(CVW)+"*(CC1mu1piXp_MC_Signal&&CC1mu1piXp_Selected&&"+tb+"&&"+rb+")").c_str(),"goff"); diagn[b]+=h3.Integral(0,2)*s.scale;
        TH1D h4("h4","",1,-0.5,1.5); t->Draw("0.5>>h4",(std::string(CVW)+"*(CC1mu1piXp_Selected&&"+rb+")").c_str(),"goff"); reco[b]+=h4.Integral(0,2)*s.scale; } }
    printf("%-8s %-26s",c.obs.c_str(),c.label.c_str()); double mind=9;
    for(int b=0;b<nb;b++){ double d=sel[b]>0?diagn[b]/sel[b]:0; mind=std::min(mind,d); printf(" [%.2f e%.0f%% n%.0f]",d,100*sel[b]/gen[b],reco[b]); }
    printf("  min diag %.2f %s\n",mind,mind>0.68?"PASS":"fail");
  }
}
