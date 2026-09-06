// sb_transfer.C -- transfer factors T_{b,c} = N_b^SR / N_b^CR per background class b and control
// region c, on the sb_-instrumented MC (CV-weighted, per-run POT-scaled). Classes are built
// from the selection's EventCategory (0 signal, 2 OOFV, 3 other numu CC, 4 nue CC, 5 NC) with the
// dominant class 3 split by true final state from the generator daughter list.
//   root -l -b -q 'macros/sb_transfer.C("fhc")'
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void sb_transfer_pf(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb_pi0/"; std::vector<Src> mc;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); }
  // Framework safe_weight rule: non-finite, negative or >30 -> 1 (one authoritative rule, 2026-09-03)
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* NPI="CC1mu1piXp_mc_n_threshold_pionpm", *NPI0="CC1mu1piXp_mc_n_threshold_pion0";  // selection truth counters
  const int NC=9; TString cls[NC]; const char* cln[NC]={"signal","numuCC 0pi","numuCC 1pi+- out-of-PS","numuCC pi0 (any)","numuCC >=2pi+-, no pi0","numuCC other","OOFV","NC + nueCC","CC+NC with a pi0 (target)"};
  cls[0]="CC1mu1piXp_EventCategory==0";
  cls[1]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+"==0 && "+NPI0+"==0";
  cls[2]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+"==1 && "+NPI0+"==0";
  cls[3]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI0+">=1";
  cls[4]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+">=2 && "+NPI0+"==0";
  cls[5]=TString("CC1mu1piXp_EventCategory==3 && !(")+NPI+"==0 && "+NPI0+"==0) && !("+NPI+"==1 && "+NPI0+"==0) && !("+NPI0+">=1) && !("+NPI+">=2 && "+NPI0+"==0)";
  cls[6]="CC1mu1piXp_EventCategory==2"; cls[7]="(CC1mu1piXp_EventCategory==4||CC1mu1piXp_EventCategory==5)";
  // ratified 2026-09-06: the pi0 region's target class is CC+NC pi0 (overlaps classes 3 and 7; not added to the total)
  cls[8]=TString("(CC1mu1piXp_EventCategory==3||CC1mu1piXp_EventCategory==5) && ")+NPI0+">=1";
  const int NR=5; const char* reg[NR]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0_final","CC1mu1piXp_sb_cosmic"}; const char* rn[NR]={"SR","CC0pi","multipi","pi0","cosmic"};
  double N[NC][NR]={{0}};
  for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str());
    for(int k=0;k<NC;k++) for(int r=0;r<NR;r++){ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",TString(CVW)+"*(("+cls[k]+") && "+reg[r]+")","goff"); N[k][r]+=h.Integral(0,2)*x.scale; } }
  printf("\n==== %s: background composition per region (CV-weighted, POT-scaled) and transfer factors T = N^SR / N^CR\n",mode);
  printf("%-26s",""); for(int r=0;r<NR;r++) printf("%11s",rn[r]); printf("   |"); for(int r=1;r<NR;r++) printf("%9s",rn[r]); printf("\n");
  for(int k=0;k<NC;k++){ printf("%-26s",cln[k]); for(int r=0;r<NR;r++) printf("%11.1f",N[k][r]); printf("   |"); for(int r=1;r<NR;r++) printf("%9.3f",N[k][r]>0?N[k][0]/N[k][r]:0.); printf("\n"); }
  printf("%-26s","total"); for(int r=0;r<NR;r++){ double s=0; for(int k=0;k<NC-1;k++) s+=N[k][r]; printf("%11.1f",s);} printf("\n");
  printf("pi0-region target purity, CC pi0 only: %.1f%%   CC+NC pi0: %.1f%%\n",100*N[3][3]/(N[0][3]+N[1][3]+N[2][3]+N[3][3]+N[4][3]+N[5][3]+N[6][3]+N[7][3]),100*N[8][3]/(N[0][3]+N[1][3]+N[2][3]+N[3][3]+N[4][3]+N[5][3]+N[6][3]+N[7][3]));
}
