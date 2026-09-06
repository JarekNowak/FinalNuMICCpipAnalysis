// sb_transfer_unc.C -- uncertainty on the transfer factor T = N_b^SR / N_b^CR of each control
// region's target class, computed as a RATIO per universe (flux: weight_ppfx_all x tune x
// norm; xsec: weight_All_UBGenie x ppfx x norm; reint: weight_reint_all x tune x ppfx x
// norm), so the numerator/denominator covariance is included, plus the MC-statistical term
// from sum(w^2) (SR and CR are disjoint, so their statistical terms are independent).
// Also prints the SR nu-MC yield under the two weight rules in use (framework safe_weight
// clamp vs unclamped) to document the ~1% difference between macros.
//   root -l -b -q 'macros/sb_transfer_unc.C("fhc")'
#include <vector>
#include <string>
#include <cmath>
struct Src { std::string file; double scale; };
void sb_transfer_unc_p0(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb_pi0/"; std::vector<Src> mc;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); }
  // (region, target class): CC0pi->0pi, multipi->>=2pi no pi0, pi0->pi0 any, plus the total nu-MC
  const int NT=5; const char* tn[NT]={"CC0pi: numuCC 0pi","multipi: numuCC >=2pi+-","pi0: numuCC pi0 (any)","all: total nu-MC","pi0: CC+NC pi0 (ratified)"};
  const char* creg[NT]={"CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_pi0"};
  const char* SY[3]={"weight_ppfx_all","weight_All_UBGenie","weight_reint_all"}; const char* syn[3]={"flux","xsec","reint"};
  std::vector<double> sSR[3][NT], sCR[3][NT]; double cvSR[NT]={0},cvCR[NT]={0},w2SR[NT]={0},w2CR[NT]={0}; int NU[3]={0,0,0};
  double srClamp=0, srRaw=0;
  auto safe=[](double w){ return (std::isfinite(w)&&w>=0&&w<=30)?w:1.0; };
  for(auto& x:mc){
    TFile f(x.file.c_str()); TTree* t=(TTree*)f.Get("stv_tree"); t->SetBranchStatus("*",0);
    for(auto b:{"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_EventCategory","CC1mu1piXp_mc_n_threshold_pionpm","CC1mu1piXp_mc_n_threshold_pion0","tuned_cv_weight","ppfx_cv_weight","normalisation_weight","weight_ppfx_all","weight_All_UBGenie","weight_reint_all"}) t->SetBranchStatus(b,1);
    bool sel,c0,cm,cp; int cat,npi,npi0; float tune,ppfx,norm; std::vector<double>* wv[3]={nullptr,nullptr,nullptr};
    t->SetBranchAddress("CC1mu1piXp_Selected",&sel); t->SetBranchAddress("CC1mu1piXp_sb_cc0pi",&c0); t->SetBranchAddress("CC1mu1piXp_sb_multipi",&cm); t->SetBranchAddress("CC1mu1piXp_sb_pi0",&cp);
    t->SetBranchAddress("CC1mu1piXp_EventCategory",&cat); t->SetBranchAddress("CC1mu1piXp_mc_n_threshold_pionpm",&npi); t->SetBranchAddress("CC1mu1piXp_mc_n_threshold_pion0",&npi0);
    t->SetBranchAddress("tuned_cv_weight",&tune); t->SetBranchAddress("ppfx_cv_weight",&ppfx); t->SetBranchAddress("normalisation_weight",&norm);
    for(int s=0;s<3;s++) t->SetBranchAddress(SY[s],&wv[s]);
    t->Draw(">>el","CC1mu1piXp_Selected||CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0","entrylist"); TEntryList* el=(TEntryList*)gDirectory->Get("el");
    for(Long64_t k=0;k<el->GetN();k++){ t->GetEntry(el->GetEntry(k));
      bool cls[NT]={ cat==3&&npi==0&&npi0==0, cat==3&&npi>=2&&npi0==0, cat==3&&npi0>=1, true, (cat==3||cat==5)&&npi0>=1 };
      bool inCR[NT]={c0,cm,cp,c0,cp};
      double raw=tune*ppfx*norm, w0=safe(raw);
      if(sel){ srClamp+=w0*x.scale; srRaw+=(std::isfinite(raw)&&raw>=0?raw:0)*x.scale; }
      for(int s=0;s<3;s++){ if(!wv[s]) continue; size_t nu=wv[s]->size(); if(NU[s]==0){NU[s]=nu; for(int j=0;j<NT;j++){sSR[s][j].assign(nu,0.); sCR[s][j].assign(nu,0.);}} }
      for(int j=0;j<NT;j++){ if(!cls[j]) continue;
        if(sel){ cvSR[j]+=w0*x.scale; w2SR[j]+=pow(w0*x.scale,2); } if(inCR[j]){ cvCR[j]+=w0*x.scale; w2CR[j]+=pow(w0*x.scale,2); }
        for(int s=0;s<3;s++){ if(!wv[s]) continue; for(size_t u=0;u<(size_t)NU[s];u++){ double wu=(*wv[s])[u]; double w=(s==0)?wu*tune*norm:(s==1)?wu*ppfx*norm:wu*tune*ppfx*norm; w=safe(w)*x.scale;
          if(sel) sSR[s][j][u]+=w; if(inCR[j]) sCR[s][j][u]+=w; } } }
    }
    delete el;
  }
  printf("\n==== %s: SR nu-MC yield, framework clamp (w>30 -> 1) = %.1f ; unclamped = %.1f\n",mode,srClamp,srRaw);
  printf("==== %s: transfer factor T = N_SR/N_CR and its fractional uncertainty (%%), ratio per universe\n",mode);
  printf("%-26s %8s %8s %8s %7s %7s %7s %7s %8s\n","target class (region)","N_SR","N_CR","T","flux","xsec","reint","MCstat","tot(no det)");
  for(int j=0;j<NT;j++){ double T=cvSR[j]/cvCR[j]; double fr[3]; double tot=0;
    for(int s=0;s<3;s++){ double m=0,v=0; int nu=NU[s]; for(int u=0;u<nu;u++){ double Tu=sSR[s][j][u]/sCR[s][j][u]; v+=pow(Tu-T,2);} fr[s]=sqrt(v/nu)/T; tot+=fr[s]*fr[s]; }
    double st=sqrt(w2SR[j]/(cvSR[j]*cvSR[j])+w2CR[j]/(cvCR[j]*cvCR[j])); tot+=st*st;
    printf("%-26s %8.1f %8.1f %8.4f %6.1f%% %6.1f%% %6.1f%% %6.1f%% %7.1f%%\n",tn[j],cvSR[j],cvCR[j],T,100*fr[0],100*fr[1],100*fr[2],100*st,100*sqrt(tot)); }
}
