// sb_syst.C -- normalisation systematic uncertainty of the neutrino-MC prediction in the signal
// region and in each control region, from the multisim weights (flux: weight_ppfx_all x tune x
// norm; xsec: weight_All_UBGenie x ppfx x norm; reint: weight_reint_all x tune x ppfx x norm),
// i.e. the same weight conventions as UniverseMaker. RMS over universes of the POT-scaled
// weighted yield, divided by the CV yield. Detector variations are NOT included (no sb_-flagged
// detVar samples). Also reports the correlation coefficient between each CR and the SR yield.
//   root -l -b -q 'macros/sb_syst.C("fhc")'
#include <vector>
#include <string>
#include <cmath>
struct Src { std::string file; double scale; };
void sb_syst(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb/"; std::vector<Src> mc;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); }
  const int NR=5; const char* nm[NR]={"SR","CC0pi","multipi","pi0","cosmic"};
  const char* SY[3]={"weight_ppfx_all","weight_All_UBGenie","weight_reint_all"}; const char* syn[3]={"flux","xsec","reint"};
  std::vector<std::vector<double>> sum[3]; double cv[NR]={0}; int NU[3]={0,0,0};
  for(int s=0;s<3;s++) sum[s].assign(NR,std::vector<double>());
  auto safe=[](double w){ return (std::isfinite(w)&&w>=0&&w<=30)?w:1.0; };
  for(auto& x:mc){
    TFile f(x.file.c_str()); TTree* t=(TTree*)f.Get("stv_tree");
    t->SetBranchStatus("*",0);
    for(auto b:{"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic","tuned_cv_weight","ppfx_cv_weight","normalisation_weight","weight_ppfx_all","weight_All_UBGenie","weight_reint_all"}) t->SetBranchStatus(b,1);
    bool fl[NR]; float tune,ppfx,norm; std::vector<double> *wv[3]={nullptr,nullptr,nullptr};
    const char* fb[NR]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic"};
    for(int r=0;r<NR;r++) t->SetBranchAddress(fb[r],&fl[r]);
    t->SetBranchAddress("tuned_cv_weight",&tune); t->SetBranchAddress("ppfx_cv_weight",&ppfx); t->SetBranchAddress("normalisation_weight",&norm);
    for(int s=0;s<3;s++) t->SetBranchAddress(SY[s],&wv[s]);
    // entry list of events in any region
    t->Draw(">>el","CC1mu1piXp_Selected||CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0||CC1mu1piXp_sb_cosmic","entrylist");
    TEntryList* el=(TEntryList*)gDirectory->Get("el");
    for(Long64_t k=0;k<el->GetN();k++){ t->GetEntry(el->GetEntry(k));
      double w0=safe(tune*ppfx*norm);
      for(int r=0;r<NR;r++) if(fl[r]) cv[r]+=w0*x.scale;
      for(int s=0;s<3;s++){ if(!wv[s]) continue; size_t nu=wv[s]->size(); if(NU[s]==0){NU[s]=nu; for(int r=0;r<NR;r++) sum[s][r].assign(nu,0.);}
        for(size_t u=0;u<nu && u<(size_t)NU[s];u++){ double wu=(*wv[s])[u]; double w = (s==0)? wu*tune*norm : (s==1)? wu*ppfx*norm : wu*tune*ppfx*norm; w=safe(w);
          for(int r=0;r<NR;r++) if(fl[r]) sum[s][r][u]+=w*x.scale; } }
    }
    delete el;
  }
  printf("\n==== %s: neutrino-MC normalisation systematics per region (fractional, %%); CV yields\n",mode);
  printf("%-9s %10s","region","CV"); for(int s=0;s<3;s++) printf("%9s",syn[s]); printf("%9s %12s\n","total","corr(CR,SR)");
  for(int r=0;r<NR;r++){ printf("%-9s %10.1f",nm[r],cv[r]); double tot=0; double corr=0;
    for(int s=0;s<3;s++){ double m=0,v=0; int nu=NU[s]; for(int u=0;u<nu;u++) m+=sum[s][r][u]; m/=nu; for(int u=0;u<nu;u++) v+=pow(sum[s][r][u]-cv[r],2); double frac=sqrt(v/nu)/cv[r]; tot+=frac*frac; printf("%8.1f%%",100*frac);
      // correlation with SR for the flux term (dominant), reported once
      if(s==0){ double c=0,va=0,vb=0; for(int u=0;u<nu;u++){ double a=sum[s][r][u]-cv[r], b=sum[s][0][u]-cv[0]; c+=a*b; va+=a*a; vb+=b*b;} corr=c/sqrt(va*vb); } }
    printf("%8.1f%% %12.2f\n",100*sqrt(tot),corr); }
  printf("(flux/xsec/reint multisims; detector variations not included; corr = flux-universe correlation with the SR yield)\n");
}
