// sb_detvar.C -- detector-variation uncertainty of the neutrino-MC prediction in the signal
// region and each control region, and on the transfer factor of each region's target class,
// from the Run-4 detector-variation samples (which carry the sb_ flags). Framework
// convention: each variation is compared with the detVar CV sample, and the eight
// differences are summed in quadrature (no averaging). Yields are CV-weighted; the common
// POT scale cancels in every ratio, so none is applied.
//   root -l -b -q 'macros/sb_detvar.C("fhc")'   (rhc)
#include <vector>
#include <string>
#include <cmath>
void sb_detvar(const char* mode="fhc"){
  const char* P="/data/uboone/processed/"; std::string tag = std::string(mode)=="fhc" ? "run4fhc" : "run4rhc";
  const char* knobs[8]={"LYdown","LYrayl","Recomb2","SCE","WMAngleXZ","WMAngleYZ","WMX","WMYZ"};
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* NPI="CC1mu1piXp_mc_n_threshold_pionpm", *NPI0="CC1mu1piXp_mc_n_threshold_pion0";
  const int NC=7; TString cls[NC]; const char* cln[NC]={"all nu-MC","signal","numuCC 0pi","numuCC pi0 (any)","numuCC >=2pi+-, no pi0","OOFV+NC+nueCC","CC+NC pi0 (target)"};
  cls[0]="1"; cls[1]="CC1mu1piXp_EventCategory==0";
  cls[2]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+"==0 && "+NPI0+"==0";
  cls[3]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI0+">=1";
  cls[4]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+">=2 && "+NPI0+"==0";
  cls[5]="(CC1mu1piXp_EventCategory==2||CC1mu1piXp_EventCategory==4||CC1mu1piXp_EventCategory==5)";
  cls[6]=TString("(CC1mu1piXp_EventCategory==3||CC1mu1piXp_EventCategory==5) && ")+NPI0+">=1";
  const int NR=5; const char* reg[NR]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic"}; const char* rn[NR]={"SR","CC0pi","multipi","pi0","cosmic"};
  auto yields=[&](const char* file, double N[NC][NR]){ TFile f(file); TTree* t=(TTree*)f.Get("stv_tree");
    for(int k=0;k<NC;k++) for(int r=0;r<NR;r++){ TH1D h("h","",1,-0.5,1.5); t->Draw("0.5>>h",TString(CVW)+"*(("+cls[k]+") && "+reg[r]+")","goff"); N[k][r]=h.Integral(0,2);} };
  double CV[NC][NR]; yields(Form("%sxsec-ana-detvar_%s_CV.root",P,tag.c_str()),CV);
  double sum2[NC][NR]={{0}}, sum2T[NC][NR]={{0}}; double perknob[8][NR];
  for(int kk=0;kk<8;kk++){ double V[NC][NR]; yields(Form("%sxsec-ana-detvar_%s_%s.root",P,tag.c_str(),knobs[kk]),V);
    for(int k=0;k<NC;k++) for(int r=0;r<NR;r++){ double d=(V[k][r]-CV[k][r])/CV[k][r]; sum2[k][r]+=d*d;
      if(r>0 && CV[k][0]>0 && V[k][r]>0){ double T0=CV[k][0]/CV[k][r], T1=V[k][0]/V[k][r]; double dt=(T1-T0)/T0; sum2T[k][r]+=dt*dt; } }
    for(int r=0;r<NR;r++) perknob[kk][r]=(V[0][r]-CV[0][r])/CV[0][r]; }
  printf("\n==== %s (%s detVar set): detector-variation fractional uncertainty (%%), quadrature over 8 knobs\n",mode,tag.c_str());
  printf("%-26s","class / region: yield unc."); for(int r=0;r<NR;r++) printf("%12s",rn[r]); printf("   | T = N_SR/N_CR unc.:"); for(int r=1;r<NR;r++) printf("%9s",rn[r]); printf("\n");
  for(int k=0;k<NC;k++){ printf("%-26s",cln[k]); for(int r=0;r<NR;r++) printf("%11.1f%%",100*sqrt(sum2[k][r])); printf("   |                    "); for(int r=1;r<NR;r++) printf("%8.1f%%",100*sqrt(sum2T[k][r])); printf("\n"); }
  printf("per-knob shift of the total nu-MC yield (%%):\n"); for(int kk=0;kk<8;kk++){ printf("  %-10s",knobs[kk]); for(int r=0;r<NR;r++) printf("%8.1f",100*perknob[kk][r]); printf("\n"); }
  printf("CV yields (unscaled detVar sample):"); for(int r=0;r<NR;r++) printf(" %s=%.0f",rn[r],CV[0][r]); printf("\n");
}
