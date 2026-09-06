// sb_detvar_split.C -- MONTE-CARLO-STATISTICAL FLOOR of the detector-variation term: the same
// quantities as sb_detvar.C evaluated on two disjoint halves of the generated events (split on the
// true neutrino energy, identical across the CV and the eight variation samples). For each knob the
// difference of the two half-sample fractional shifts, divided by two, estimates the statistical
// noise of the full-sample shift; its quadrature over knobs is the floor to compare with the term.
// Original: sb_detvar.C -- detector-variation uncertainty of the neutrino-MC prediction in the signal
// region and each control region, and on the transfer factor of each region's target class,
// from the Run-4 detector-variation samples (which carry the sb_ flags). Framework
// convention: each variation is compared with the detVar CV sample, and the eight
// differences are summed in quadrature (no averaging). Yields are CV-weighted; the common
// POT scale cancels in every ratio, so none is applied.
//   root -l -b -q 'macros/sb_detvar.C("fhc")'   (rhc)
#include <vector>
#include <string>
#include <cmath>
void sb_detvar_split_pf(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb_pi0/"; std::string tag = std::string(mode)=="fhc" ? "run4fhc" : "run4rhc";
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
  const int NR=5; const char* reg[NR]={"CC1mu1piXp_Selected","CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0_final","CC1mu1piXp_sb_cosmic"}; const char* rn[NR]={"SR","CC0pi","multipi","pi0","cosmic"};
  auto yields=[&](const char* file, double N[NC][NR], int half){ TFile f(file); TTree* t=(TTree*)f.Get("stv_tree"); TString hs = half<0 ? "1" : Form("(int(mc_nu_energy*1.e6)%%2)==%d",half);
    for(int k=0;k<NC;k++) for(int r=0;r<NR;r++){ TH1D h("h","",1,-0.5,1.5); t->Draw("0.5>>h",TString(CVW)+"*(("+cls[k]+") && "+reg[r]+" && "+hs+")","goff"); N[k][r]=h.Integral(0,2);} };
  double CVh[2][NC][NR]; yields(Form("%sxsec-ana-detvar_%s_CV.root",P,tag.c_str()),CVh[0],0); yields(Form("%sxsec-ana-detvar_%s_CV.root",P,tag.c_str()),CVh[1],1);
  double s2[NC][NR]={{0}}, f2[NC][NR]={{0}}, s2T[NC][NR]={{0}}, f2T[NC][NR]={{0}};
  for(int kk=0;kk<8;kk++){ double Vh[2][NC][NR]; for(int h=0;h<2;h++) yields(Form("%sxsec-ana-detvar_%s_%s.root",P,tag.c_str(),knobs[kk]),Vh[h],h);
    for(int k=0;k<NC;k++) for(int r=0;r<NR;r++){ double d[2],dt[2]; for(int h=0;h<2;h++){ d[h]=(Vh[h][k][r]-CVh[h][k][r])/CVh[h][k][r];
        dt[h]=0; if(r>0&&CVh[h][k][0]>0&&Vh[h][k][r]>0){ double T0=CVh[h][k][0]/CVh[h][k][r], T1=Vh[h][k][0]/Vh[h][k][r]; dt[h]=(T1-T0)/T0; } }
      double dfull=0.5*(d[0]+d[1]), dnoise=0.5*(d[0]-d[1]); s2[k][r]+=dfull*dfull; f2[k][r]+=dnoise*dnoise;
      double tfull=0.5*(dt[0]+dt[1]), tnoise=0.5*(dt[0]-dt[1]); s2T[k][r]+=tfull*tfull; f2T[k][r]+=tnoise*tnoise; } }
  printf("\n==== %s (%s detVar set): detector term (quadrature over 8 knobs, mean of halves) and its MC-statistical floor (%%)\n",mode,tag.c_str());
  printf("%-26s","class / region: yield");  for(int r=0;r<NR;r++) printf("%16s",rn[r]); printf("   | T=N_SR/N_CR:"); for(int r=1;r<NR;r++) printf("%16s",rn[r]); printf("\n");
  for(int k=0;k<NC;k++){ printf("%-26s",cln[k]); for(int r=0;r<NR;r++) printf("%8.1f/%6.1f%%",100*sqrt(s2[k][r]),100*sqrt(f2[k][r])); printf("   |             "); for(int r=1;r<NR;r++) printf("%8.1f/%6.1f%%",100*sqrt(s2T[k][r]),100*sqrt(f2T[k][r])); printf("\n"); }
  printf("(term/floor; a term close to its floor is dominated by the finite statistics of the variation samples)\n");
  printf("CV half-sample yields:"); for(int r=0;r<NR;r++) printf(" %s=%.0f/%.0f",rn[r],CVh[0][0][r],CVh[1][0][r]); printf("\n");
}
