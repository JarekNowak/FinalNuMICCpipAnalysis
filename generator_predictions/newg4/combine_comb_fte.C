// Build the COMBINED (FHC+RHC) generator FTE prediction from the per-mode FTE
// files. The flux-averaged xsec over the joint exposure is the fluence-weighted
// average of the two per-mode flux-averaged xsecs:
//   sigma_comb = (E_FHC*sigma_FHC + E_RHC*sigma_RHC)/(E_FHC+E_RHC)
// with per-mode fluence  E_mode = Phi_tot,mode * POT_mode.  Since the FTE bin
// content is sigma_bin*width, the SAME weighted average applies bin-by-bin to
// the FTE histograms directly (widths identical between modes).
//
//   FHC: Phi_numu=4.43515e-10, Phi_numubar=2.37644e-10  (>60 MeV, /POT/cm2)
//   RHC: Phi_numu=2.92348e-10, Phi_numubar=3.52298e-10
//   POT_FHC = 8.857e20  (nue-matched FHC exposure)
//   POT_RHC = 1.1082e21 (full RHC exposure)
void combine_comb_fte(const char* fhc_fte,const char* rhc_fte,const char* out_fte){
  const double PhiF = 4.43515e-10 + 2.37644e-10;   // 6.81159e-10
  const double PhiR = 2.92348e-10 + 3.52298e-10;   // 6.44646e-10
  const double POTF = 8.857e20, POTR = 1.1082e21;
  const double wF = PhiF*POTF, wR = PhiR*POTR, WT = wF+wR;
  printf("  weights: E_FHC=%.4e E_RHC=%.4e  frac(FHC)=%.4f frac(RHC)=%.4f\n",
         wF,wR,wF/WT,wR/WT);
  // NOTE: thetamu was added to the per-mode FTE files later than the rest and was
  // missing from this list, so the combined file had no thetamu_fte and the combined
  // theta_mu unfold could not run at all. Keep this list in step with make_fte.C.
  const char* obs[6]={"pmu","ppi","costhmu","costhpi","thmupi","thetamu"};
  TFile fF(fhc_fte), fR(rhc_fte), fo(out_fte,"recreate");
  for(int o=0;o<6;o++){
    TString nm=Form("%s_fte",obs[o]);
    TH1D* a=(TH1D*)fF.Get(nm); TH1D* b=(TH1D*)fR.Get(nm);
    if(!a||!b){printf("  missing %s (FHC:%p RHC:%p)\n",nm.Data(),(void*)a,(void*)b);continue;}
    TH1D* h=(TH1D*)a->Clone(nm); h->Reset(); h->SetDirectory(&fo);
    for(int k=1;k<=a->GetNbinsX();++k){
      h->SetBinContent(k,(wF*a->GetBinContent(k)+wR*b->GetBinContent(k))/WT);
      h->SetBinError(k, sqrt(pow(wF*a->GetBinError(k),2)+pow(wR*b->GetBinError(k),2))/WT);
    }
    h->Write();
  }
  fo.Close(); printf("wrote %s\n",out_fte);
}
