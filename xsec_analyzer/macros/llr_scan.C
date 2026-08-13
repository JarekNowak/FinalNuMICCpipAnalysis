void llr_scan(){
  TFile f("/data/uboone/processed/wtest/xsec-ana-Run1_fhc_llr.root");
  TTree*t=(TTree*)f.Get("stv_tree");
  if(!t){printf("RESULT no tree\n");return;}
  double Ngen=t->GetEntries("CC1mu1pi1p_MC_Signal");
  printf("RESULT Ngen=%.0f\n",Ngen); fflush(stdout);
  double cuts[6]={0.2,0.1,0.05,0.0,-0.1,-0.2};
  printf("RESULT %-7s %6s %6s %7s %7s %8s\n","LLR<","Nsig","Nbkg","eff%","pur%","S/sqrtB"); fflush(stdout);
  for(int i=0;i<6;i++){ double x=cuts[i];
    TString base=Form("CC1mu1pi1p_Selected && CC1mu1pi1p_proton_llr_reco<%f",x);
    double ns=t->GetEntries(base+" && CC1mu1pi1p_MC_Signal");
    double nb=t->GetEntries(base+" && !CC1mu1pi1p_MC_Signal");
    printf("RESULT %-7.2f %6.0f %6.0f %7.1f %7.1f %8.1f\n",x,ns,nb,100*ns/Ngen,100*ns/(ns+nb),nb>0?ns/sqrt(nb):0); fflush(stdout);
  }
  for(int k=0;k<2;k++){ double x=k?0.0:0.2; TH1D h("h","",100,-1,1);
    t->Draw("(CC1mu1pi1p_W_pipr_reco-CC1mu1pi1p_W_pipr_true)/CC1mu1pi1p_W_pipr_true>>h",
      Form("CC1mu1pi1p_Selected && CC1mu1pi1p_MC_Signal && CC1mu1pi1p_W_pipr_true>0 && CC1mu1pi1p_proton_llr_reco<%f",x),"goff");
    printf("RESULT Wpipr_res(LLR<%.1f) mean=%.3f RMS=%.3f\n",x,h.GetMean(),h.GetRMS()); fflush(stdout); }
  printf("RESULT DONE\n");
}
