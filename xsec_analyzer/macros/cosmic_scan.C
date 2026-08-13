void cosmic_scan(){
  TChain c("stv_tree");
  c.Add("/data/uboone/processed/sb/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  c.Add("/data/uboone/processed/sb/xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  printf("RES cosmic sideband signal fraction vs opening angle (FHC MC trend)\n");
  double ths[5]={2.6,2.7,2.8,2.9,3.0};
  for(int i=0;i<5;i++){ double th=ths[i];
    TString base=Form("CC1mu1piXp_sb_cosmic && CC1mu1piXp_mu_pi_opening_angle>%f",th);
    double s=c.GetEntries(base+" && CC1mu1piXp_MC_Signal");
    double b=c.GetEntries(base+" && !CC1mu1piXp_MC_Signal");
    printf("RES theta>%.2f  Nsig=%.0f Nbkg=%.0f sig%%=%.1f\n",th,s,b,(s+b)>0?100*s/(s+b):0); fflush(stdout);
  }
  printf("RES DONE\n");
}
