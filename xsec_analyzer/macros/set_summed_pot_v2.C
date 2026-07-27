// Corrected MC normalisation: use {Run1, Run2, Run4}; Run4 = Run4c+Run4d (same
// events) so 4c/4d are dropped to avoid double-counting. Sigma-POT = 7.655e21.
// Restore 4c/4d to their true POT (now unused). True per-run POTs from the
// original processing (backed up in pot_backup_run15.txt).
void set_summed_pot_v2(){
  const double p1=2.328199e21, p2=2.493374e21, p4=2.833741e21;
  const double SIG=p1+p2+p4;                 // 7.655314e21
  printf("SIGMA(Run1,Run2,Run4) = %.6e\n", SIG);
  struct{const char* r; double pot;} set[]={
    {"Run1_fhc",SIG},{"Run2_fhc",SIG},{"Run4_fhc",SIG},          // combined MC POT
    {"Run4c_fhc",1.104199e21},{"Run4d_fhc",1.729525e21}};        // restore true (unused)
  for(auto s:set){
    TString p=Form("/data/uboone/processed/xsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",s.r);
    TFile f(p,"UPDATE"); TParameter<float> sp("summed_pot",(float)s.pot);
    sp.Write("summed_pot",TObject::kOverwrite); f.Close();
    printf("  %-10s summed_pot=%.6e\n",s.r,s.pot);
  }
  printf("data/SigmaPOT scale = %.6f (onBNB 8.817e20)\n", 8.817e20/SIG);
}
