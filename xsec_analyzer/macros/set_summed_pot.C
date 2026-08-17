// set_summed_pot.C — normalise the multi-run FHC numuMC overlays for the
// statistics boost. The framework scales each MC file by data_POT/summed_pot and
// then SUMS the files, so to combine N runs as one sample at Run-1 exposure (not
// N x the exposure) every file's stored summed_pot must be the TOTAL Sigma-POT.
// Then per-file scaling data/Sigma-POT summed over files -> the Run-1 rate, built
// from ~4.5x the events. Backs up the true per-run POTs to pot_backup_run15.txt.
void set_summed_pot(){
  const char* files[]={
    "/data/uboone/processed/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/processed/xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/processed/xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/processed/xsec-ana-Run4c_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/processed/xsec-ana-Run4d_fhc_new_numi_flux_fhc_pandora_ntuple.root"};
  const int N=5;
  double truepot[N]={0}; double sum=0; bool ok=true;
  for(int i=0;i<N;i++){
    TFile f(files[i]);
    if(f.IsZombie()){ printf("ZOMBIE %s\n",files[i]); ok=false; continue; }
    TParameter<float>* p=nullptr; f.GetObject("summed_pot",p);
    truepot[i]= p? p->GetVal() : 0.0;
    if(truepot[i]<=0){ printf("BAD summed_pot in %s\n",files[i]); ok=false; }
    printf("%-70s true_pot=%.6e\n",files[i],truepot[i]);
    sum+=truepot[i]; f.Close();
  }
  if(!ok){ printf("ABORT: a file was missing/bad; not writing.\n"); return; }
  printf("\nSIGMA_POT (Runs 1,2,4,4c,4d) = %.6e\n",sum);
  // backup true POTs
  FILE* bk=fopen("pot_backup_run15.txt","w");
  for(int i=0;i<N;i++) fprintf(bk,"%s %.8e\n",files[i],truepot[i]);
  fclose(bk);
  // overwrite summed_pot = Sigma in each file
  for(int i=0;i<N;i++){
    TFile f(files[i],"UPDATE");
    TParameter<float> np("summed_pot",(float)sum);
    np.Write("summed_pot",TObject::kOverwrite);
    f.Close();
  }
  // verify
  printf("\n--- verify ---\n");
  for(int i=0;i<N;i++){ TFile f(files[i]); TParameter<float>* p=nullptr;
    f.GetObject("summed_pot",p); printf("%-70s summed_pot=%.6e\n",files[i],p?p->GetVal():0); f.Close(); }
  printf("\nDONE. data/SigmaPOT scale = %.6f (was Run1-only %.6f)\n",
         3.283e20/sum, 3.283e20/truepot[0]);
}
