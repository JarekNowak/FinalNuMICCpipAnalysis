// throw_fakedata_mares.C — build Poisson-THROWN M_A^RES fake data with a
// realistic covariance (unlike the weighted-MC version, whose 2.32M events gave
// a degenerate near-Asimov covariance). Each MC event contributes
// Poisson(mu) weight-1 copies, mu = (tuned_cv*ppfx_cv*norm)*w_MaRES*POTSCALE,
// so each reco bin gets N ~ Poisson(expected count at the real FHC data POT
// 6.626e20). Weight-1 events => the framework's stat covariance is genuine
// Poisson sqrt(N). Signal RES events are reweighted by the axial dipole ratio at
// M_A^RES = 1.12 +/- 0.22 GeV. The bulky systematic-universe maps are dropped
// (fake DATA needs only CV reco + truth). gRandom seeded for reproducibility.
void throw_fakedata_mares(){
  const char* PROC="/data/uboone/processed/";
  const char* runs[]={"Run1_fhc","Run2_fhc","Run4_fhc"}; const int NR=3;
  const double SIGPOT=7.655314e21, DATAPOT=6.626e20, POTSCALE=DATAPOT/SIGPOT;
  const double mmu=0.10566, MA=1.12, MAp=1.34, MAm=0.90;
  const double bx=0.47417, by=0.06980, bz=0.87766;
  gRandom->SetSeed( 20260727 );

  for(int sign=0; sign<2; ++sign){
    double MAv = sign==0 ? MAp : MAm;
    const char* tag = sign==0 ? "plus" : "minus";
    TString outname = Form("%sxsec-ana-fakedata_MaRES_%s_run15.root",PROC,tag);
    printf("=== throwing %s (MA=%.2f) ===\n", outname.Data(), MAv);

    TChain cin("stv_tree");
    for(int r=0;r<NR;r++) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,runs[r]));
    // keep everything EXCEPT the bulky multisim maps (not needed for fake data)
    cin.SetBranchStatus("*",1);
    for(auto b:{"weight_All_UBGenie","weight_ppfx_all","weight_reint_all"}) cin.SetBranchStatus(b,0);

    float tcv,pcv,nw,nuE; int ccnc,inter; bool ismc;
    std::vector<int>*pdg=0; std::vector<float>*E=0,*px=0,*py=0,*pz=0;
    cin.SetBranchAddress("tuned_cv_weight",&tcv);cin.SetBranchAddress("ppfx_cv_weight",&pcv);
    cin.SetBranchAddress("normalisation_weight",&nw);cin.SetBranchAddress("is_mc",&ismc);
    cin.SetBranchAddress("mc_ccnc",&ccnc);cin.SetBranchAddress("mc_interaction",&inter);
    cin.SetBranchAddress("mc_nu_energy",&nuE);
    cin.SetBranchAddress("mc_pdg",&pdg);cin.SetBranchAddress("mc_E",&E);
    cin.SetBranchAddress("mc_px",&px);cin.SetBranchAddress("mc_py",&py);cin.SetBranchAddress("mc_pz",&pz);

    TFile* out=new TFile(outname,"recreate"); out->SetCompressionLevel(1);
    TTree* ot=cin.CloneTree(0);
    float one=1.0f;
    ot->SetBranchAddress("tuned_cv_weight",&one);       // weight-1 pseudo-data
    ot->SetBranchAddress("ppfx_cv_weight",&one);
    ot->SetBranchAddress("normalisation_weight",&one);

    Long64_t N=cin.GetEntries(); long kept=0; double expct=0;
    for(Long64_t i=0;i<N;i++){
      cin.GetEntry(i);
      double cv=tcv*pcv*nw; if(!std::isfinite(cv)||cv<0) continue;
      double w=1.0;
      if(ismc && ccnc==0 && inter==1 && pdg){
        int mi=-1; for(size_t j=0;j<pdg->size();j++) if(abs((*pdg)[j])==13){mi=j;break;}
        if(mi>=0){ double Emu=(*E)[mi], pl=(*px)[mi]*bx+(*py)[mi]*by+(*pz)[mi]*bz;
          double Q2=2.0*nuE*(Emu-pl)-mmu*mmu;
          if(Q2>0){ double b0=1.0+Q2/(MA*MA); w=pow(b0/(1.0+Q2/(MAv*MAv)),4); } }
      }
      double mu = cv*w*POTSCALE; expct += mu;
      int ncopy = gRandom->Poisson(mu);        // Poisson throw
      for(int c=0;c<ncopy;c++){ one=1.0f; ot->Fill(); kept++; }
    }
    printf("  kept %ld pseudo-events (expected %.0f) from %lld MC entries\n",kept,expct,N);
    ot->Write("",TObject::kOverwrite);
    TParameter<float> sp("summed_pot",(float)DATAPOT); sp.Write("summed_pot",TObject::kOverwrite);
    out->Close();
  }
  printf("DONE throw_fakedata_mares\n");
}
