// build_fakedata_mares.C — build two high-statistics fake-data ntuples from the
// combined 5-run FHC GENIE overlay, with the CC-RES cross section reweighted by
// M_A^RES = 1.12 +/- 0.22 GeV (+/-1 sigma). The reweight is the axial dipole
// |F_A|^2 ratio, w(Q2) = [(1+Q2/MA^2)/(1+Q2/MA'^2)]^4, applied to true CCRES
// events (mc_ccnc==0 && mc_interaction==1). Q2 is the true interaction momentum
// transfer, Q2 = 2 E_nu (E_mu - p_mu . nhat) - m_mu^2, with nhat the mean NuMI
// beam direction and the true muon taken from the mc_ particle arrays.
//
// Implementation: fast-clone every branch EXCEPT tuned_cv_weight (bulk copy, no
// decompression of the big multisim maps), then add a fresh tuned_cv_weight =
// orig * w+/- by looping over only the ~10 branches needed for Q2. summed_pot is
// set to Sigma-POT so the fake data represents the full combined FHC exposure.
void build_fakedata_mares(){
  const char* PROC="/data/uboone/processed/";
  // Corrected set: Run4 = Run4c+Run4d (same events), so use {Run1,Run2,Run4}.
  const char* runs[]={"Run1_fhc","Run2_fhc","Run4_fhc"};
  const int NR=3;
  const double SIGPOT=7.655314e21;      // combined MC POT (Run1+Run2+Run4)
  const double DATAPOT=6.626e20;        // total FHC data POT, nue convention:
                                        // Run1=3.283 (std trigger) + Run2 1.268 + Run4 2.075 (excl Run5)
  const double POTSCALE=DATAPOT/SIGPOT; // scale fake data to the real data exposure
  const double mmu=0.10566, MA=1.12, MAp=1.34, MAm=0.90;
  const double bx=0.47417, by=0.06980, bz=0.87766;   // mean NuMI beam dir (from raw)

  for(int sign=0; sign<2; ++sign){
    double MAv = sign==0 ? MAp : MAm;
    const char* tag = sign==0 ? "plus" : "minus";
    TString outname = Form("%sxsec-ana-fakedata_MaRES_%s_run15.root",PROC,tag);
    printf("=== building %s (MA=%.2f) ===\n", outname.Data(), MAv);

    // ---- fast-clone all-but-tuned_cv_weight from the 5 inputs ----
    TChain cin("stv_tree");
    for(int r=0;r<NR;r++) cin.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,runs[r]));
    cin.SetBranchStatus("*",1);
    cin.SetBranchStatus("tuned_cv_weight",0);
    TFile* out=new TFile(outname,"recreate");
    out->SetCompressionLevel(1);
    TTree* ot=cin.CloneTree(-1,"fast");
    printf("  fast-cloned %lld entries\n", ot->GetEntries());

    // ---- add reweighted tuned_cv_weight ----
    TChain cw("stv_tree");
    for(int r=0;r<NR;r++) cw.Add(Form("%sxsec-ana-%s_new_numi_flux_fhc_pandora_ntuple.root",PROC,runs[r]));
    cw.SetBranchStatus("*",0);
    for(auto b:{"mc_pdg","mc_E","mc_px","mc_py","mc_pz","mc_ccnc","mc_interaction","mc_nu_energy","is_mc","tuned_cv_weight"})
      cw.SetBranchStatus(b,1);
    std::vector<int>*pdg=0; std::vector<float>*E=0,*px=0,*py=0,*pz=0;
    int ccnc,inter; bool ismc; float nuE,otcv_in;
    cw.SetBranchAddress("mc_pdg",&pdg);cw.SetBranchAddress("mc_E",&E);
    cw.SetBranchAddress("mc_px",&px);cw.SetBranchAddress("mc_py",&py);cw.SetBranchAddress("mc_pz",&pz);
    cw.SetBranchAddress("mc_ccnc",&ccnc);cw.SetBranchAddress("mc_interaction",&inter);
    cw.SetBranchAddress("mc_nu_energy",&nuE);cw.SetBranchAddress("is_mc",&ismc);
    cw.SetBranchAddress("tuned_cv_weight",&otcv_in);

    float otcv; TBranch* nb=ot->Branch("tuned_cv_weight",&otcv,"tuned_cv_weight/F");
    Long64_t Nt=ot->GetEntries();
    long nres=0; double sumcv=0,sumw=0;
    for(Long64_t i=0;i<Nt;i++){
      cw.GetEntry(i);
      double w=1.0;
      if(ismc && ccnc==0 && inter==1 && pdg){
        int mi=-1; for(size_t j=0;j<pdg->size();j++) if(abs((*pdg)[j])==13){mi=j;break;}
        if(mi>=0){
          double Emu=(*E)[mi], pl=(*px)[mi]*bx+(*py)[mi]*by+(*pz)[mi]*bz;
          double Q2=2.0*nuE*(Emu-pl)-mmu*mmu;
          if(Q2>0){ double b0=1.0+Q2/(MA*MA); w=pow(b0/(1.0+Q2/(MAv*MAv)),4); }
          nres++;
        }
      }
      otcv = otcv_in * (float)w * (float)POTSCALE;
      if(ismc && std::isfinite(otcv_in)){ sumcv+=otcv_in; sumw+=otcv_in*w; }  // guarded, POT-scale-free ratio
      nb->Fill();
    }
    printf("  reweighted RES-CC=%ld  M_A^RES rate change=%+.2f%%  (POTscale=%.5f)\n",nres,100*(sumw/sumcv-1),POTSCALE);
    ot->Write("",TObject::kOverwrite);
    TParameter<float> sp("summed_pot",(float)DATAPOT); sp.Write("summed_pot",TObject::kOverwrite);
    out->Close();
    printf("  wrote %s\n", outname.Data());
  }
  printf("DONE build_fakedata_mares\n");
}
