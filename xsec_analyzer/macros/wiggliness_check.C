// wiggliness_check.C — does trk_avg_deflection_stdev_v separate stopping ("golden")
// pions from ones that interacted hadronically? Truth proxy for golden: the mass-corrected
// range momentum recovers the true pion momentum to within 20%. Runs on the RAW PeLEE
// ntuple (needs backtracked_* and the deflection branches).
//   usage: root -l -b -q macros/wiggliness_check.C
void wig3(){
  TFile f("/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t=(TTree*)f.Get("nuselection/NeutrinoSelectionFilter");
  t->SetBranchStatus("*",0);
  for(auto b:{"backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
              "trk_range_muon_mom_v","trk_avg_deflection_stdev_v","trk_score_v",
              "trk_len_v","backtracked_purity"}) t->SetBranchStatus(b,1);
  std::vector<int>*pdg=0; std::vector<float>*px=0,*py=0,*pz=0,*rm=0,*ws=0,*ts=0,*len=0,*pur=0;
  t->SetBranchAddress("backtracked_pdg",&pdg); t->SetBranchAddress("backtracked_px",&px);
  t->SetBranchAddress("backtracked_py",&py);   t->SetBranchAddress("backtracked_pz",&pz);
  t->SetBranchAddress("trk_range_muon_mom_v",&rm);
  t->SetBranchAddress("trk_avg_deflection_stdev_v",&ws);
  t->SetBranchAddress("trk_score_v",&ts); t->SetBranchAddress("trk_len_v",&len);
  t->SetBranchAddress("backtracked_purity",&pur);
  const double mmu=0.10566,mpi=0.13957;
  TH1D hg("hg","",2000,0,0.1), hn("hn","",2000,0,0.1);
  for(Long64_t i=0,N=t->GetEntries();i<N;i++){
    t->GetEntry(i); if(!pdg) continue;
    for(size_t j=0;j<pdg->size();j++){
      if(abs(pdg->at(j))!=211) continue;
      if(pur->at(j)<0.5||ts->at(j)<0.5||len->at(j)<5) continue;
      double pr=rm->at(j), w=ws->at(j); if(pr<=0||w<0) continue;
      double pt=sqrt(px->at(j)*px->at(j)+py->at(j)*py->at(j)+pz->at(j)*pz->at(j));
      if(pt<0.1) continue;
      double E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
      double pc=sqrt(std::max(0.,E*E-mpi*mpi)), r=pc/pt;
      if(r>0.8&&r<1.2) hg.Fill(w); else if(r<0.6) hn.Fill(w);
    }
  }
  double tg=hg.Integral(), tn=hn.Integral();
  printf("  golden-like %.0f tracks, interacting %.0f tracks\n",tg,tn);
  printf("  starting golden fraction: %.1f%%\n\n",100*tg/(tg+tn));
  printf("  cut w > X      gold-eff   inter-rej   golden fraction\n");
  for(double c : {0.0005,0.001,0.002,0.003,0.004,0.005,0.0075}){
    int b=hg.FindBin(c);
    double eg=hg.Integral(b,2001)/tg, en=hn.Integral(b,2001)/tn;
    printf("   %.4f       %6.1f%%     %6.1f%%      %.1f%% -> %.1f%%\n",
      c,100*eg,100*(1-en),100*tg/(tg+tn),100*(eg*tg)/(eg*tg+en*tn));
  }
}
