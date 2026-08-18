// mcs_vs_range_pion.C -- compare MCS and range momentum estimators for charged pions.
// Conclusion: MCS is unbiased where range saturates (1.011 vs 0.319 in the top bin) but
// carries a 210% RMS there and 1400% at low momentum, because our pion tracks are 25-84 cm
// and MCS needs about a metre. Unusable as a per-event estimator. See the note, Sec. ppi_bias.
//   usage: root -l -b -q "macros/mcs_vs_range_pion.C(0)"   // optional min track length in cm
// Does MCS beat range for charged pions?
// Range saturates because pions interact before stopping, truncating the track. MCS infers
// momentum from scattering along whatever track exists, so in principle it does not care how
// the track ended. The risk is resolution: MCS needs length, and our pions are short.
#include <vector>
void mcs_pion() {
  const double mmu=0.10566, mpi=0.13957;
  const char* files[3]={
    "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"};

  const int NB=6; double edge[NB+1]={0.175,0.250,0.320,0.420,0.550,0.800,1.500};
  std::vector<double> n(NB,0), sr(NB,0), sm(NB,0), sr2(NB,0), sm2(NB,0), len(NB,0);
  long ngood=0, nmcs_bad=0;

  for (int fi=0; fi<3; ++fi) {
    TFile* f=TFile::Open(files[fi]); if(!f||f->IsZombie()) continue;
    TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    t->SetBranchStatus("*",0);
    for (auto b:{"backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
                 "backtracked_purity","mc_pdg","mc_px","mc_py","mc_pz",
                 "trk_score_v","trk_len_v","trk_range_muon_mom_v","trk_mcs_muon_mom_v"})
      t->SetBranchStatus(b,1);
    std::vector<int> *bpdg=0,*mpdg=0;
    std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0,*ts=0,*tl=0,*rm=0,*mm=0;
    t->SetBranchAddress("backtracked_pdg",&bpdg); t->SetBranchAddress("backtracked_px",&bpx);
    t->SetBranchAddress("backtracked_py",&bpy);   t->SetBranchAddress("backtracked_pz",&bpz);
    t->SetBranchAddress("backtracked_purity",&bpur);
    t->SetBranchAddress("mc_pdg",&mpdg); t->SetBranchAddress("mc_px",&mpx);
    t->SetBranchAddress("mc_py",&mpy);   t->SetBranchAddress("mc_pz",&mpz);
    t->SetBranchAddress("trk_score_v",&ts); t->SetBranchAddress("trk_len_v",&tl);
    t->SetBranchAddress("trk_range_muon_mom_v",&rm);
    t->SetBranchAddress("trk_mcs_muon_mom_v",&mm);

    Long64_t N=t->GetEntries();
    for (Long64_t i=0;i<N;++i){ t->GetEntry(i);
      if(!bpdg||!mpdg) continue;
      for (size_t j=0;j<bpdg->size();++j){
        if (abs(bpdg->at(j))!=211) continue;
        if (bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.) continue;
        double pr=rm->at(j), pm=mm->at(j);
        if (pr<=0.) continue;
        int best=-1;
        for (size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(j)) continue;
          if (fabs(mpx->at(k)-bpx->at(j))+fabs(mpy->at(k)-bpy->at(j))+fabs(mpz->at(k)-bpz->at(j))<1e-4){best=k;break;} }
        if (best<0) continue;
        double pt=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
        if (pt<edge[0]||pt>=edge[NB]) continue;
        // range: kinetic-energy-preserving mu->pi swap, exactly as the bin config does it
        double E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
        double p_range=sqrt(std::max(0.,E*E-mpi*mpi));
        ++ngood; if (pm<=0.) ++nmcs_bad;
        int b=-1; for(int k=0;k<NB;++k) if(pt>=edge[k]&&pt<edge[k+1]) b=k;
        if(b<0) continue;
        n[b]++; len[b]+=tl->at(j);
        sr[b]+=p_range/pt; sr2[b]+=pow(p_range/pt,2);
        if (pm>0.){ sm[b]+=pm/pt; sm2[b]+=pow(pm/pt,2); }
      } }
    f->Close();
  }
  printf("\n  truth-matched pion tracks: %ld   (MCS unfilled for %.1f%%)\n",
         ngood, 100.*nmcs_bad/std::max(1L,ngood));
  printf("\n  true p_pi bin      N     <len>   <range/true>  RMS    <MCS/true>  RMS\n");
  for(int b=0;b<NB;++b){ if(!n[b]) continue;
    double mr=sr[b]/n[b], rr=sqrt(std::max(0.,sr2[b]/n[b]-mr*mr));
    double mm_=sm[b]/n[b], rm_=sqrt(std::max(0.,sm2[b]/n[b]-mm_*mm_));
    printf("  %5.3f-%5.3f  %7.0f  %5.0fcm   %6.3f  %6.3f   %6.3f  %6.3f\n",
           edge[b],edge[b+1],n[b],len[b]/n[b],mr,rr,mm_,rm_); }
}
