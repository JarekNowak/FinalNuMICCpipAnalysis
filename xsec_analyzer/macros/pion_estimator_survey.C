// pion_estimator_survey.C -- range vs calorimetry, split by golden status.
// Calorimetry saturates too (1.136 -> 0.383) and is 3-5x noisier. Golden pions read
// 0.873 -> 0.327, and above 0.8 GeV/c are indistinguishable from non-golden.
// Calorimetric energy as a pion momentum estimator, against range and MCS.
// For a STOPPING pion the deposited charge measures the kinetic energy directly; the
// question is whether it degrades more gracefully than range when the pion interacts.
#include <vector>
void calo_pion() {   // now split by golden/non-golden
  const double mmu=0.10566, mpi=0.13957;
  const char* files[3]={
    "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"};
  const int NB=6; double edge[NB+1]={0.175,0.250,0.320,0.420,0.550,0.800,1.500};
  std::vector<double> n(NB,0), sr(NB,0),sr2(NB,0), sc(NB,0),sc2(NB,0), nc(NB,0), gold(NB,0);
  std::vector<double> gr(NB,0),gr2(NB,0),ur(NB,0),ur2(NB,0),nu(NB,0);
  long bad=0, tot=0;
  for(int fi=0;fi<3;++fi){
    TFile* f=TFile::Open(files[fi]); if(!f||f->IsZombie())continue;
    TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    t->SetBranchStatus("*",0);
    for(auto b:{"backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
                "backtracked_purity","mc_pdg","mc_px","mc_py","mc_pz","mc_end_p",
                "trk_score_v","trk_len_v","trk_range_muon_mom_v","trk_calo_energy_y_v"})
      t->SetBranchStatus(b,1);
    std::vector<int>*bpdg=0,*mpdg=0;
    std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0,*mep=0,*ts=0,*tl=0,*rm=0,*ce=0;
    t->SetBranchAddress("backtracked_pdg",&bpdg);t->SetBranchAddress("backtracked_px",&bpx);
    t->SetBranchAddress("backtracked_py",&bpy);  t->SetBranchAddress("backtracked_pz",&bpz);
    t->SetBranchAddress("backtracked_purity",&bpur);
    t->SetBranchAddress("mc_pdg",&mpdg);t->SetBranchAddress("mc_px",&mpx);
    t->SetBranchAddress("mc_py",&mpy);  t->SetBranchAddress("mc_pz",&mpz);
    t->SetBranchAddress("mc_end_p",&mep);
    t->SetBranchAddress("trk_score_v",&ts);t->SetBranchAddress("trk_len_v",&tl);
    t->SetBranchAddress("trk_range_muon_mom_v",&rm);
    t->SetBranchAddress("trk_calo_energy_y_v",&ce);
    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
      if(!bpdg||!mpdg)continue;
      for(size_t j=0;j<bpdg->size();++j){
        if(abs(bpdg->at(j))!=211)continue;
        if(bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.)continue;
        double pr=rm->at(j); if(pr<=0.)continue;
        int best=-1;
        for(size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(j))continue;
          if(fabs(mpx->at(k)-bpx->at(j))+fabs(mpy->at(k)-bpy->at(j))+fabs(mpz->at(k)-bpz->at(j))<1e-4){best=k;break;} }
        if(best<0)continue;
        double pt=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
        if(pt<edge[0]||pt>=edge[NB])continue;
        int b=-1; for(int k=0;k<NB;++k) if(pt>=edge[k]&&pt<edge[k+1]) b=k;
        if(b<0)continue;
        ++tot;
        double E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
        double p_range=sqrt(std::max(0.,E*E-mpi*mpi));
        // the branch carries occasional garbage (full-range RMS 3e4 MeV); keep physical values
        double raw = ce? ce->at(j) : -1.;
        double KE = ( raw > 0. && raw < 2000. ) ? raw/1000. : -1.;   // MeV -> GeV
        double p_calo = KE>0 ? sqrt(std::max(0.,pow(KE+mpi,2)-mpi*mpi)) : -1;
        if(p_calo<=0) ++bad;
        n[b]++; sr[b]+=p_range/pt; sr2[b]+=pow(p_range/pt,2);
        if(p_calo>0){ sc[b]+=p_calo/pt; sc2[b]+=pow(p_calo/pt,2); nc[b]++; }
        bool isg = mep->at(best)<0.01;
        if(isg){ gold[b]++; gr[b]+=p_range/pt; gr2[b]+=pow(p_range/pt,2); }
        else   { nu[b]++;  ur[b]+=p_range/pt; ur2[b]+=pow(p_range/pt,2); }
      } }
    f->Close();
  }
  printf("\n  %ld pion tracks (calorimetric energy unfilled for %.1f%%)\n",tot,100.*bad/tot);
  printf("\n  true p_pi        N   golden%%   GOLDEN <range/true> RMS   NON-GOLDEN <range/true> RMS\n");
  for(int b=0;b<NB;++b){ if(!n[b])continue;
    double mr=sr[b]/n[b], rr=sqrt(std::max(0.,sr2[b]/n[b]-mr*mr));
    double mc=nc[b]?sc[b]/nc[b]:0, rc=nc[b]?sqrt(std::max(0.,sc2[b]/nc[b]-mc*mc)):0;
    double mg=gold[b]?gr[b]/gold[b]:0, rg=gold[b]?sqrt(std::max(0.,gr2[b]/gold[b]-mg*mg)):0;
    double mu=nu[b]?ur[b]/nu[b]:0,   ru=nu[b]?sqrt(std::max(0.,ur2[b]/nu[b]-mu*mu)):0;
    printf("  %5.3f-%5.3f %6.0f   %5.1f       %6.3f %6.3f          %6.3f %6.3f\n",
      edge[b],edge[b+1],n[b],100*gold[b]/n[b],mg,rg,mu,ru); }
}
