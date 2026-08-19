// pion_daughter_recovery.C -- can the pion interaction products be summed to recover its
// momentum? No. With hierarchy association ~10% of the missing energy is recovered at high
// momentum, ~25% mid-range: most of it leaves as neutrons or sub-threshold nucleons.
// Step 1: is there recoverable energy near the pion endpoint at all?
// PeLEE has no parent index, so association is geometric: PFPs whose start lies within R of
// the pion track end. First just characterise what is there, before trying to use it.
#include <vector>
void daughter1() {
  const double mpi=0.13957, mmu=0.10566;
  TFile* f=TFile::Open("/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
  std::vector<int>*bpdg=0,*mpdg=0; std::vector<unsigned int>*gen=0;
  std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0,*mep=0;
  std::vector<float> *ts=0,*tl=0,*rm=0,*ce=0,*se=0;
  std::vector<float> *sx=0,*sy=0,*sz=0,*ex=0,*ey=0,*ez=0;
  t->SetBranchAddress("backtracked_pdg",&bpdg);t->SetBranchAddress("backtracked_px",&bpx);
  t->SetBranchAddress("backtracked_py",&bpy);  t->SetBranchAddress("backtracked_pz",&bpz);
  t->SetBranchAddress("backtracked_purity",&bpur);
  t->SetBranchAddress("mc_pdg",&mpdg);t->SetBranchAddress("mc_px",&mpx);
  t->SetBranchAddress("mc_py",&mpy);  t->SetBranchAddress("mc_pz",&mpz);
  t->SetBranchAddress("mc_end_p",&mep);
  t->SetBranchAddress("trk_score_v",&ts);t->SetBranchAddress("trk_len_v",&tl);
  t->SetBranchAddress("trk_range_muon_mom_v",&rm);
  t->SetBranchAddress("trk_calo_energy_y_v",&ce);
  t->SetBranchAddress("shr_energy_y_v",&se);
  t->SetBranchAddress("trk_sce_start_x_v",&sx);t->SetBranchAddress("trk_sce_start_y_v",&sy);
  t->SetBranchAddress("trk_sce_start_z_v",&sz);
  t->SetBranchAddress("trk_sce_end_x_v",&ex);  t->SetBranchAddress("trk_sce_end_y_v",&ey);
  t->SetBranchAddress("trk_sce_end_z_v",&ez);
  t->SetBranchAddress("pfp_generation_v",&gen);

  const int NR=4; double R[NR]={5.,10.,20.,40.};
  const int NB=3; double edge[NB+1]={0.175,0.420,0.800,1.500};
  double n[NB]={0}, ndau[NB][NR]={{0}}, edau[NB][NR]={{0}}, miss[NB]={0}, gold[NB]={0};
  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
    if(!bpdg||!mpdg) continue;
    for(size_t j=0;j<bpdg->size();++j){
      if(abs(bpdg->at(j))!=211) continue;
      if(bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.) continue;
      if(rm->at(j)<=0.) continue;
      int best=-1;
      for(size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(j))continue;
        if(fabs(mpx->at(k)-bpx->at(j))+fabs(mpy->at(k)-bpy->at(j))+fabs(mpz->at(k)-bpz->at(j))<1e-4){best=k;break;} }
      if(best<0) continue;
      double pt=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
      if(pt<edge[0]||pt>=edge[NB]) continue;
      int b=-1; for(int k=0;k<NB;++k) if(pt>=edge[k]&&pt<edge[k+1]) b=k;
      if(b<0) continue;
      n[b]++; if(mep->at(best)<0.01) gold[b]++;
      double KEtrue=sqrt(pt*pt+mpi*mpi)-mpi;
      double own = (ce&&ce->at(j)>0.&&ce->at(j)<2000.) ? ce->at(j)/1000. : 0.;
      miss[b] += KEtrue-own;                       // energy unaccounted for by the pion track
      TVector3 end(ex->at(j),ey->at(j),ez->at(j));
      int gpi = (gen&&j<gen->size())? gen->at(j) : -1;
      for(size_t k=0;k<bpdg->size();++k){
        if(k==j) continue;
        // hierarchy cut: only PFPs one generation below the pion count as its products
        if(gen && k<gen->size() && gpi>0 && (int)gen->at(k) != gpi+1) continue;
        TVector3 st(sx->at(k),sy->at(k),sz->at(k));
        double d=(st-end).Mag();
        double e = 0.;
        if(ce&&ce->at(k)>0.&&ce->at(k)<2000.) e=ce->at(k)/1000.;
        else if(se&&se->at(k)>0.&&se->at(k)<5.) e=se->at(k);
        for(int ri=0;ri<NR;++ri) if(d<R[ri]){ ndau[b][ri]++; edau[b][ri]+=e; }
      }
    } }
  printf("\n  Run 1: same, but restricted to PFPs one generation BELOW the pion\n");
  printf("  true p_pi      N   golden%%  <KE missing>   <n PFP within R>  <sum E within R> [GeV]\n");
  printf("                                            5cm  10cm 20cm 40cm     5cm   10cm   20cm   40cm\n");
  for(int b=0;b<NB;++b){ if(!n[b])continue;
    printf("  %5.3f-%5.3f %6.0f  %5.1f    %6.3f    ",edge[b],edge[b+1],n[b],100*gold[b]/n[b],miss[b]/n[b]);
    for(int ri=0;ri<NR;++ri) printf("%5.2f ",ndau[b][ri]/n[b]);
    printf("  ");
    for(int ri=0;ri<NR;++ri) printf("%6.3f ",edau[b][ri]/n[b]);
    printf("\n"); }
  f->Close();
}
