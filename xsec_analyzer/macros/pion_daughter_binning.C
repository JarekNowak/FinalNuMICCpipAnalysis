// pion_daughter_binning.C -- does adding the interaction-product energy change the p_pi
// binning verdict? No. Upper diagonals triple (5.5% -> 15.7%) but need 68%; the two-bin
// scheme trades a degraded first bin (97.7 -> 82.4) for a better second (61.4 -> 71.0).
// Does adding the interaction-product energy change the p_pi BINNING assessment?
// Estimator: KE from the mass-corrected range momentum, plus the calorimetric energy of the
// pion's hierarchy daughters within R. Reports the migration diagonals for the original
// five-bin scheme and the adopted two-bin one, against plain range.
#include <vector>
void daughter_bin(double R=20.) {
  const double mpi=0.13957, mmu=0.10566;
  const char* files[1]={
    "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root"};


  const int N5=5; double e5[N5+1]={0.175,0.250,0.320,0.420,0.550,1e9};
  long n5[N5]={0}, dR5[N5]={0}, dD5[N5]={0};
  long n2[2]={0}, dR2[2]={0}, dD2[2]={0};
  double sR[N5]={0}, sD[N5]={0};
  auto bin=[&](double v,double* e,int n){ for(int k=0;k<n;++k) if(v>=e[k]&&v<e[k+1]) return k; return v<e[0]?0:-1; };

  for(int fi=0;fi<1;++fi){
    TFile* f=TFile::Open(files[fi]); if(!f||f->IsZombie())continue;
    TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    std::vector<int>*bpdg=0,*mpdg=0; std::vector<unsigned int>*gen=0;
    std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0;
    std::vector<float> *ts=0,*tl=0,*rm=0,*ce=0,*se=0,*sx=0,*sy=0,*sz=0,*ex=0,*ey=0,*ez=0;
    t->SetBranchAddress("backtracked_pdg",&bpdg);t->SetBranchAddress("backtracked_px",&bpx);
    t->SetBranchAddress("backtracked_py",&bpy);  t->SetBranchAddress("backtracked_pz",&bpz);
    t->SetBranchAddress("backtracked_purity",&bpur);
    t->SetBranchAddress("mc_pdg",&mpdg);t->SetBranchAddress("mc_px",&mpx);
    t->SetBranchAddress("mc_py",&mpy);  t->SetBranchAddress("mc_pz",&mpz);
    t->SetBranchAddress("trk_score_v",&ts);t->SetBranchAddress("trk_len_v",&tl);
    t->SetBranchAddress("trk_range_muon_mom_v",&rm);
    t->SetBranchAddress("trk_calo_energy_y_v",&ce);t->SetBranchAddress("shr_energy_y_v",&se);
    t->SetBranchAddress("trk_sce_start_x_v",&sx);t->SetBranchAddress("trk_sce_start_y_v",&sy);
    t->SetBranchAddress("trk_sce_start_z_v",&sz);
    t->SetBranchAddress("trk_sce_end_x_v",&ex);t->SetBranchAddress("trk_sce_end_y_v",&ey);
    t->SetBranchAddress("trk_sce_end_z_v",&ez);
    t->SetBranchAddress("pfp_generation_v",&gen);
    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
      if(!bpdg||!mpdg)continue;
      for(size_t j=0;j<bpdg->size();++j){
        if(abs(bpdg->at(j))!=211)continue;
        if(bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.)continue;
        if(rm->at(j)<=0.)continue;
        int best=-1;
        for(size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(j))continue;
          if(fabs(mpx->at(k)-bpx->at(j))+fabs(mpy->at(k)-bpy->at(j))+fabs(mpz->at(k)-bpz->at(j))<1e-4){best=k;break;} }
        if(best<0)continue;
        double pt=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
        if(pt<e5[0])continue;
        double pr=rm->at(j), E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
        double p_range=sqrt(std::max(0.,E*E-mpi*mpi));
        // add the hierarchy daughters' visible energy to the kinetic energy
        double KE=sqrt(p_range*p_range+mpi*mpi)-mpi, add=0.;
        int gpi=(gen&&j<gen->size())?gen->at(j):-1;
        TVector3 end(ex->at(j),ey->at(j),ez->at(j));
        for(size_t k=0;k<bpdg->size();++k){
          if(k==j) continue;
          if(gen&&k<gen->size()&&gpi>0&&(int)gen->at(k)!=gpi+1) continue;
          TVector3 st(sx->at(k),sy->at(k),sz->at(k));
          if((st-end).Mag()>R) continue;
          if(ce&&ce->at(k)>0.&&ce->at(k)<2000.) add+=ce->at(k)/1000.;
          else if(se&&se->at(k)>0.&&se->at(k)<5.) add+=se->at(k);
        }
        double p_dau=sqrt(std::max(0.,pow(KE+add+mpi,2)-mpi*mpi));
        int tb=bin(pt,e5,N5); if(tb<0)continue;
        int rb=bin(p_range,e5,N5), db=bin(p_dau,e5,N5);
        n5[tb]++; if(rb==tb)dR5[tb]++; if(db==tb)dD5[tb]++;
        sR[tb]+=p_range/pt; sD[tb]+=p_dau/pt;
        int t2=pt<0.205?0:1, r2=p_range<0.205?0:1, d2=p_dau<0.205?0:1;
        n2[t2]++; if(r2==t2)dR2[t2]++; if(d2==t2)dD2[t2]++;
      } }
    f->Close();
  }
  printf("\n  daughters within %.0f cm added to the kinetic energy\n", R);
  printf("\n  true p_pi        N    <range/true> diag    <range+dau/true> diag\n");
  for(int b=0;b<N5;++b){ if(!n5[b])continue;
    printf("  %5.3f-%-7.3f %6ld   %6.3f  %5.1f%%      %6.3f  %5.1f%%\n",
      e5[b], e5[b+1]>1e8?1.5:e5[b+1], n5[b], sR[b]/n5[b],100.*dR5[b]/n5[b],
      sD[b]/n5[b],100.*dD5[b]/n5[b]); }
  printf("\n  two-bin scheme (0.205):  range %5.1f%% / %5.1f%%   range+daughters %5.1f%% / %5.1f%%\n",
    100.*dR2[0]/n2[0],100.*dR2[1]/n2[1],100.*dD2[0]/n2[0],100.*dD2[1]/n2[1]);
}
