// w_with_regression.C -- rebuild W_pipr with range, the BDTG regression, and true p_pi.
// Bias -0.127 -> -0.005 -> +0.003; RMS 0.144 -> 0.119 -> 0.046.
// Rebuild W_pipr with three pion-momentum estimators and compare against true W.
//   (a) range, mass-corrected -- the analysis estimator
//   (b) the BDTG regression
//   (c) TRUE pion momentum -- the ceiling, showing what a perfect p_pi would buy
// Proton is always the reco proton (range momentum under the proton hypothesis); directions
// are always reco. So the comparison isolates the pion-momentum contribution.
#include <vector>
void w_regr() {
  const double mmu=0.10566, mpi=0.13957, mpr=0.93827;
  float p_range,p_mcs,calo,len,wstd,wmean,wsep;
  float bragg_pion,bragg_mip,bragg_mu,bragg_p,tscore,ndau;
  TMVA::Reader r("!Color:!Silent");
  r.AddVariable("p_range",&p_range); r.AddVariable("p_mcs",&p_mcs); r.AddVariable("calo",&calo);
  r.AddVariable("len",&len); r.AddVariable("wstd",&wstd); r.AddVariable("wmean",&wmean);
  r.AddVariable("wsep",&wsep); r.AddVariable("bragg_pion",&bragg_pion);
  r.AddVariable("bragg_mip",&bragg_mip); r.AddVariable("bragg_mu",&bragg_mu);
  r.AddVariable("bragg_p",&bragg_p); r.AddVariable("tscore",&tscore); r.AddVariable("ndau",&ndau);
  r.BookMVA("BDTG","pion_regr/weights/PionMom_BDTG.weights.xml");

  // held-out run only
  TFile* f=TFile::Open("/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
  std::vector<int>*bpdg=0,*mpdg=0;
  std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0;
  std::vector<float> *ts=0,*tl=0,*rm=0,*mcs=0,*ce=0,*ws=0,*wm=0,*wsp=0;
  std::vector<float> *bpi=0,*bmip=0,*bmu=0,*bp=0,*dx=0,*dy=0,*dz=0,*ep=0;
  std::vector<unsigned int> *dtr=0,*dsh=0;
  t->SetBranchAddress("backtracked_pdg",&bpdg);t->SetBranchAddress("backtracked_px",&bpx);
  t->SetBranchAddress("backtracked_py",&bpy);  t->SetBranchAddress("backtracked_pz",&bpz);
  t->SetBranchAddress("backtracked_purity",&bpur);
  t->SetBranchAddress("mc_pdg",&mpdg);t->SetBranchAddress("mc_px",&mpx);
  t->SetBranchAddress("mc_py",&mpy);  t->SetBranchAddress("mc_pz",&mpz);
  t->SetBranchAddress("trk_score_v",&ts); t->SetBranchAddress("trk_len_v",&tl);
  t->SetBranchAddress("trk_range_muon_mom_v",&rm); t->SetBranchAddress("trk_mcs_muon_mom_v",&mcs);
  t->SetBranchAddress("trk_calo_energy_y_v",&ce);
  t->SetBranchAddress("trk_avg_deflection_stdev_v",&ws);
  t->SetBranchAddress("trk_avg_deflection_mean_v",&wm);
  t->SetBranchAddress("trk_avg_deflection_separation_mean_v",&wsp);
  t->SetBranchAddress("trk_bragg_pion_v",&bpi);t->SetBranchAddress("trk_bragg_mip_v",&bmip);
  t->SetBranchAddress("trk_bragg_mu_v",&bmu);  t->SetBranchAddress("trk_bragg_p_v",&bp);
  t->SetBranchAddress("pfp_trk_daughters_v",&dtr);t->SetBranchAddress("pfp_shr_daughters_v",&dsh);
  t->SetBranchAddress("trk_dir_x_v",&dx);t->SetBranchAddress("trk_dir_y_v",&dy);
  t->SetBranchAddress("trk_dir_z_v",&dz);t->SetBranchAddress("trk_energy_proton_v",&ep);

  TH1D hA("hA","",200,-1,1), hB("hB","",200,-1,1), hC("hC","",200,-1,1);
  // migration diagonals for the current six-bin scheme and for candidate coarser ones
  double e6[7]={1.08,1.19,1.23,1.27,1.34,1.47,2.90};
  double e3[4]={1.08,1.15,1.25,2.90};
  double e2[3]={1.08,1.15,2.90};
  long n6[6]={0},dA6[6]={0},dB6[6]={0}, n3[3]={0},dA3[3]={0},dB3[3]={0}, n2[2]={0},dA2[2]={0},dB2[2]={0};
  auto bn=[](double v,double* e,int n){ for(int k=0;k<n;++k) if(v>=e[k]&&v<e[k+1]) return k; return v<e[0]?0:n-1; };
  long nev=0;
  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
    if(!bpdg||!mpdg) continue;
    int ip=-1, ipr=-1;
    for(size_t j=0;j<bpdg->size();++j){
      if(bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.) continue;
      if(abs(bpdg->at(j))==211 && ip<0) ip=j;
      if(bpdg->at(j)==2212 && ipr<0) ipr=j;
    }
    if(ip<0||ipr<0) continue;
    if(rm->at(ip)<=0.||ep->at(ipr)<=0.) continue;
    // truth match the pion
    int best=-1;
    for(size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(ip))continue;
      if(fabs(mpx->at(k)-bpx->at(ip))+fabs(mpy->at(k)-bpy->at(ip))+fabs(mpz->at(k)-bpz->at(ip))<1e-4){best=k;break;} }
    if(best<0) continue;
    double ptrue=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
    if(ptrue<0.175) continue;
    // pion estimators
    double pr=rm->at(ip);
    double E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
    p_range=sqrt(std::max(0.,E*E-mpi*mpi));
    p_mcs=(mcs&&mcs->at(ip)>0.&&mcs->at(ip)<5.)?mcs->at(ip):-1.f;
    double raw=ce?ce->at(ip):-1.; calo=(raw>0.&&raw<2000.)?raw/1000.:-1.f;
    len=tl->at(ip); wstd=ws->at(ip); wmean=wm->at(ip); wsep=wsp->at(ip);
    auto san=[](float v){return std::isfinite(v)?v:-1.f;};
    bragg_pion=san(bpi->at(ip)); bragg_mip=san(bmip->at(ip));
    bragg_mu=san(bmu->at(ip)); bragg_p=san(bp->at(ip));
    tscore=ts->at(ip); ndau=dtr->at(ip)+dsh->at(ip);
    double p_reg=r.EvaluateRegression("BDTG")[0];
    // proton: range momentum under the proton hypothesis (trk_energy_proton is KE in GeV)
    double KEp=ep->at(ipr), Ep=KEp+mpr, ppr=sqrt(std::max(0.,Ep*Ep-mpr*mpr));
    TVector3 upi(dx->at(ip),dy->at(ip),dz->at(ip)); if(upi.Mag()>0) upi=upi.Unit();
    TVector3 upr(dx->at(ipr),dy->at(ipr),dz->at(ipr)); if(upr.Mag()>0) upr=upr.Unit();
    TVector3 tpi(mpx->at(best),mpy->at(best),mpz->at(best));
    auto W=[&](double ppi, TVector3 dirpi){
      TLorentzVector a(dirpi*ppi, sqrt(ppi*ppi+mpi*mpi));
      TLorentzVector b(upr*ppr,   sqrt(ppr*ppr+mpr*mpr));
      return (a+b).M(); };
    double Wtrue = W(ptrue, tpi.Unit());
    double WA=W(p_range,upi), WB=W(p_reg,upi);
    hA.Fill(WA-Wtrue); hB.Fill(WB-Wtrue); hC.Fill(W(ptrue,upi)-Wtrue);
    int t6=bn(Wtrue,e6,6); n6[t6]++; if(bn(WA,e6,6)==t6)dA6[t6]++; if(bn(WB,e6,6)==t6)dB6[t6]++;
    int t3=bn(Wtrue,e3,3); n3[t3]++; if(bn(WA,e3,3)==t3)dA3[t3]++; if(bn(WB,e3,3)==t3)dB3[t3]++;
    int t2=bn(Wtrue,e2,2); n2[t2]++; if(bn(WA,e2,2)==t2)dA2[t2]++; if(bn(WB,e2,2)==t2)dB2[t2]++;
    ++nev;
  }
  printf("\n  W_pipr rebuilt on %ld pion+proton events (held-out Run 2)\n", nev);
  printf("    pion estimator      W bias    W RMS\n");
  printf("    range (analysis)   %+7.3f  %7.3f\n", hA.GetMean(), hA.GetRMS());
  printf("    BDTG regression    %+7.3f  %7.3f\n", hB.GetMean(), hB.GetRMS());
  printf("    TRUE p_pi (ceiling)%+7.3f  %7.3f\n", hC.GetMean(), hC.GetRMS());
  printf("\n  migration diagonals, range -> regression\n");
  printf("    six-bin  (current):"); for(int b=0;b<6;++b) printf(" %4.0f%%",n6[b]?100.*dA6[b]/n6[b]:0); printf("\n");
  printf("                       "); for(int b=0;b<6;++b) printf(" %4.0f%%",n6[b]?100.*dB6[b]/n6[b]:0); printf("\n");
  printf("    three-bin (1.15,1.25):"); for(int b=0;b<3;++b) printf(" %4.0f%%",n3[b]?100.*dA3[b]/n3[b]:0); printf("\n");
  printf("                         "); for(int b=0;b<3;++b) printf(" %4.0f%%",n3[b]?100.*dB3[b]/n3[b]:0); printf("\n");
  printf("    two-bin  (1.15):    "); for(int b=0;b<2;++b) printf(" %4.0f%%",n2[b]?100.*dA2[b]/n2[b]:0);
  printf("   N %ld/%ld\n", n2[0],n2[1]);
  printf("                        "); for(int b=0;b<2;++b) printf(" %4.0f%%",n2[b]?100.*dB2[b]/n2[b]:0); printf("\n");
}
