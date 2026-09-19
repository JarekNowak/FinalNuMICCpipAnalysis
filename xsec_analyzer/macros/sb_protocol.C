// sb_protocol.C -- the statistics of the frozen control-region procedure (analysis note,
// "Frozen procedure for opening the control regions"), exercised on the full-stack pseudo-data:
//  (a) per-region normalisation pulls, against the full pre-fit covariance and against the
//      statistical + detector terms alone;
//  (b) normalisation-marginalised shape chi2 of every frozen distribution
//      (configs/sideband_plots.txt), with the covariance projected onto shape;
//  (c) the GLOBAL normalisation chi2 over the eight (region, mode) totals with the full 8x8
//      pre-fit covariance -- the single pass/fail statistic -- with its power: the per-region
//      shift and the common normalisation shift that fail it at 95%.
// Prediction = neutrino MC (per-run POT scaled, CV weight, safe-weight rule) + per-run beam-off
// (gate-scaled) + dirt. Covariance: flux / cross-section / re-interaction multisims (index-matched
// across FHC and RHC, so cross-mode correlations are carried), the eight detector-variation knobs
// as fully correlated unisims from the Run-4 detVar samples, MC, beam-off and dirt statistics, and
// the data Poisson term. Pseudo-data: the released per-run fake-data throws of the neutrino MC in
// each region plus Poisson throws of the beam-off and dirt expectations (seed 20260903, as
// sideband_compare.C mode fakestack).
//   root -l -b -q macros/sb_protocol.C
#include "sb_guard.h"   // blinding guard on every beam-on file (2026-09-06)
#include <vector>
#include <string>
#include <cmath>
#include "TMatrixD.h"
#include "TDecompSVD.h"
struct Src { std::string file; double scale; };
static void add_ext(std::vector<Src>& v, const char* m){
  const char* E="/data/uboone/processed/ext_perrun/xsec-ana-";
  const char* R1="neutrinoselection_filt_run1_beamoff.root", *R3B="neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root","numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5="numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCCX=0.98;
  if(std::string(m)=="fhc"){ v.push_back({std::string(E)+R1,OCCX*5748692./G1}); v.push_back({std::string(E)+R1,OCCX*3535129./G1}); for(auto f:R4) v.push_back({std::string(E)+f,OCCX*4131149./G4}); v.push_back({std::string(E)+R5,OCCX*5154196./G5}); }
  else { v.push_back({std::string(E)+R1,OCCX*1458253./G1}); v.push_back({std::string(E)+R1,OCCX*5422907./G1}); v.push_back({std::string(E)+R3B,OCCX*10349610./G3}); for(auto f:R4) v.push_back({std::string(E)+f,OCCX*6304167./G4}); }
}
struct Var { int reg; std::string name, br; std::vector<double> edges; bool overflowTop; };
void sb_protocol(){
  const char* P="/data/uboone/processed/sb/"; const char* DV="/data/uboone/processed/";
  const int NM=2, NR=4; const char* modes[NM]={"fhc","rhc"}; const char* dvtag[NM]={"run4fhc","run4rhc"};
  const char* rb[NR]={"CC1mu1piXp_sb_cc0pi","CC1mu1piXp_sb_multipi","CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_cosmic"}; const char* rn[NR]={"CC0pi","multipi","pi0","cosmic"};
  const char* knobs[8]={"LYdown","LYrayl","Recomb2","SCE","WMAngleXZ","WMAngleYZ","WMX","WMYZ"};
  const char* SY[3]={"weight_ppfx_all","weight_All_UBGenie","weight_reint_all"};
  std::vector<double> PMU={0.15,0.35,0.55,0.75,0.95,1.25,1.75,3.0}, CTH={-1,0.45,0.65,0.8,0.9,1};
  auto uni=[](double lo,double hi,int n){ std::vector<double> e; for(int i=0;i<=n;i++) e.push_back(lo+(hi-lo)*i/n); return e; };
  std::vector<Var> vars={{0,"cc0pi_pmu","CC1mu1piXp_candidate_muon_mom_reco",PMU,true},{0,"cc0pi_costh","CC1mu1piXp_candidate_muon_costh_reco",CTH,false},
    {1,"multipi_pmu","CC1mu1piXp_candidate_muon_mom_reco",PMU,true},{1,"multipi_costh","CC1mu1piXp_candidate_muon_costh_reco",CTH,false},{1,"multipi_oa","CC1mu1piXp_mu_pi_opening_angle",uni(0,3.15,6),false},
    {2,"pi0_pmu","CC1mu1piXp_candidate_muon_mom_reco",PMU,true},{2,"pi0_costh","CC1mu1piXp_candidate_muon_costh_reco",CTH,false},{2,"pi0_shr","shr_energy_cali",uni(0,1.0,10),true},
    {3,"cosmic_pmu","CC1mu1piXp_candidate_muon_mom_reco",PMU,true},{3,"cosmic_costh","CC1mu1piXp_candidate_muon_costh_reco",CTH,false},{3,"cosmic_oa","CC1mu1piXp_mu_pi_opening_angle",uni(2.6,3.15,4),false}};
  const int NV=vars.size();
  auto binof=[&](const Var& v,double x)->int{ int nb=v.edges.size()-1; if(x<v.edges[0]) return -1; for(int b=0;b<nb;b++) if(x<v.edges[b+1]) return b; return v.overflowTop? nb-1 : -1; };
  // accumulators: index k=m*NR+r for totals; per var: [m][v][bin]
  const int NK=NM*NR; double cvT[NK]={0}, w2T[NK]={0}, extT[NK]={0}, ext2T[NK]={0}, dirtT[NK]={0}, dirt2T[NK]={0}, dataT[NK]={0};
  std::vector<std::vector<double>> UT[3]; int NU[3]={0,0,0};            // UT[s][u][k]
  std::vector<std::vector<double>> cvH(NM), w2H(NM), extH(NM), ext2H(NM), dirtH(NM), dirt2H(NM), dataH(NM); // [m][v-bin flat]
  std::vector<int> voff(NV+1,0); for(int v=0;v<NV;v++) voff[v+1]=voff[v]+vars[v].edges.size()-1; const int NB=voff[NV];
  for(int m=0;m<NM;m++){ cvH[m].assign(NB,0); w2H[m].assign(NB,0); extH[m].assign(NB,0); ext2H[m].assign(NB,0); dirtH[m].assign(NB,0); dirt2H[m].assign(NB,0); dataH[m].assign(NB,0); }
  std::vector<std::vector<double>> UH[3];                               // UH[s][u][m*NB+i]
  double knT[8][NK]={{0}}, knCVT[NK]={0}; std::vector<std::vector<double>> knH[8]; std::vector<std::vector<double>> knCVH(NM);
  for(int kk=0;kk<8;kk++){ knH[kk].assign(NM,std::vector<double>(NB,0)); } for(int m=0;m<NM;m++) knCVH[m].assign(NB,0);
  auto safe=[](double w){ return (std::isfinite(w)&&w>=0&&w<=30)?w:1.0; };
  // generic loop over a file: fills totals and hists for the regions, weighted by w (CV or unit) and scale
  auto loop=[&](const char* file, double scale, bool weighted, bool universes, int m,
                double* T, double* T2, std::vector<double>& H, std::vector<double>& H2){
    TFile f(file); TTree* t=(TTree*)f.Get("stv_tree"); if(!t){ printf("missing %s\n",file); return; }
    t->Draw(">>el","CC1mu1piXp_sb_cc0pi||CC1mu1piXp_sb_multipi||CC1mu1piXp_sb_pi0||CC1mu1piXp_sb_cosmic","entrylist"); TEntryList* el=(TEntryList*)gDirectory->Get("el");
    // read only what is needed: the trees carry hundreds of branches including the universe vectors
    t->SetBranchStatus("*",0); for(int r=0;r<NR;r++) t->SetBranchStatus(rb[r],1);
    if(weighted) for(auto b:{"tuned_cv_weight","ppfx_cv_weight","normalisation_weight"}) t->SetBranchStatus(b,1);
    if(universes) for(int s=0;s<3;s++) t->SetBranchStatus(SY[s],1);
    for(int v=0;v<NV;v++) t->SetBranchStatus(vars[v].br.c_str(),1);
    bool fl[NR]; for(int r=0;r<NR;r++) t->SetBranchAddress(rb[r],&fl[r]);
    float tune=1,ppfx=1,norm=1; if(weighted){ t->SetBranchAddress("tuned_cv_weight",&tune); t->SetBranchAddress("ppfx_cv_weight",&ppfx); t->SetBranchAddress("normalisation_weight",&norm); }
    std::vector<double>* wv[3]={nullptr,nullptr,nullptr}; if(universes) for(int s=0;s<3;s++) t->SetBranchAddress(SY[s],&wv[s]);
    std::vector<TLeaf*> lf(NV); for(int v=0;v<NV;v++) lf[v]=t->GetLeaf(vars[v].br.c_str());
    for(Long64_t e=0;e<el->GetN();e++){ t->GetEntry(el->GetEntry(e));
      double w0=weighted? safe(tune*ppfx*norm):1.0; double w=w0*scale;
      std::vector<int> bins(NV); for(int v=0;v<NV;v++) bins[v]= fl[vars[v].reg] ? binof(vars[v],lf[v]->GetValue()) : -2;
      for(int r=0;r<NR;r++) if(fl[r]){ int k=m*NR+r; T[k]+=w; T2[k]+=w*w; }
      for(int v=0;v<NV;v++) if(bins[v]>=0){ int i=voff[v]+bins[v]; H[i]+=w; H2[i]+=w*w; }
      if(universes) for(int s=0;s<3;s++){ if(!wv[s]) continue; size_t nu=wv[s]->size(); if(NU[s]==0){ NU[s]=nu; UT[s].assign(nu,std::vector<double>(NK,0.)); UH[s].assign(nu,std::vector<double>(NM*NB,0.)); }
        for(size_t u=0;u<nu&&u<(size_t)NU[s];u++){ double wu=(*wv[s])[u]; double ww=(s==0)?wu*tune*norm:(s==1)?wu*ppfx*norm:wu*tune*ppfx*norm; ww=safe(ww)*scale;
          for(int r=0;r<NR;r++) if(fl[r]) UT[s][u][m*NR+r]+=ww; for(int v=0;v<NV;v++) if(bins[v]>=0) UH[s][u][m*NB+voff[v]+bins[v]]+=ww; } }
    }
    delete el; };
  for(int m=0;m<NM;m++){
    std::vector<Src> mc, ext; std::vector<std::string> data; double sc_dirt;
    if(m==0){ const char* rnm[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
      for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rnm[i]+".root",sc[i]}); for(auto r:{"run1","run2","run4","run5"}) data.push_back(std::string(P)+"xsec-ana-fakedata_fhc_"+r+".root"); sc_dirt=0.092402*0.65; }
    else { const char* rnm[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
      for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rnm[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]}); for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066});
      for(auto r:{"run1","run2","run3","run4"}) data.push_back(std::string(P)+"xsec-ana-fakedata_rhc_"+r+".root"); sc_dirt=0.071666*0.65; }
    add_ext(ext,modes[m]);
    for(auto& x:mc) loop(x.file.c_str(),x.scale,true,true,m,cvT,w2T,cvH[m],w2H[m]);
    for(auto& x:ext) loop(x.file.c_str(),x.scale,false,false,m,extT,ext2T,extH[m],ext2H[m]);
    loop(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P),sc_dirt,true,false,m,dirtT,dirt2T,dirtH[m],dirt2H[m]);
    double dummy2[NK]={0}; std::vector<double> dh2(NB,0); for(auto& d:data){ sb_guard_data(d); loop(d.c_str(),1.0,false,false,m,dataT,dummy2,dataH[m],dh2); }
    // detector variations: CV + 8 knobs (unscaled; only ratios are used)
    double kT[NK]={0}, k2[NK]={0}; std::vector<double> kH(NB,0), kH2(NB,0);
    loop(Form("%sxsec-ana-detvar_%s_CV.root",DV,dvtag[m]),1.0,true,false,m,knCVT,k2,knCVH[m],kH2);
    for(int kk=0;kk<8;kk++){ std::fill(kH.begin(),kH.end(),0.); std::fill(kH2.begin(),kH2.end(),0.); double t1[NK]={0},t2[NK]={0}; loop(Form("%sxsec-ana-detvar_%s_%s.root",DV,dvtag[m],knobs[kk]),1.0,true,false,m,t1,t2,kH,kH2);
      for(int r=0;r<NR;r++) knT[kk][m*NR+r]=t1[m*NR+r]; knH[kk][m]=kH; }
  }
  // ---- pseudo-data: neutrino throw (files) + Poisson throws of beam-off and dirt (seed as sideband_compare fakestack)
  TRandom3 rng(20260903); double Nd[NK], pred[NK]; std::vector<double> dH(NM*NB), pH(NM*NB);
  for(int k=0;k<NK;k++){ Nd[k]=dataT[k]+rng.Poisson(extT[k])+rng.Poisson(dirtT[k]); pred[k]=cvT[k]+extT[k]+dirtT[k]; }
  for(int m=0;m<NM;m++) for(int i=0;i<NB;i++){ dH[m*NB+i]=dataH[m][i]+rng.Poisson(extH[m][i])+rng.Poisson(dirtH[m][i]); pH[m*NB+i]=cvH[m][i]+extH[m][i]+dirtH[m][i]; }
  // ---- covariance builders
  auto cov_from=[&](int n, auto getcv, auto getU, auto getknob, auto getstat, TMatrixD& Vsys, TMatrixD& Vstatdet){
    Vsys.ResizeTo(n,n); Vsys.Zero(); Vstatdet.ResizeTo(n,n); Vstatdet.Zero();
    for(int s=0;s<3;s++){ int nu=NU[s]; for(int u=0;u<nu;u++){ std::vector<double> d(n); for(int i=0;i<n;i++) d[i]=getU(s,u,i)-getcv(i); for(int i=0;i<n;i++) for(int j=0;j<n;j++) Vsys(i,j)+=d[i]*d[j]/nu; } }
    for(int kk=0;kk<8;kk++){ std::vector<double> d(n); for(int i=0;i<n;i++) d[i]=getknob(kk,i); for(int i=0;i<n;i++) for(int j=0;j<n;j++){ Vsys(i,j)+=d[i]*d[j]; Vstatdet(i,j)+=d[i]*d[j]; } }
    for(int i=0;i<n;i++){ double st=getstat(i); Vsys(i,i)+=st; Vstatdet(i,i)+=st; } };
  // totals (8x8)
  TMatrixD Vsys, Vsd;
  cov_from(NK,[&](int i){return cvT[i];},[&](int s,int u,int i){return UT[s][u][i];},
    [&](int kk,int i){ return knCVT[i]>0? (knT[kk][i]-knCVT[i])/knCVT[i]*cvT[i] : 0.; },
    [&](int i){ return w2T[i]+ext2T[i]+dirt2T[i]+Nd[i]; }, Vsys, Vsd);
  TMatrixD V(Vsys); TMatrixD Vinv(V); Vinv.Invert();
  printf("\n==== CONTROL-REGION PROTOCOL STATISTICS on the full-stack pseudo-data (seed 20260903) ====\n");
  printf("%-5s %-8s %9s %9s %8s %8s | %10s %10s %14s | %8s\n","mode","region","N_data","N_pred","sig%","EXT%","pull(full)","pull(st+det)","sigma(full)/N","fail@95%");
  double chi2=0; for(int i=0;i<NK;i++) for(int j=0;j<NK;j++) chi2+=(Nd[i]-pred[i])*Vinv(i,j)*(Nd[j]-pred[j]);
  for(int k=0;k<NK;k++){ int m=k/NR, r=k%NR; double d=Nd[k]-pred[k]; double sfail=sqrt(15.507/(pred[k]*pred[k]*Vinv(k,k)));
    printf("%-5s %-8s %9.1f %9.1f %8s %7.1f%% | %10.2f %10.2f %13.1f%% | %7.1f%%\n",modes[m],rn[r],Nd[k],pred[k],"",100*extT[k]/pred[k],d/sqrt(V(k,k)),d/sqrt(Vsd(k,k)),100*sqrt(V(k,k))/pred[k],100*sfail); }
  printf("GLOBAL normalisation chi2 = %.2f / %d dof, p = %.3f (pass if p>0.05)\n",chi2,NK,TMath::Prob(chi2,NK));
  { double q=0; for(int i=0;i<NK;i++) for(int j=0;j<NK;j++) q+=pred[i]*Vinv(i,j)*pred[j]; printf("blind spot: a COMMON normalisation shift of all eight totals fails the global test only beyond %.1f%%\n",100*sqrt(15.507/q)); }
  // correlation matrix of the prediction (systematic part)
  printf("prediction correlation (full pre-fit covariance):\n"); for(int i=0;i<NK;i++){ printf("  %-4s%-8s",modes[i/NR],rn[i%NR]); for(int j=0;j<NK;j++) printf("%6.2f",V(i,j)/sqrt(V(i,i)*V(j,j))); printf("\n"); }
  // ---- shape tests
  printf("\n(b) normalisation-marginalised shape tests of the frozen distributions:\n%-5s %-14s %5s %9s %8s\n","mode","distribution","bins","chi2/ndf","p");
  for(int m=0;m<NM;m++) for(int v=0;v<NV;v++){ int nb=vars[v].edges.size()-1; int o=m*NB+voff[v];
    std::vector<double> d(nb),p(nb); double sd=0,sp=0; for(int b=0;b<nb;b++){ d[b]=dH[o+b]; p[b]=pH[o+b]; sd+=d[b]; sp+=p[b]; } if(sp<=0||sd<=0) continue; double sc=sd/sp;
    TMatrixD Vs,Vsd2; int mm=m; cov_from(nb,[&](int i){return cvH[mm][voff[v]+i];},[&](int s,int u,int i){return UH[s][u][o+i];},
      [&](int kk,int i){ double c=knCVH[mm][voff[v]+i]; return c>0? (knH[kk][mm][voff[v]+i]-c)/c*cvH[mm][voff[v]+i]:0.; },
      [&](int i){ return w2H[mm][voff[v]+i]+ext2H[mm][voff[v]+i]+dirt2H[mm][voff[v]+i]+dH[o+i]/(sc*sc); }, Vs, Vsd2);
    // scale prediction (and its covariance) to the data total; project out the normalisation mode
    TMatrixD Vn(nb,nb); for(int i=0;i<nb;i++) for(int j=0;j<nb;j++) Vn(i,j)=Vs(i,j)*sc*sc; std::vector<double> ps(nb); double sps=0; for(int b=0;b<nb;b++){ ps[b]=p[b]*sc; sps+=ps[b]; }
    TMatrixD Pm(nb,nb); for(int i=0;i<nb;i++) for(int j=0;j<nb;j++) Pm(i,j)=(i==j?1.:0.)-ps[i]/sps; TMatrixD PT(TMatrixD::kTransposed,Pm); TMatrixD Vp=Pm*Vn*PT;
    TDecompSVD svd(Vp); TMatrixD Vpi(nb,nb); { TMatrixD U=svd.GetU(), Vv=svd.GetV(); TVectorD S=svd.GetSig(); double smax=S(0); Vpi.Zero(); for(int a=0;a<nb;a++){ if(S(a)<1e-10*smax) continue; for(int i=0;i<nb;i++) for(int j=0;j<nb;j++) Vpi(i,j)+=Vv(i,a)*U(j,a)/S(a); } }
    double c2=0; for(int i=0;i<nb;i++) for(int j=0;j<nb;j++) c2+=(d[i]-ps[i])*Vpi(i,j)*(d[j]-ps[j]);
    printf("%-5s %-14s %5d %6.2f/%-2d %8.3f\n",modes[m],vars[v].name.c_str(),nb,c2,nb-1,TMath::Prob(c2,nb-1)); }
  printf("\n(a) pulls: full = pre-fit prediction covariance + data Poisson; st+det = MC/EXT/dirt/data statistics + detector knobs only.\n(c) fail@95%%: fractional shift of that region's total that alone brings the global chi2 to its 95%% point (15.51/8).\n");
}
