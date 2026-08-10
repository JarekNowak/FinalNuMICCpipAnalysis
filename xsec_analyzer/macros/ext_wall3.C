// ext_wall3.C — backward-projection distance to the active-volume wall for
//   (a) the muon candidate, (b) the pion candidate, and
//   (c) the COMBINED mu+pi momentum direction (p_mu*dir_mu + p_pi*dir_pi),
// projected backwards from the reconstructed vertex (muon start). Plotted+scanned
// for Signal / nu-bkg / EXT at the final selection.
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
static const double AVx0=0.0,AVx1=256.35,AVy0=-116.5,AVy1=116.5,AVz0=0.0,AVz1=1036.8;
static double back_wall(double sx,double sy,double sz,double dx,double dy,double dz){
  double vx=-dx,vy=-dy,vz=-dz,t=1e9,n=std::sqrt(dx*dx+dy*dy+dz*dz);
  if(n<1e-6) return 200; vx/=n; vy/=n; vz/=n; // ensure unit
  if(vx> 1e-9)t=std::min(t,(AVx1-sx)/vx); else if(vx<-1e-9)t=std::min(t,(AVx0-sx)/vx);
  if(vy> 1e-9)t=std::min(t,(AVy1-sy)/vy); else if(vy<-1e-9)t=std::min(t,(AVy0-sy)/vy);
  if(vz> 1e-9)t=std::min(t,(AVz1-sz)/vz); else if(vz<-1e-9)t=std::min(t,(AVz0-sz)/vz);
  return t;
}
// single pass: MC fills *s if signal else *b; EXT fills the *b slot.
// loose=true -> "early" stage: mu+pi candidate found, NO FV/containment/final cut.
static void fill(TChain&c,bool isMC,bool loose,TH1D&hms,TH1D&hmb,TH1D&hps,TH1D&hpb,TH1D&hcs,TH1D&hcb){
  bool sel=0,sig=0; int mI=0,pI=0; double mMom=0,pMom=0;
  std::vector<float>*sx=0,*sy=0,*sz=0,*dx=0,*dy=0,*dz=0;
  c.SetBranchStatus("*",0);
  for(auto b:{"CC1mu1piXp_Selected","CC1mu1piXp_CandidateMuonIndex","CC1mu1piXp_CandidatePionIndex","CC1mu1piXp_candidate_muon_mom_reco","CC1mu1piXp_candidate_pion_mom_reco","trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v","trk_dir_x_v","trk_dir_y_v","trk_dir_z_v"}) c.SetBranchStatus(b,1);
  if(isMC)c.SetBranchStatus("CC1mu1piXp_MC_Signal",1);
  c.SetBranchAddress("CC1mu1piXp_Selected",&sel);c.SetBranchAddress("CC1mu1piXp_CandidateMuonIndex",&mI);c.SetBranchAddress("CC1mu1piXp_CandidatePionIndex",&pI);
  c.SetBranchAddress("CC1mu1piXp_candidate_muon_mom_reco",&mMom);c.SetBranchAddress("CC1mu1piXp_candidate_pion_mom_reco",&pMom);
  if(isMC)c.SetBranchAddress("CC1mu1piXp_MC_Signal",&sig);
  c.SetBranchAddress("trk_sce_start_x_v",&sx);c.SetBranchAddress("trk_sce_start_y_v",&sy);c.SetBranchAddress("trk_sce_start_z_v",&sz);
  c.SetBranchAddress("trk_dir_x_v",&dx);c.SetBranchAddress("trk_dir_y_v",&dy);c.SetBranchAddress("trk_dir_z_v",&dz);
  Long64_t N=c.GetEntries();
  for(Long64_t i=0;i<N;i++){c.GetEntry(i); if(!loose && !sel)continue;
    if(!sx||!dx) continue; int ns=sx->size(); if(mI<0||mI>=ns||pI<0||pI>=ns)continue;
    bool s = isMC ? sig : false;  // EXT counted as "bkg" side
    double mx=dx->at(mI),my=dy->at(mI),mz=dz->at(mI), px=dx->at(pI),py=dy->at(pI),pz=dz->at(pI);
    double vx=sx->at(mI),vy=sy->at(mI),vz=sz->at(mI);
    double cx=mMom*mx+pMom*px, cy=mMom*my+pMom*py, cz=mMom*mz+pMom*pz;
    double dm=std::min(back_wall(vx,vy,vz,mx,my,mz),199.9);
    double dp=std::min(back_wall(sx->at(pI),sy->at(pI),sz->at(pI),px,py,pz),199.9);
    double dc=std::min(back_wall(vx,vy,vz,cx,cy,cz),199.9);
    if(s){ hms.Fill(dm); hps.Fill(dp); hcs.Fill(dc); } else { hmb.Fill(dm); hpb.Fill(dp); hcb.Fill(dc); }
  }
  c.ResetBranchAddresses();
}
static void scan(const char* tag,TH1D&hs,TH1D&hb,TH1D&he){
  printf("--- %s (backward-to-wall) ---  signal=%.0f nu-bkg=%.0f EXT=%.0f\n",tag,hs.Integral(),hb.Integral(),he.Integral());
  printf("  %-12s %8s %8s %8s\n","cut d>X","effS","effB","effE");
  for(double x:{10,20,30,50,80}){ double eS=hs.Integral(hs.FindBin(x),41)/hs.Integral(),eB=hb.Integral(hb.FindBin(x),41)/hb.Integral(),eE=he.Integral(he.FindBin(x),41)/he.Integral();
    printf("  d > %-8.0f %7.1f%% %7.1f%% %7.1f%%\n",x,100*eS,100*eB,100*eE);} printf("\n");
}
static void plot(const char* nm,const char* xt,TH1D&hs,TH1D&hb,TH1D&he){
  gStyle->SetOptStat(0);int cS=TColor::GetColor("#0072B2"),cB=TColor::GetColor("#E69F00"),cE=TColor::GetColor("#999999");
  for(TH1D*h:{&hs,&hb,&he}) if(h->Integral()>0)h->Scale(1.0/h->Integral());
  hs.SetLineColor(cS);hb.SetLineColor(cB);he.SetLineColor(cE);hs.SetLineWidth(3);hb.SetLineWidth(2);he.SetLineWidth(3);
  hs.SetFillColorAlpha(cS,0.15);he.SetFillColorAlpha(cE,0.25);
  hs.SetMaximum(std::max({hs.GetMaximum(),hb.GetMaximum(),he.GetMaximum()})*1.3);hs.SetTitle(Form(";%s;area-normalised",xt));
  TCanvas c("c","",800,600);hs.Draw("hist");hb.Draw("hist same");he.Draw("hist same");
  TLegend lg(0.55,0.72,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
  lg.AddEntry(&hs,"Signal","l");lg.AddEntry(&hb,"#nu background","l");lg.AddEntry(&he,"EXT (cosmic)","l");lg.Draw();
  c.SaveAs(Form("unfold_output/%s.pdf",nm));
}
void ext_wall3(bool loose=false){
  const char* tag = loose ? "early" : "final";
  const char* P="/data/uboone/processed/"; TChain mc("stv_tree"),ext("stv_tree");
  const char* mcf[]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc","Run1_rhc_new_numi_flux_rhc_pandora_ntuple","Run2_rhc_new_numi_flux_rhc_pandora_ntuple","Run4a_rhc_new_numi_flux_rhc_pandora_ntuple","Run4b_rhc_new_numi_flux_rhc_pandora_ntuple","Run4c_rhc_new_numi_flux_rhc_pandora_ntuple"};
  for(auto f:mcf)mc.Add(Form("%sxsec-ana-%s.root",P,f));
  for(auto s:{"aa","ab","ac","ad","ae"})mc.Add(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",P,s));
  ext.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
  TH1D ms("ms","",40,0,200),mb("mb","",40,0,200),me("me","",40,0,200);
  TH1D ps("ps","",40,0,200),pb("pb","",40,0,200),pe("pe","",40,0,200);
  TH1D cs("cs","",40,0,200),cb("cb","",40,0,200),ce("ce","",40,0,200);
  fill(mc,true,loose, ms,mb, ps,pb, cs,cb);
  fill(ext,false,loose, me,me, pe,pe, ce,ce);
  printf("STAGE = %s\n",tag);
  scan("MUON",ms,mb,me); scan("PION",ps,pb,pe); scan("COMBINED mu+pi",cs,cb,ce);
  plot(Form("ext_wall_%s_muon",tag),"muon backward dist to AV wall [cm]",ms,mb,me);
  plot(Form("ext_wall_%s_pion",tag),"pion backward dist to AV wall [cm]",ps,pb,pe);
  plot(Form("ext_wall_%s_comb",tag),"combined #mu+#pi backward dist to AV wall [cm]",cs,cb,ce);
  printf("wrote ext_wall_%s_{muon,pion,comb}.pdf\n",tag); fflush(stdout); gSystem->Exit(0);
}
