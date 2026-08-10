// ext_wallgap.C — new cosmic-discriminating variable: for the muon candidate,
// project the track BACKWARDS from its start point (opposite to its direction) and
// find the distance to where it exits the active volume (AV) wall. A cosmic muon
// enters through a wall, so its backward projection exits almost immediately (small
// distance); a neutrino muon starts at an interior vertex (larger distance).
// Computed per event (ray-box intersection); plotted+scanned for Signal/nu-bkg/EXT.
#include <vector>
#include <string>
#include <algorithm>
// Active volume (TPC) bounds, cm
static const double AVx0=0.0, AVx1=256.35, AVy0=-116.5, AVy1=116.5, AVz0=0.0, AVz1=1036.8;
static double back_wall_dist(double sx,double sy,double sz,double dx,double dy,double dz){
  // backward travel velocity = -dir (dir is unit). exit distance from S inside the box.
  double vx=-dx, vy=-dy, vz=-dz, t=1e9;
  if(vx> 1e-9) t=std::min(t,(AVx1-sx)/vx); else if(vx<-1e-9) t=std::min(t,(AVx0-sx)/vx);
  if(vy> 1e-9) t=std::min(t,(AVy1-sy)/vy); else if(vy<-1e-9) t=std::min(t,(AVy0-sy)/vy);
  if(vz> 1e-9) t=std::min(t,(AVz1-sz)/vz); else if(vz<-1e-9) t=std::min(t,(AVz0-sz)/vz);
  return t; // cm (dir unit)
}
static void fill(TChain& c, bool isMC, int want_sig, TH1D& h){
  bool sel=0, sig=0; int mIdx=0; std::vector<float> *sx=0,*sy=0,*sz=0,*dx=0,*dy=0,*dz=0;
  c.SetBranchStatus("*",0);
  for(auto b:{"CC1mu1piXp_Selected","CC1mu1piXp_CandidateMuonIndex","trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v","trk_dir_x_v","trk_dir_y_v","trk_dir_z_v"}) c.SetBranchStatus(b,1);
  if(isMC) c.SetBranchStatus("CC1mu1piXp_MC_Signal",1);
  c.SetBranchAddress("CC1mu1piXp_Selected",&sel); c.SetBranchAddress("CC1mu1piXp_CandidateMuonIndex",&mIdx);
  if(isMC) c.SetBranchAddress("CC1mu1piXp_MC_Signal",&sig);
  c.SetBranchAddress("trk_sce_start_x_v",&sx); c.SetBranchAddress("trk_sce_start_y_v",&sy); c.SetBranchAddress("trk_sce_start_z_v",&sz);
  c.SetBranchAddress("trk_dir_x_v",&dx); c.SetBranchAddress("trk_dir_y_v",&dy); c.SetBranchAddress("trk_dir_z_v",&dz);
  Long64_t N=c.GetEntries();
  for(Long64_t i=0;i<N;i++){ c.GetEntry(i);
    if(!sel) continue; if(isMC && (sig!=want_sig)) continue;
    if(mIdx<0 || mIdx>=(int)sx->size()) continue;
    double d=back_wall_dist(sx->at(mIdx),sy->at(mIdx),sz->at(mIdx),dx->at(mIdx),dy->at(mIdx),dz->at(mIdx));
    h.Fill(std::min(d,199.9));
  }
}
void ext_wallgap(){
  const char* P="/data/uboone/processed/";
  TChain mc("stv_tree"), ext("stv_tree");
  const char* mcf[]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
    "Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc",
    "Run1_rhc_new_numi_flux_rhc_pandora_ntuple","Run2_rhc_new_numi_flux_rhc_pandora_ntuple",
    "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple","Run4b_rhc_new_numi_flux_rhc_pandora_ntuple",
    "Run4c_rhc_new_numi_flux_rhc_pandora_ntuple"};
  for(auto f:mcf) mc.Add(Form("%sxsec-ana-%s.root",P,f));
  for(auto s:{"aa","ab","ac","ad","ae"}) mc.Add(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",P,s));
  ext.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
  gStyle->SetOptStat(0);
  TH1D hs("hs","",40,0,200), hb("hb","",40,0,200), he("he","",40,0,200);
  fill(mc,true,1,hs); fill(mc,true,0,hb); fill(ext,false,0,he);
  printf("backward-to-wall distance: signal=%.0f nu-bkg=%.0f EXT=%.0f\n",hs.Integral(),hb.Integral(),he.Integral());
  // discrimination scan: fraction kept for cut d > X
  printf("%-14s %8s %8s %8s\n","cut d>X[cm]","effS","effB","effE");
  for(double x:{5,10,15,20,30,50}){
    double eS=hs.Integral(hs.FindBin(x),41)/hs.Integral(), eB=hb.Integral(hb.FindBin(x),41)/hb.Integral(), eE=he.Integral(he.FindBin(x),41)/he.Integral();
    printf("d > %-9.0f %7.1f%% %7.1f%% %7.1f%%\n",x,100*eS,100*eB,100*eE);
  }
  int cSig=TColor::GetColor("#0072B2"),cBkg=TColor::GetColor("#E69F00"),cExt=TColor::GetColor("#999999");
  for(TH1D* h:{&hs,&hb,&he}) if(h->Integral()>0) h->Scale(1.0/h->Integral());
  hs.SetLineColor(cSig);hb.SetLineColor(cBkg);he.SetLineColor(cExt);
  hs.SetLineWidth(3);hb.SetLineWidth(2);he.SetLineWidth(3);hs.SetFillColorAlpha(cSig,0.15);he.SetFillColorAlpha(cExt,0.25);
  hs.SetMaximum(std::max({hs.GetMaximum(),hb.GetMaximum(),he.GetMaximum()})*1.3);
  hs.SetTitle(";muon backward distance to AV wall [cm];area-normalised");
  TCanvas c("c","",800,600); hs.Draw("hist");hb.Draw("hist same");he.Draw("hist same");
  TLegend lg(0.55,0.72,0.88,0.88);lg.SetBorderSize(0);lg.SetFillStyle(0);
  lg.AddEntry(&hs,"Signal","l");lg.AddEntry(&hb,"#nu background","l");lg.AddEntry(&he,"EXT (cosmic)","l");lg.Draw();
  c.SaveAs("unfold_output/ext_wallgap.pdf"); printf("wrote unfold_output/ext_wallgap.pdf\n");
}
