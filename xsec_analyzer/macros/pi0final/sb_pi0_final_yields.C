// sb_pi0_final_yields.C -- the pi0 control region with and without the final multiplicity
// cut (pre-ratification comparison, 2026-09-06), on the sb_pi0/ reprocess (MC, dirt, per-run
// beam-off; no beam-on file). For each definition: nu-MC yield by class, target purity
// (CC+NC pi0, of the nu-MC), signal contamination, beam-off and dirt (full-stack shares),
// overlap with multi-pi, and the frozen-distribution shapes (leading shower energy) plus the
// multiplicity distributions that the cut acts on. Scales as in sb_protocol.C.
//   root -l -b -q 'macros/pi0final/sb_pi0_final_yields.C("fhc")'
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
struct Src { std::string file; double scale; };
void sb_pi0_final_yields(const char* mode="fhc"){
  const char* P="/data/uboone/processed/sb_pi0/"; std::vector<Src> mc, ext; double sc_dirt;
  if(std::string(mode)=="fhc"){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); sc_dirt=0.092402*0.65; }
  else { const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); sc_dirt=0.071666*0.65; }
  // per-run beam-off, scaled by gates as in sb_protocol.C add_ext (same numbers)
  const char* E="/data/uboone/processed/sb_pi0/xsec-ana-";
  const char* R1="neutrinoselection_filt_run1_beamoff.root", *R3B="neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root","numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5="numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCCX=0.98;
  if(std::string(mode)=="fhc"){ ext.push_back({std::string(E)+R1,OCCX*9846635./G1}); ext.push_back({std::string(E)+R1,OCCX*3535129./G1}); for(auto f:R4) ext.push_back({std::string(E)+f,OCCX*4131149./G4}); ext.push_back({std::string(E)+R5,OCCX*5154196./G5}); }
  else { ext.push_back({std::string(E)+R1,OCCX*1458253./G1}); ext.push_back({std::string(E)+R1,OCCX*5422907./G1}); ext.push_back({std::string(E)+R3B,OCCX*10349610./G3}); for(auto f:R4) ext.push_back({std::string(E)+f,OCCX*6304167./G4}); }
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* NPI="CC1mu1piXp_mc_n_threshold_pionpm", *NPI0="CC1mu1piXp_mc_n_threshold_pion0";
  const int NC=5; TString cls[NC]; const char* cln[NC]={"signal","CC+NC pi0 (target)","numuCC 0pi","other nu (non-target)","all nu-MC"};
  cls[0]="CC1mu1piXp_EventCategory==0";
  cls[1]=TString("(CC1mu1piXp_EventCategory==3||CC1mu1piXp_EventCategory==5) && ")+NPI0+">=1";
  cls[2]=TString("CC1mu1piXp_EventCategory==3 && ")+NPI+"==0 && "+NPI0+"==0";
  cls[3]=TString("!(")+cls[0]+") && !("+cls[1]+")";
  cls[4]="1";
  const int ND=2; const char* def[ND]={"CC1mu1piXp_sb_pi0","CC1mu1piXp_sb_pi0_final"}; const char* dn[ND]={"without final cut (frozen)","with final cut"};
  double N[ND][NC]={{0}}, NSR[NC]={0}, Nov[ND]={0}, Next[ND]={0}, Ndirt[ND]={0};
  const int NV=3; const char* var[NV]={"CC1mu1piXp_shr_energy_cali","CC1mu1piXp_sb_nprimtrk","CC1mu1piXp_sb_nnonproton"}; const char* vn[NV]={"E_{shr} [GeV]","N_{primary tracks}","N_{non-proton}"}; int nb[NV]={10,8,8}; double lo[NV]={0,-0.5,-0.5}, hi[NV]={1.0,7.5,7.5};
  TH1D* H[ND][NV]; for(int d=0;d<ND;d++) for(int v=0;v<NV;v++){ H[d][v]=new TH1D(Form("h_%d_%d",d,v),"",nb[v],lo[v],hi[v]); H[d][v]->Sumw2(); }
  // the shower-energy branch name: use the frozen plotting-list variable if the selection stores it under another name
  for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str());
    if(&x==&mc[0]){ if(!c.GetBranch(var[0])) { var[0]="shr_energy_cali"; } }
    for(int d=0;d<ND;d++){ for(int k=0;k<NC;k++){ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",TString(CVW)+"*(("+cls[k]+") && "+def[d]+")","goff"); N[d][k]+=h.Integral(0,2)*x.scale; }
      { TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",TString(CVW)+"*("+def[d]+" && CC1mu1piXp_sb_multipi)","goff"); Nov[d]+=h.Integral(0,2)*x.scale; }
      for(int v=0;v<NV;v++){ TH1D h("hv","",nb[v],lo[v],hi[v]); c.Draw(TString(var[v])+">>hv",TString(CVW)+"*("+def[d]+")","goff"); H[d][v]->Add(&h,x.scale); } }
    for(int k=0;k<NC;k++){ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",TString(CVW)+"*(("+cls[k]+") && CC1mu1piXp_Selected)","goff"); NSR[k]+=h.Integral(0,2)*x.scale; } }
  for(auto&x:ext){ TChain c("stv_tree"); c.Add(x.file.c_str()); for(int d=0;d<ND;d++){ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",def[d],"goff"); Next[d]+=h.Integral(0,2)*x.scale; } }
  { TChain c("stv_tree"); c.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P)); for(int d=0;d<ND;d++){ TH1D h("h","",1,-0.5,1.5); c.Draw("0.5>>h",TString(CVW)+"*("+def[d]+")","goff"); Ndirt[d]+=h.Integral(0,2)*sc_dirt; } }
  printf("\n==== %s: pi0 region, two definitions (CV-weighted, POT/gate-scaled)\n",mode);
  printf("%-34s %18s %18s\n","", dn[0], dn[1]);
  for(int k=0;k<NC;k++) printf("%-34s %18.1f %18.1f\n",cln[k],N[0][k],N[1][k]);
  for(int d=0;d<ND;d++){ double stack=N[d][4]+Next[d]+Ndirt[d]; printf("%s: beam-off %.1f  dirt %.1f  full stack %.1f\n",dn[d],Next[d],Ndirt[d],stack); }
  printf("%-34s %17.1f%% %17.1f%%\n","target purity (of nu-MC)",100*N[0][1]/N[0][4],100*N[1][1]/N[1][4]);
  printf("%-34s %17.1f%% %17.1f%%\n","target share (of full stack)",100*N[0][1]/(N[0][4]+Next[0]+Ndirt[0]),100*N[1][1]/(N[1][4]+Next[1]+Ndirt[1]));
  printf("%-34s %17.1f%% %17.1f%%\n","signal contamination (of nu-MC)",100*N[0][0]/N[0][4],100*N[1][0]/N[1][4]);
  printf("%-34s %17.1f%% %17.1f%%\n","beam-off share (of full stack)",100*Next[0]/(N[0][4]+Next[0]+Ndirt[0]),100*Next[1]/(N[1][4]+Next[1]+Ndirt[1]));
  printf("%-34s %18.1f %18.1f\n","overlap with multi-pi (nu-MC events)",Nov[0],Nov[1]);
  printf("%-34s %18.4f %18.4f\n","T(target) = N_SR/N_CR",NSR[1]/N[0][1],NSR[1]/N[1][1]);
  printf("%-34s %18.3f %18.3f\n","leverage non-target/target",N[0][3]/N[0][1],N[1][3]/N[1][1]);
  printf("%-34s %18.3f %18.3f\n","leverage signal/target",N[0][0]/N[0][1],N[1][0]/N[1][1]);
  printf("kept fraction with the final cut: nu-MC %.3f, target %.3f, signal %.3f, beam-off %.3f\n",N[1][4]/N[0][4],N[1][1]/N[0][1],N[1][0]/N[0][0],Next[1]/Next[0]);
  // shape comparison: the frozen pi0 distribution and the multiplicities, normalised to unit area
  gStyle->SetOptStat(0); TCanvas cv("cv","",1500,500); cv.Divide(3,1);
  for(int v=0;v<NV;v++){ cv.cd(v+1); double m=0; for(int d=0;d<ND;d++){ H[d][v]->SetLineColor(d?kRed+1:kBlue+1); H[d][v]->SetLineWidth(2); H[d][v]->GetXaxis()->SetTitle(vn[v]); H[d][v]->GetYaxis()->SetTitle("events (full exposure)"); m=std::max(m,H[d][v]->GetMaximum()); }
    H[0][v]->SetMaximum(1.3*m); H[0][v]->Draw("hist"); H[1][v]->Draw("hist same");
    if(v==0){ TLegend* L=new TLegend(0.45,0.7,0.88,0.88); L->AddEntry(H[0][v],"#pi^{0} region, frozen (no final cut)","l"); L->AddEntry(H[1][v],"with final multiplicity cut","l"); L->SetBorderSize(0); L->Draw(); }
    double chi2=0; int nb2=0; for(int b=1;b<=nb[v];b++){ double a=H[0][v]->GetBinContent(b)/H[0][v]->Integral(), bb=H[1][v]->GetBinContent(b)/H[1][v]->Integral(); double ea=H[0][v]->GetBinError(b)/H[0][v]->Integral(), eb=H[1][v]->GetBinError(b)/H[1][v]->Integral(); if(ea+eb>0){ chi2+=pow(a-bb,2)/(ea*ea+eb*eb); nb2++; } }
    printf("shape %s: unit-area chi2/nbins = %.2f/%d (MC-stat only; the two samples are nested, so this is indicative)\n",var[v],chi2,nb2); }
  cv.SaveAs(Form("../report/figures/sb_pi0_final_%s.pdf",mode));
}
