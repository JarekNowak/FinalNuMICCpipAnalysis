// sb_plots_data.C -- the four control regions drawn with REAL beam-on data from the
// signal-region-stripped skim (/data/uboone/processed/beamon_skim).
//
// Differences from sb_plots.C (pseudo-data), all of them mandatory for real data:
//   (1) the data are the beam-on skim, which contains every control-region event and no
//       signal-region event, so control-region distributions are complete;
//   (2) the Poisson throw of EXT+dirt that sb_plots.C ADDS to pseudo-data is REMOVED --
//       real data already contains cosmics and dirt, which here belong to the prediction;
//   (3) the prediction is scaled to ONLY the run periods for which beam-on data exists.
//       No Run-2 beam-on sample is staged, so Run 2 is dropped from the MC, the beam-off
//       and the dirt on both sides. Retained: FHC runs 1,4,5 = 7.589/8.857 = 85.7% of the
//       FHC exposure; RHC runs 1,3,4 = 8.491/11.082 = 76.6% of the RHC exposure.
//   usage: root -l -b -q 'macros/sb_plots_data.C("fhc")'
#include "sb_guard.h"
#include <vector>
#include <string>
struct Src { std::string file; double scale; };
void sb_plots_data(const char* mode="fhc", bool compact=false){
  gStyle->SetOptStat(0);
  // sb_pi0/ (2026-09-06, current selection). sb/ (2026-08-11) predates the beam-frame correction: its cos(theta_mu)
  // is about detector z while the beam-on skim is beam-frame (the 2026-09-09 muon-angle "shape failure").
  const char* P="/data/uboone/processed/sb_pi0/";       // MC with the region flags
  const char* B="/data/uboone/processed/beamon_skim/";  // beam-on control-region skim
  const char* E="/data/uboone/processed/ext_perrun/xsec-ana-";
  const char* R1="neutrinoselection_filt_run1_beamoff.root", *R3B="neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root","numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5="numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCCX=0.98;
  std::vector<Src> mc, ext; std::vector<std::string> data; double sc_dirt; double pot_frac; TString runs;
  if(std::string(mode)=="fhc"){
    // Run 2 dropped: no beam-on file.  Scales are those of the exposure table.
    mc.push_back({std::string(P)+"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",0.14101});
    mc.push_back({std::string(P)+"xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",0.07323});
    mc.push_back({std::string(P)+"xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root",0.11560});
    ext.push_back({std::string(E)+R1,OCCX*9846635./G1});
    for(auto f:R4) ext.push_back({std::string(E)+f,OCCX*4131149./G4});
    ext.push_back({std::string(E)+R5,OCCX*5154196./G5});
    for(auto r:{"run1","run4c","run4d","run5"}) data.push_back(std::string(B)+"xsec-ana-beamon_fhc_"+r+".root");
    pot_frac=(3.283+2.075+2.231)/8.857; sc_dirt=0.092402*0.65*pot_frac; runs="Runs 1, 4, 5 (no Run-2 beam-on)";
  } else {
    mc.push_back({std::string(P)+"xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root",0.06728});
    for(auto s:{"Run4a_rhc","Run4b_rhc","Run4c_rhc"}) mc.push_back({std::string(P)+"xsec-ana-"+std::string(s)+"_new_numi_flux_rhc_pandora_ntuple.root",0.08847});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066});
    ext.push_back({std::string(E)+R1,OCCX*1458253./G1});
    ext.push_back({std::string(E)+R3B,OCCX*10349610./G3});
    for(auto f:R4) ext.push_back({std::string(E)+f,OCCX*6304167./G4});
    for(auto r:{"run1","run3b","run4a","run4b"}) data.push_back(std::string(B)+"xsec-ana-beamon_rhc_"+r+".root");
    pot_frac=(0.6053+5.003+2.883)/11.082; sc_dirt=0.071666*0.65*pot_frac; runs="Runs 1, 3, 4 (no Run-2 beam-on)";
  }
  for(auto&d:data){ if(gSystem->AccessPathName(d.c_str())){ printf("MISSING %s -- run not yet processed, aborting\n",d.c_str()); return; } }
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  struct V { const char* reg; const char* name; const char* br; int nb; double lo,hi; const char* ttl; };
  std::vector<V> vars={{"sb_cc0pi","cc0pi_pmu","CC1mu1piXp_candidate_muon_mom_reco",7,0.15,1.75,"p_{#mu} [GeV/c]"},{"sb_cc0pi","cc0pi_costh","CC1mu1piXp_candidate_muon_costh_reco",5,-1,1,"cos#theta_{#mu}"},
    {"sb_multipi","multipi_pmu","CC1mu1piXp_candidate_muon_mom_reco",7,0.15,1.75,"p_{#mu} [GeV/c]"},{"sb_multipi","multipi_costh","CC1mu1piXp_candidate_muon_costh_reco",5,-1,1,"cos#theta_{#mu}"},{"sb_multipi","multipi_oa","CC1mu1piXp_mu_pi_opening_angle",6,0,3.15,"#theta_{#mu#pi} [rad]"},
    {"sb_pi0","pi0_pmu","CC1mu1piXp_candidate_muon_mom_reco",7,0.15,1.75,"p_{#mu} [GeV/c]"},{"sb_pi0","pi0_costh","CC1mu1piXp_candidate_muon_costh_reco",5,-1,1,"cos#theta_{#mu}"},{"sb_pi0","pi0_shr","shr_energy_cali",10,0,1.0,"leading shower energy [GeV]"},
    {"sb_cosmic","cosmic_pmu","CC1mu1piXp_candidate_muon_mom_reco",7,0.15,1.75,"p_{#mu} [GeV/c]"},{"sb_cosmic","cosmic_costh","CC1mu1piXp_candidate_muon_costh_reco",5,-1,1,"cos#theta_{#mu}"},{"sb_cosmic","cosmic_oa","CC1mu1piXp_mu_pi_opening_angle",4,2.6,3.15,"#theta_{#mu#pi} [rad]"}};
  if(compact){ std::vector<V> v4; for(auto&v:vars) if(!strcmp(v.name,"cc0pi_pmu")||!strcmp(v.name,"multipi_pmu")||!strcmp(v.name,"pi0_shr")||!strcmp(v.name,"cosmic_oa")) v4.push_back(v); vars=v4; }
  double pmu_edges[8]={0.15,0.35,0.55,0.75,0.95,1.25,1.75,3.0}, cth_edges[6]={-1,0.45,0.65,0.8,0.9,1};
  TCanvas c("c","",compact?1400:1800,compact?900:1200); if(compact) c.Divide(2,2); else c.Divide(4,3); int pad=0;
  std::vector<TObject*> keep;
  printf("\n==== %s control regions, BEAM-ON DATA (%s; %.1f%% of the mode exposure) ====\n", mode, runs.Data(), 100*pot_frac);
  printf("%-14s %10s %10s %8s %8s %8s %8s   %8s\n","region/var","data","pred","signal","nu bkg","beam-off","dirt","data/pred");
  for(auto& v: vars){
    TString F=TString("CC1mu1piXp_")+v.reg;
    auto mk=[&](const char* n)->TH1D*{ TH1D* h; if(TString(v.br).Contains("mom_reco")) h=new TH1D(n,"",7,pmu_edges); else if(TString(v.br).Contains("costh")) h=new TH1D(n,"",5,cth_edges); else h=new TH1D(n,"",v.nb,v.lo,v.hi); h->SetDirectory(0); keep.push_back(h); return h; };
    TH1D *hsig=mk(Form("sig_%s",v.name)), *hbkg=mk(Form("bkg_%s",v.name)), *hext=mk(Form("ext_%s",v.name)), *hdirt=mk(Form("dirt_%s",v.name)), *hdat=mk(Form("dat_%s",v.name));
    auto fill=[&](TH1D* h,const std::string& file,double scale,const TString& cut,bool w){ TChain ch("stv_tree"); ch.Add(file.c_str());
      TH1D* tmp=(TH1D*)h->Clone("sbtmp"); tmp->Reset(); tmp->SetDirectory(gROOT);
      ch.Draw(Form("min(max(%s,%g),%g)>>sbtmp",v.br,h->GetXaxis()->GetXmin()+1e-6,h->GetXaxis()->GetXmax()-1e-6),(w?TString(CVW):TString("1"))+"*("+cut+")","goff");
      h->Add(tmp,scale); tmp->SetDirectory(0); delete tmp; };
    for(auto&x:mc){ fill(hsig,x.file,x.scale,F+" && CC1mu1piXp_MC_Signal",true); fill(hbkg,x.file,x.scale,F+" && !CC1mu1piXp_MC_Signal",true); }
    for(auto&x:ext) fill(hext,x.file,x.scale,F,false);
    fill(hdirt,std::string(P)+"xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",sc_dirt,F,true);
    for(auto&d:data){ sb_guard_data(d); fill(hdat,d,1.0,F,false); }   // NO Poisson EXT/dirt added: real data
    for(int b=1;b<=hdat->GetNbinsX();b++) hdat->SetBinError(b,sqrt(std::max(0.,hdat->GetBinContent(b))));
    double D=hdat->Integral(0,hdat->GetNbinsX()+1), S=hsig->Integral(0,hsig->GetNbinsX()+1), Bk=hbkg->Integral(0,hbkg->GetNbinsX()+1),
           X=hext->Integral(0,hext->GetNbinsX()+1), Dt=hdirt->Integral(0,hdirt->GetNbinsX()+1), Pr=S+Bk+X+Dt;
    printf("%-14s %10.0f %10.1f %8.1f %8.1f %8.1f %8.1f   %8.3f\n", v.name, D, Pr, S, Bk, X, Dt, Pr>0?D/Pr:0.);
    c.cd(++pad); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    THStack* st=new THStack(Form("st_%s",v.name),""); keep.push_back(st);
    hdirt->SetFillColor(TColor::GetColor("#009E73")); hext->SetFillColor(kGray+1); hbkg->SetFillColor(TColor::GetColor("#E69F00")); hsig->SetFillColor(TColor::GetColor("#0072B2"));
    for(auto h:{hdirt,hext,hbkg,hsig}){ h->SetLineColor(kBlack); h->SetLineWidth(1); st->Add(h); }
    double ymax=std::max(st->GetMaximum(),hdat->GetMaximum()+hdat->GetBinError(hdat->GetMaximumBin()));
    st->SetMaximum(1.45*ymax); st->Draw("hist"); st->GetXaxis()->SetTitle(v.ttl); st->GetYaxis()->SetTitle("events"); st->GetXaxis()->SetTitleSize(0.05); st->GetYaxis()->SetTitleSize(0.05);
    hdat->SetMarkerStyle(20); hdat->SetMarkerSize(0.8); hdat->Draw("E1 same");
    TLegend* lg=new TLegend(0.42,0.60,0.89,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036); keep.push_back(lg);
    lg->SetHeader(Form("%s  %s  [%s]", mode, v.reg+3, runs.Data()));
    lg->AddEntry(hdat,"beam-on data","lep"); lg->AddEntry(hsig,"signal (contamination)","f"); lg->AddEntry(hbkg,"#nu background","f"); lg->AddEntry(hext,"beam-off","f"); lg->AddEntry(hdirt,"dirt","f"); lg->Draw();
  }
  c.SaveAs(Form("../report/figures/sideband_data_%s%s.pdf",mode,compact?"_compact":"")); printf("wrote sideband_data_%s%s.pdf\n",mode,compact?"_compact":"");
}
