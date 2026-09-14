// eff_turnon.C -- selection efficiency across the measured-phase-space thresholds (MC only).
// Signal WITHOUT the kinematic thresholds (the sig_ component flags); each turn-on applies the
// other two phase-space cuts, so it shows where its own threshold sits. FHC, per-run POT-scaled,
// full CV weight. One pass per file (the ntuples are ~10 GB each).
//     root -l -b -q macros/eff_turnon.C
#include <vector>
#include <string>
void eff_turnon(){
  gStyle->SetOptStat(0);
  const char* P="/data/uboone/processed/xsec-ana-";
  std::vector<std::pair<std::string,double>> mc={{"Run1_fhc_new_numi_flux_fhc_pandora_ntuple",0.14101},{"Run2_fhc_new_numi_flux_fhc_pandora_ntuple",0.05085},
    {"Run4_fhc_new_numi_flux_fhc_pandora_ntuple",0.07323},{"reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc",0.11560}};
  const char* TF="CC1mu1piXp_sig_truevertex_in_fv && CC1mu1piXp_sig_ccnc && CC1mu1piXp_sig_is_numu && CC1mu1piXp_sig_one_muon_above_thresh && CC1mu1piXp_sig_one_charged_pion && CC1mu1piXp_sig_no_pions && CC1mu1piXp_sig_no_heavy_mesons";
  const char* W="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  const char* PMU="CC1mu1piXp_candidate_muon_mom_true", *PPI="CC1mu1piXp_candidate_pion_mom_true", *OA="CC1mu1piXp_true_mu_pi_opening_angle";
  struct V { const char* name; int nb; double lo,hi,thr; bool upper; const char* ttl; };
  std::vector<V> vars={{"pmu",40,0.0,1.0,0.15,false,"true p_{#mu} [GeV/c]"},{"ppi",24,0.0,0.6,0.175,false,"true p_{#pi} [GeV/c]"},{"thmupi",32,0.0,3.2,2.6,true,"true #theta_{#mu#pi} [rad]"}};
  std::vector<TH1D*> den, num;
  for(auto& v: vars){ den.push_back(new TH1D(Form("den_%s",v.name),"",v.nb,v.lo,v.hi)); num.push_back(new TH1D(Form("num_%s",v.name),"",v.nb,v.lo,v.hi)); den.back()->Sumw2(); num.back()->Sumw2(); }
  double sig_in=0, sel_in=0;
  for(auto& m: mc){
    TFile f((std::string(P)+m.first+".root").c_str()); TTree* t=(TTree*)f.Get("stv_tree");
    t->SetBranchStatus("*",0);
    for(auto b:{"CC1mu1piXp_Selected","CC1mu1piXp_MC_Signal","CC1mu1piXp_sig_truevertex_in_fv","CC1mu1piXp_sig_ccnc","CC1mu1piXp_sig_is_numu","CC1mu1piXp_sig_one_muon_above_thresh",
      "CC1mu1piXp_sig_one_charged_pion","CC1mu1piXp_sig_no_pions","CC1mu1piXp_sig_no_heavy_mesons",PMU,PPI,OA,"tuned_cv_weight","ppfx_cv_weight","normalisation_weight"}) t->SetBranchStatus(b,1);
    TTreeFormula fTF("tf",TF,t), fSel("sel","CC1mu1piXp_Selected",t), fSig("sig","CC1mu1piXp_MC_Signal",t), fW("w",W,t), fMu("mu",PMU,t), fPi("pi",PPI,t), fOA("oa",OA,t);
    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;i++){
      Long64_t le=t->LoadTree(i); if(le<0) break;
      for(auto fm:{&fTF,&fSel,&fSig,&fW,&fMu,&fPi,&fOA}) fm->UpdateFormulaLeaves();
      bool sig=fSig.EvalInstance()!=0, tf=fTF.EvalInstance()!=0; if(!sig && !tf) continue;
      double w=m.second*fW.EvalInstance(); bool sel=fSel.EvalInstance()!=0;
      if(sig){ sig_in+=w; if(sel) sel_in+=w; }
      if(!tf) continue;
      double mu=fMu.EvalInstance(), pi=fPi.EvalInstance(), oa=fOA.EvalInstance();
      double x[3]={mu,pi,oa}; bool others[3]={pi>0.175&&oa<2.6, mu>0.15&&oa<2.6, mu>0.15&&pi>0.175};
      for(int k=0;k<3;k++){ if(!others[k]) continue; double xc=std::min(std::max(x[k],vars[k].lo+1e-6),vars[k].hi-1e-6);
        den[k]->Fill(xc,w); if(sel) num[k]->Fill(xc,w); }
    }
    printf("done %s (%lld entries)\n",m.first.c_str(),N); fflush(stdout);
  }
  printf("closure: efficiency in the measured phase space = %.2f%% (signal %.1f, selected %.1f)\n",100*sel_in/sig_in,sig_in,sel_in);
  for(size_t k=0;k<vars.size();k++){ auto& v=vars[k];
    TCanvas c(Form("c_%s",v.name),"",700,500); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.1);
    TH1D* e=(TH1D*)num[k]->Clone(Form("eff_%s",v.name)); e->Divide(num[k],den[k],1,1,"B");
    printf("\n== %s (threshold %.3f) ==\n",v.name,v.thr);
    for(int b=1;b<=e->GetNbinsX();b++) printf("  [%.3f,%.3f] eff %5.1f%%  gen %.1f\n",e->GetBinLowEdge(b),e->GetBinLowEdge(b+1),100*e->GetBinContent(b),den[k]->GetBinContent(b));
    double ymax=1.3*e->GetMaximum(); e->SetMaximum(ymax); e->SetMinimum(0);
    e->GetXaxis()->SetTitle(v.ttl); e->GetYaxis()->SetTitle("selection efficiency"); e->GetXaxis()->SetTitleSize(0.055); e->GetYaxis()->SetTitleSize(0.055);
    e->GetXaxis()->SetLabelSize(0.045); e->GetYaxis()->SetLabelSize(0.045);
    e->SetLineColor(TColor::GetColor("#0072B2")); e->SetLineWidth(2); e->SetMarkerStyle(20); e->SetMarkerSize(0.7); e->SetMarkerColor(TColor::GetColor("#0072B2"));
    e->Draw("E1");
    TBox box(v.upper?v.thr:v.lo,0,v.upper?v.hi:v.thr,ymax); box.SetFillColorAlpha(kGray+1,0.35); box.Draw();
    TLine ln(v.thr,0,v.thr,ymax); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw();
    e->Draw("E1 same");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045); tx.DrawLatex(0.13,0.93,"MicroBooNE NuMI FHC simulation");
    c.SaveAs(Form("../report/figures/eff_turnon_%s_fhc.pdf",v.name));
  }
}
