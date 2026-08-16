void sat(){
  gStyle->SetOptStat(0);
  TFile f("/data/uboone/processed/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
  TTree* t=(TTree*)f.Get("stv_tree");
  const char* COR="sqrt(pow(sqrt(pow(CC1mu1piXp_candidate_pion_mom_reco,2)+0.011164)-0.10566+0.13957,2)-0.019480)";
  const char* RAW="CC1mu1piXp_candidate_pion_mom_reco";
  TString S="CC1mu1piXp_Selected && CC1mu1piXp_MC_Signal && CC1mu1piXp_candidate_pion_mom_true>0.175";
  const int N=8; double e[N+1]={0.175,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.20};
  double x[N],yc[N],yr[N],ex[N],eyc[N];
  for(int i=0;i<N;i++){
    TString c=Form("%s && CC1mu1piXp_candidate_pion_mom_true>%g && CC1mu1piXp_candidate_pion_mom_true<%g",S.Data(),e[i],e[i+1]);
    TH1D h1("h1","",1,-1e9,1e9),h2("h2","",1,-1e9,1e9),h3("h3","",1,-1e9,1e9);
    t->Draw("0>>h1",Form("(%s)*(%s)",c.Data(),COR),"goff");
    t->Draw("0>>h2",Form("(%s)*(CC1mu1piXp_candidate_pion_mom_true)",c.Data()),"goff");
    t->Draw("0>>h3",Form("(%s)*(%s)",c.Data(),RAW),"goff");
    Long64_t n=t->GetEntries(c);
    x[i]=h2.Integral(0,2)/n; yc[i]=h1.Integral(0,2)/n; yr[i]=h3.Integral(0,2)/n;
    ex[i]=0; eyc[i]=0;
  }
  TCanvas c("c","",800,600);
  c.SetLeftMargin(0.13); c.SetBottomMargin(0.13);
  TH2D fr("fr",";true p_{#pi} [GeV/c];reconstructed p_{#pi} [GeV/c]",10,0.15,1.05,10,0.15,1.05);
  fr.GetXaxis()->SetTitleSize(0.045); fr.GetYaxis()->SetTitleSize(0.045);
  fr.Draw();
  TLine* id=new TLine(0.15,0.15,1.05,1.05); id->SetLineStyle(2); id->SetLineColor(kGray+2); id->SetLineWidth(2); id->Draw();
  TGraph* gc=new TGraph(N,x,yc); gc->SetMarkerStyle(20); gc->SetMarkerSize(1.3);
  gc->SetMarkerColor(TColor::GetColor("#D55E00")); gc->SetLineColor(TColor::GetColor("#D55E00")); gc->SetLineWidth(3);
  TGraph* gr=new TGraph(N,x,yr); gr->SetMarkerStyle(24); gr->SetMarkerSize(1.2);
  gr->SetMarkerColor(TColor::GetColor("#0072B2")); gr->SetLineColor(TColor::GetColor("#0072B2")); gr->SetLineWidth(2); gr->SetLineStyle(7);
  gr->Draw("PL same"); gc->Draw("PL same");
  TLine* s1=new TLine(0.15,0.34,1.05,0.34); s1->SetLineStyle(3); s1->SetLineColor(kRed+1); s1->SetLineWidth(2); s1->Draw();
  TLatex tx; tx.SetTextSize(0.033); tx.SetTextColor(kRed+1);
  tx.DrawLatex(0.62,0.355,"saturation #approx 0.34 GeV/c");
  tx.SetTextColor(kGray+3); tx.SetTextAngle(38); tx.DrawLatex(0.60,0.64,"ideal (reco = true)");
  TLegend* lg=new TLegend(0.16,0.66,0.64,0.82);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
  lg->AddEntry(gc,"analysis estimator (#mu#rightarrow#pi mass corrected)","pl");
  lg->AddEntry(gr,"raw branch (muon hypothesis)","pl");
  lg->Draw();
  tx.SetTextAngle(0); tx.SetTextColor(kBlack); tx.SetTextSize(0.030);
  tx.SetNDC(); tx.SetTextSize(0.034); tx.DrawLatex(0.13,0.94,"MicroBooNE NuMI simulation: selected CC1#mu1#pi^{#pm} signal"); tx.SetNDC(kFALSE);
  c.SaveAs("/data/uboone/processed/validate/ppi_saturation.png"); c.SaveAs("../report/figures/ppi_saturation.png");
}
