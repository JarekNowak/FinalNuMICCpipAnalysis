void mares_band(){
  TFile fn("genie_mares_nom.root"),fm("genie_mares_minus.root"),fp("genie_mares_plus.root");
  TH1D*hn=(TH1D*)fn.Get("costhmu"),*hm=(TH1D*)fm.Get("costhmu"),*hp=(TH1D*)fp.Get("costhmu");
  int nb=hn->GetNbinsX();
  TGraphAsymmErrors* band=new TGraphAsymmErrors(nb);
  for(int b=1;b<=nb;b++){ double x=hn->GetBinCenter(b),y=hn->GetBinContent(b);
    double lo=TMath::Min(hm->GetBinContent(b),hp->GetBinContent(b));
    double hi=TMath::Max(hm->GetBinContent(b),hp->GetBinContent(b));
    band->SetPoint(b-1,x,y); band->SetPointError(b-1,hn->GetBinWidth(b)/2,hn->GetBinWidth(b)/2,y-lo,hi-y);}
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","",900,650); c->SetLeftMargin(0.13); c->SetBottomMargin(0.12);
  hn->SetTitle(""); hn->GetXaxis()->SetTitle("cos#theta_{#mu}"); hn->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/Ar]");
  hn->GetYaxis()->SetTitleOffset(1.1); hn->SetMaximum(hp->GetMaximum()*1.15); hn->SetMinimum(0);
  hn->SetLineColor(kBlack); hn->SetLineWidth(3); hn->Draw("hist");
  band->SetFillColorAlpha(kOrange+1,0.45); band->SetLineWidth(0); band->Draw("3 same");
  hn->Draw("hist same");
  hm->SetLineColor(kBlue+1); hm->SetLineStyle(2); hm->SetLineWidth(2); hm->Draw("hist same");
  hp->SetLineColor(kRed+1); hp->SetLineStyle(2); hp->SetLineWidth(2); hp->Draw("hist same");
  TLegend* lg=new TLegend(0.16,0.63,0.55,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
  lg->AddEntry(hn,"GENIE nominal","l"); lg->AddEntry(band,"M_{A}^{RES} #pm1#sigma band","f");
  lg->AddEntry(hp,"M_{A}^{RES} +1#sigma (+32.5%)","l"); lg->AddEntry(hm,"M_{A}^{RES} #minus1#sigma (#minus5.9%)","l");
  lg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.038); t.DrawLatex(0.55,0.83,"CC1#mu1#pi Xp, NuMI FHC");
  c->SaveAs("/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/unfold_output/mares_band_costhmu.pdf");
  printf("wrote mares_band_costhmu.pdf\n");
}
