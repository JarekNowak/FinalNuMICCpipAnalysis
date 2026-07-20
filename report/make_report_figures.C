// Regenerate all RESULTS figures for the analysis note from the current pipeline:
//   custom xsec  : xsec_analyzer/unfold_output/xsec_<obs>_wiener_svd_fw.root
//   framework    : xsec_analyzer/unfold_output/ccpi_Run1_<obs>_xsec.root
//   reco spectra : unfolding_inputs.root (per-component sig/bkg/ext/dirt/data)
// Writes PNGs into report/figures/ (overwrites). Run from the repo root:
//   root -l -b -q report/make_report_figures.C
#include <vector>
#include <string>

static const char* OBS[5]   = {"pmu","ppi","costhmu","costhpi","thmupi"};
static const char* GNM[5]   = {"xsec_dsdpmu","xsec_dsdppi","xsec_dsdcosthmu","xsec_dsdcosthpi","xsec_dsdthmupi"};
static const char* XLAB[5]  = {"p_{#mu} [GeV/c]","p_{#pi} [GeV/c]","cos#theta_{#mu}","cos#theta_{#pi}","#theta_{#mu#pi} [rad]"};
static const char* DLAB[5]  = {"d#sigma/dp_{#mu}","d#sigma/dp_{#pi}","d#sigma/dcos#theta_{#mu}","d#sigma/dcos#theta_{#pi}","d#sigma/d#theta_{#mu#pi}"};
static const char* FIG      = "report/figures/";

// Framework dsigma/dx histogram (stored = sigma-per-bin; dsdx = content*0.1/width)
TH1D* fwDsdx(const char* obs) {
  TFile* ff = TFile::Open(Form("xsec_analyzer/unfold_output/ccpi_Run1_%s_xsec.root",obs));
  TDirectory* xu=(TDirectory*)ff->Get("XsecUnits"); TString sl="";
  for(auto key:*xu->GetListOfKeys()){TString nm=key->GetName(); if(nm!="Covariances"&&nm!="binnumber"){sl=nm;break;}}
  TH1D* hf=(TH1D*)xu->Get(Form("%s/%s_total",sl.Data(),sl.Data()));
  TH1D* h=(TH1D*)hf->Clone(Form("fw_%s",obs)); h->SetDirectory(0);
  for(int b=1;b<=h->GetNbinsX();b++){double w=h->GetBinWidth(b); h->SetBinContent(b,hf->GetBinContent(b)*0.1/w); h->SetBinError(b,0);}
  ff->Close(); return h;
}

void make_report_figures() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(80);

  TFile* fin = TFile::Open("unfolding_inputs.root");

  std::vector<double> sc(5), sf(5);

  for (int k=0;k<5;k++) {
    const char* o = OBS[k];
    TFile* fx = TFile::Open(Form("xsec_analyzer/unfold_output/xsec_%s_wiener_svd_fw.root",o));
    TH2D* hsm = (TH2D*)fx->Get("h_smear_sel");
    TH1D* hef = (TH1D*)fx->Get("h_efficiency");
    TH1D* htg = (TH1D*)fx->Get("h_true_gen");
    TH1D* hmeas=(TH1D*)fx->Get("h_measurement");
    TH1D* hbkg =(TH1D*)fx->Get("h_bkg_total");
    TH1D* hsig =(TH1D*)fx->Get("h_reco_sig");
    TGraphErrors* g=(TGraphErrors*)fx->Get(GNM[k]);
    if(!g){TIter it(fx->GetListOfKeys());TKey*kk;while((kk=(TKey*)it())){if(TString(kk->GetClassName())=="TGraphErrors"){g=(TGraphErrors*)fx->Get(kk->GetName());break;}}}

    // ---------- reco spectrum (stacked MC components + data) ----------
    {
      TH1D* rsig=(TH1D*)fin->Get(Form("h_reco_sig_%s",o));
      TH1D* rbkg=(TH1D*)fin->Get(Form("h_reco_bkg_%s",o));
      TH1D* rext=(TH1D*)fin->Get(Form("h_reco_ext_%s",o));
      TH1D* rdrt=(TH1D*)fin->Get(Form("h_reco_dirt_%s",o));
      TH1D* rdat=(TH1D*)fin->Get(Form("h_reco_data_%s",o));
      TCanvas c("c","",700,550);
      THStack* hs=new THStack("hs",Form(";%s;Events (3.283#times10^{20} POT)",XLAB[k]));
      rsig->SetFillColor(kAzure+1); rbkg->SetFillColor(kOrange+1);
      rext->SetFillColor(kGray+1);  rdrt->SetFillColor(kGreen+2);
      for(auto h:{rsig,rbkg,rext,rdrt}){h->SetLineColor(kBlack);h->SetLineWidth(1);}
      hs->Add(rsig);hs->Add(rbkg);hs->Add(rdrt);hs->Add(rext);
      rdat->SetMarkerStyle(20);rdat->SetMarkerSize(0.9);rdat->SetLineColor(kBlack);
      double ymax=std::max(hs->GetMaximum(),rdat->GetMaximum())*1.35;
      hs->SetMaximum(ymax); hs->Draw("hist"); rdat->Draw("E1 same");
      TLegend* l=new TLegend(0.60,0.62,0.88,0.88);
      l->AddEntry(rdat,"NuWro fake data","lep");
      l->AddEntry(rsig,"Signal (CC#pi^{+})","f");
      l->AddEntry(rbkg,"Beam bkg","f");
      l->AddEntry(rdrt,"Dirt","f");
      l->AddEntry(rext,"EXT (cosmic)","f");
      l->SetBorderSize(0);l->SetFillStyle(0);l->Draw();
      c.SaveAs(Form("%sreco_%s.png",FIG,o));
    }

    // ---------- response (smearceptance) matrix: S(reco,true) ----------
    {
      int nt=htg->GetNbinsX();
      TH2D* S=(TH2D*)hsm->Clone(Form("S_%s",o)); S->SetDirectory(0); S->Reset();
      for(int t=1;t<=nt;t++){double col=0;for(int r=1;r<=hsm->GetNbinsY();r++)col+=hsm->GetBinContent(t,r);
        if(col>0)for(int r=1;r<=hsm->GetNbinsY();r++) S->SetBinContent(t,r,hsm->GetBinContent(t,r)/col);}
      S->SetTitle(Form("Response (column-normalised);%s (true);%s (reco)",XLAB[k],XLAB[k]));
      TCanvas c("c","",680,600); c.SetRightMargin(0.15);
      S->Draw("colz"); c.SaveAs(Form("%sresp_%s.png",FIG,o));
    }

    // ---------- efficiency ----------
    {
      TCanvas c("c","",700,520);
      TH1D* he=(TH1D*)hef->Clone(); he->SetDirectory(0);
      he->SetTitle(Form("Selection efficiency;%s (true);#varepsilon",XLAB[k]));
      he->SetLineColor(kAzure+2);he->SetLineWidth(2);he->SetMarkerStyle(20);he->SetMarkerColor(kAzure+2);
      he->SetMinimum(0);he->SetMaximum(he->GetMaximum()*1.3);
      he->Draw("E1"); c.SaveAs(Form("%seff_%s.png",FIG,o));
    }

    // ---------- dsigma/dx (custom, measured phase space) ----------
    {
      TCanvas c("c","",720,560);
      g->SetMarkerStyle(20);g->SetMarkerColor(kBlack);g->SetLineColor(kBlack);g->SetLineWidth(2);
      g->SetTitle(Form("%s;%s;%s [10^{-38} cm^{2}/Ar]",DLAB[k],XLAB[k],DLAB[k]));
      double ymax=0; for(int i=0;i<g->GetN();i++){double x,y;g->GetPoint(i,x,y);ymax=std::max(ymax,y+g->GetErrorY(i));}
      ymax*=1.3;
      g->SetMinimum(0); g->SetMaximum(ymax);
      g->Draw("AP");
      TLegend* l=new TLegend(0.50,0.78,0.88,0.88);
      l->AddEntry(g,"Custom (Wiener-SVD)","lep");
      l->SetBorderSize(0);l->SetFillStyle(0);l->Draw();
      c.SaveAs(Form("%sxsec_%s.png",FIG,o));
    }

    // ---------- background subtraction diagnostic ----------
    {
      TCanvas c("c","",700,540);
      TH1D* hbsub=(TH1D*)hmeas->Clone(Form("bsub_%s",o)); hbsub->SetDirectory(0); hbsub->Add(hbkg,-1);
      hmeas->SetLineColor(kBlack);hmeas->SetMarkerStyle(20);hmeas->SetMarkerSize(0.8);
      hsig->SetLineColor(kAzure+1);hsig->SetLineWidth(2);
      hbkg->SetLineColor(kRed+1);hbkg->SetLineWidth(2);
      hbsub->SetLineColor(kGreen+2);hbsub->SetLineWidth(2);hbsub->SetMarkerStyle(24);hbsub->SetMarkerColor(kGreen+2);
      hmeas->SetTitle(Form("Background subtraction;%s;Events",XLAB[k]));
      hmeas->SetMinimum(0);hmeas->SetMaximum(hmeas->GetMaximum()*1.35);
      hmeas->Draw("E1");hsig->Draw("hist same");hbkg->Draw("hist same");hbsub->Draw("E1 same");
      TLegend* l=new TLegend(0.55,0.64,0.88,0.88);
      l->AddEntry(hmeas,"Data (+EXT)","lep");l->AddEntry(hsig,"Signal MC","l");
      l->AddEntry(hbkg,"Total bkg","l");l->AddEntry(hbsub,"Bkg-subtracted","lep");
      l->SetBorderSize(0);l->SetFillStyle(0);l->Draw();
      c.SaveAs(Form("%sbkgsub_%s.png",FIG,o));
    }

    // integrated
    double s=0; for(int i=0;i<g->GetN();i++){double x,y;g->GetPoint(i,x,y);s+=y*2*g->GetErrorX(i);} sc[k]=s;
    fx->Close();
  }

  // ---------- integrated cross section (custom, measured phase space) ----------
  {
    TCanvas c("c","",760,560);
    TH1D* hc=new TH1D("hc",";;#sigma_{int} [10^{-38} cm^{2}/Ar]",5,0,5);
    const char* lab[5]={"p_{#mu}","p_{#pi}","cos#theta_{#mu}","cos#theta_{#pi}","#theta_{#mu#pi}"};
    for(int k=0;k<5;k++){hc->SetBinContent(k+1,sc[k]);hc->GetXaxis()->SetBinLabel(k+1,lab[k]);}
    hc->SetBarWidth(0.6);hc->SetBarOffset(0.2);hc->SetFillColor(kAzure+1);hc->SetLineColor(kBlack);
    hc->SetMaximum(*std::max_element(sc.begin(),sc.end())*1.25);
    hc->SetMinimum(0);hc->GetXaxis()->SetLabelSize(0.05);hc->SetStats(0);
    hc->Draw("bar");
    c.SaveAs(Form("%sxsec_integrated.png",FIG));
  }

  fin->Close();
  printf("Wrote report figures to %s\n",FIG);
}
