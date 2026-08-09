// systbreak_fig.C — colorblind-safe systematic-breakdown figure per configuration,
// from the [SYSTDUMP] output of UnfolderNuMI (../logs/systdump/<cfg>_<obs>.dump).
// One figure per config: bin-averaged fractional uncertainty (%) by source, grouped
// by observable. Okabe-Ito palette. Replaces the orphaned Run-1 fw_systbreak_*.
//   usage: root -l -b -q 'macros/systbreak_fig.C("fhc5")'  // or "rhcfull","comb"
#include <map>
#include <string>
#include <fstream>
void systbreak_fig(const char* cfg="fhc5") {
  const char* D="../logs/systdump/";
  const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* obsT[5]={"p_{#mu}","p_{#pi}","cos#theta_{#mu}","cos#theta_{#pi}","#theta_{#mu#pi}"};
  // sources (key in dump -> label, Okabe-Ito colour)
  struct Src{const char* key; const char* lab; const char* hex;};
  std::vector<Src> src = {
    {"total","Total","#000000"}, {"flux_total","Flux","#0072B2"},
    {"detVar_total","Detector","#009E73"}, {"xsec_total","Cross section","#D55E00"},
    {"reint","Reinteraction","#CC79A7"}, {"DataStats","Data stat","#E69F00"} };
  auto val=[&](const char* o,const char* key)->double{
    std::ifstream f(Form("%s%s_%s.dump",D,cfg,o)); std::string tag,name; double v;
    while(f>>tag>>name>>v){ if(tag=="[SYSTDUMP]" && name==key) return v; } return 0.; };
  gStyle->SetOptStat(0);
  TCanvas c("c","",1000,600); c.SetBottomMargin(0.12);
  TH1D* frame=new TH1D("fr",Form("Systematic breakdown (%s);;bin-averaged fractional uncertainty [%%]",cfg),5,0,5);
  for(int o=0;o<5;o++) frame->GetXaxis()->SetBinLabel(o+1,obsT[o]);
  double ymax=0; for(int o=0;o<5;o++) ymax=std::max(ymax,val(obs[o],"total"));
  frame->SetMaximum(ymax*1.25); frame->SetMinimum(0); frame->Draw();
  TLegend lg(0.62,0.60,0.88,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0);
  int ns=src.size();
  std::vector<TH1D*> keep;
  for(int s=0;s<ns;s++){
    TH1D* h=new TH1D(Form("h_%d",s),"",5,0,5); keep.push_back(h);
    for(int o=0;o<5;o++) h->SetBinContent(o+1,val(obs[o],src[s].key));
    int col=TColor::GetColor(src[s].hex);
    h->SetMarkerColor(col); h->SetLineColor(col); h->SetMarkerStyle(20+s); h->SetMarkerSize(1.3); h->SetLineWidth(2);
    h->Draw("P L same"); lg.AddEntry(h,src[s].lab,"pl");
  }
  lg.Draw();
  TString out=Form("../report/figures/fw_systbreak_%s.pdf",
    std::string(cfg)=="fhc5"?"fhc":std::string(cfg)=="rhcfull"?"rhc":"comb");
  c.SaveAs(out); printf("wrote %s\n",out.Data());
}
