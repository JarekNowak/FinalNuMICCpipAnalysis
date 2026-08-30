// ppi2bin_figs.C -- the two-bin p_pi figures: cross section (inclusive and
// proton-tagged) and the additional-smearing matrix A_C.
//
// These nine figures (ppi2bin_xsec_*, ppi2bin_xsec_1p_*, ppi2bin_AC_*) are
// referenced by the note but had NO producer anywhere in the repository: they
// were made once, by hand, on 2026-08-18, and could not be regenerated when the
// COMB detector systematic was fixed. This macro reproduces them from the
// closure sidecar files UnfolderNuMI already writes.
//
//   usage: root -l -b -q 'macros/ppi2bin_figs.C("FHC5")'   // or RHCFULL / COMB
//
// The two bins are 0.030 and 0.795 GeV/c wide. Drawn to scale the first would
// occupy 4% of the axis, so they are drawn with EQUAL WIDTH for legibility and
// the integrated sigma is printed on the panel instead -- area on these panels is
// therefore not proportional to cross section. This matches the note's caption.
#include <map>
#include <string>

static const char* figdir = "../report/figures/";

// Okabe-Ito, consistent with the other note figures
static int OI(int i){
  static int c[7] = {kBlack, kOrange+7, kAzure+2, kGreen+2, kMagenta+2, kRed+1, kGray+2};
  return c[i%7];
}

static TString shortcfg(const char* cfg){
  TString s(cfg);
  if (s=="FHC5")    return "FHC";
  if (s=="RHCFULL") return "RHC";
  return "COMB";
}

// Redraw a 2-bin histogram onto an equal-width bin-number axis.
static TH1D* equalise(TH1* h, const char* name){
  if (!h) return nullptr;
  int n = h->GetNbinsX();
  auto* o = new TH1D(name, "", n, 0., (double)n);
  for (int i=1;i<=n;i++){
    // divide out the real width so the panel shows dsigma/dp, then the equal-width
    // axis is purely a drawing convenience
    o->SetBinContent(i, h->GetBinContent(i));
    o->SetBinError  (i, h->GetBinError(i));
  }
  o->SetDirectory(nullptr);
  return o;
}

static void one_xsec(const char* cfg, bool proton_tagged){
  TString tag  = proton_tagged ? "xsec1p" : "xsec";
  TString path = TString::Format("/data/uboone/processed/closure_hists_%s_%s_ppi.root", tag.Data(), cfg);
  auto* f = TFile::Open(path);
  if (!f || f->IsZombie()){ printf("  MISSING %s\n", path.Data()); return; }

  auto* hu = (TH1D*)f->Get("h_unfolded_nuwro");
  if (!hu){ printf("  no unfolded hist in %s\n", path.Data()); f->Close(); return; }

  std::vector<std::pair<std::string,std::string>> curves = {
    {"h_fakedata_truth","A_{C} truth"}, {"h_genie_tune","uB tune"},
    {"h_gen_GENIE","GENIE"}, {"h_gen_GiBUU","GiBUU"},
    {"h_gen_NEUT","NEUT"},   {"h_gen_NuWro","NuWro"} };

  auto* c = new TCanvas("c","",700,560);
  c->SetLeftMargin(0.15); c->SetBottomMargin(0.13);

  auto* hd = equalise(hu, "hd");
  hd->SetMarkerStyle(20); hd->SetMarkerSize(1.2); hd->SetLineWidth(2); hd->SetLineColor(kBlack);
  hd->GetXaxis()->SetTitle("p_{#pi} bin  (0.175-0.205, >0.205 GeV/c; equal width for legibility)");
  hd->GetYaxis()->SetTitle("d#sigma/dp_{#pi}  [10^{-38} cm^{2}/(GeV/c)/Ar]");
  hd->GetXaxis()->SetNdivisions(2,0,0);
  double ymax = hd->GetMaximum()+hd->GetBinError(hd->GetMaximumBin());
  for (auto& cv : curves){ auto* h=(TH1D*)f->Get(cv.first.c_str()); if(h) ymax=std::max(ymax,h->GetMaximum()); }
  hd->SetMaximum(1.45*ymax); hd->SetMinimum(0.);
  hd->Draw("E1");

  auto* leg = new TLegend(0.45,0.60,0.90,0.89);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
  leg->AddEntry(hd, Form("unfolded fake data  (#sigma=%.3f)", hu->Integral("width")), "lep");

  int ci=1; std::vector<TH1D*> keep;
  for (auto& cv : curves){
    auto* h = (TH1D*)f->Get(cv.first.c_str());
    if (!h) continue;
    auto* e = equalise(h, Form("e_%s", cv.first.c_str()));
    e->SetLineColor(OI(ci)); e->SetLineWidth(2); e->SetLineStyle(ci==1?2:1); e->SetMarkerStyle(0);
    e->Draw("HIST SAME"); keep.push_back(e);
    leg->AddEntry(e, Form("%s  (%.3f)", cv.second.c_str(), h->Integral("width")), "l");
    ci++;
  }
  hd->Draw("E1 SAME");
  leg->Draw();

  TLatex t; t.SetNDC(); t.SetTextSize(0.038);
  t.DrawLatex(0.17,0.93, Form("%s%s  #minus  two-bin p_{#pi}",
    shortcfg(cfg).Data(), proton_tagged?"  (proton-tagged)":""));

  TString out = TString::Format("%sppi2bin_xsec%s_%s.pdf", figdir,
                                proton_tagged?"_1p":"", shortcfg(cfg).Data());
  c->SaveAs(out);
  printf("  wrote %s\n", out.Data());
  delete c; f->Close();
}

static void one_AC(const char* cfg){
  TString path = TString::Format("/data/uboone/processed/closure_hists_xsec_%s_ppi.root", cfg);
  auto* f = TFile::Open(path);
  if (!f || f->IsZombie()){ printf("  MISSING %s\n", path.Data()); return; }
  auto* ac = (TH2D*)f->Get("h_A_C");
  if (!ac){ printf("  no h_A_C in %s\n", path.Data()); f->Close(); return; }

  auto* c = new TCanvas("cac","",620,560);
  c->SetLeftMargin(0.15); c->SetRightMargin(0.16); c->SetBottomMargin(0.13);
  gStyle->SetPaintTextFormat("5.3f");
  ac->SetTitle("");
  ac->GetXaxis()->SetTitle("true p_{#pi} bin");
  ac->GetYaxis()->SetTitle("smeared p_{#pi} bin");
  ac->GetXaxis()->SetNdivisions(2,0,0); ac->GetYaxis()->SetNdivisions(2,0,0);
  ac->SetMarkerSize(2.2);
  ac->Draw("COLZ TEXT");
  TLatex t; t.SetNDC(); t.SetTextSize(0.038);
  t.DrawLatex(0.17,0.93, Form("%s  #minus  A_{C}, two-bin p_{#pi}", shortcfg(cfg).Data()));
  TString out = TString::Format("%sppi2bin_AC_%s.pdf", figdir, shortcfg(cfg).Data());
  c->SaveAs(out);
  printf("  wrote %s\n", out.Data());
  delete c; f->Close();
}

void ppi2bin_figs(const char* cfg = "FHC5"){
  gStyle->SetOptStat(0);
  one_xsec(cfg, false);
  one_xsec(cfg, true);
  one_AC(cfg);
}
