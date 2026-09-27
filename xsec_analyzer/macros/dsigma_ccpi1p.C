// dsigma_ccpi1p.C — differential cross-section montage for the proton-tagged
// CC1mu1pi1p (W/TKI) observables. Reads the UnfolderNuMI closure dumps
// closure_hists_xsec_ccpi1p_<CFG>_<obs>.root (h_unfolded_nuwro = unfolded fake data,
// h_fakedata_truth = A_C-smeared truth, h_genie_tune = MicroBooNE tune) and draws one
// multi-panel figure per config. The A_C-smeared generator predictions are drawn wherever the
// sidecar carries them (h_gen_*). Auto-scales the y-axis to include every curve (no clipping).
//   usage: root -l -b -q 'macros/dsigma_ccpi1p.C("FHC5")'  (or RHCFULL / COMB)
// 2026-09-27 (review 6.1-6.3): physics notation and units from dsigma_style.h, larger text, one
// legend per figure, open-ended top bins hatched and labelled; the "wtki_noW" set is also written
// as two three-panel figures, _a (W_had, p_p, p_n) and _b (delta p_T, delta alpha_T, delta phi_T).
// Review terminology: "CV pseudo-data"; a panel without generator curves carries a note.
#include <vector>
#include "dsigma_style.h"
// set = "wtki" (the six W/TKI observables, default) or "incl" (the five inclusive-kinematic
// observables measured on the proton-tagged selection, which were re-run 2026-08-19)
void dsigma_ccpi1p(const char* cfg = "FHC5", const char* set = "wtki") {
  const bool incl = ( std::string(set) == "incl" );
  // "wtki_noW": the five proton-tagged observables that are results (W_had + the four TKI);
  // the differential W_pipr shape is withdrawn and its panel belongs in the supplement only.
  const bool noW  = ( std::string(set) == "wtki_noW" );
  const char* PROC = "/data/uboone/processed/";
  // The four TKI observables are RELEASED in the coarsened two-bin scheme; their closure
  // files are keyed dpt2bin/dalphat2bin/dphit2bin/pn2bin. The un-suffixed names are the
  // withdrawn fine-binned extraction of 2026-08-18, which predates the beam-axis fix
  // (TKI about detector z) -- reading them drew a "data" curve that was a different
  // observable from the generator curves. Those sidecars are quarantined.
  // 2026-09-26: delta phi_T in three bins (dphit3bin) and W_had in two, under the 0.50 criterion
  const char* obs_wtki[6] = {"Wpipr","Whad","dpt2bin","dalphat2bin","dphit3bin","pn2bin"};
  const char* obs_noW[6]  = {"Whad","dpt2bin","dalphat2bin","dphit3bin","pn2bin","pp"};   // pp: proton momentum (2026-09-26)
  const char* obs_incl[6] = {"pmu","ppi2bin","costhmu","costhpi","thmupi","thmupi"};
  const char* obs[6];
  for ( int i = 0; i < 6; ++i ) obs[i] = incl ? obs_incl[i] : ( noW ? obs_noW[i] : obs_wtki[i] );
  const int n_panel = incl ? 5 : 6;
  // panel titles and axis labels (physics notation, units) come from ds_obs() in dsigma_style.h
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  std::vector<TObject*> keep;

  std::vector<DsPanel> P(n_panel);
  for (int o = 0; o < n_panel; ++o) {
    P[o].key = obs[o];
    TFile* f = TFile::Open(Form("%sclosure_hists_xsec_ccpi1p_%s_%s.root", PROC, cfg, obs[o]));
    if (!f || f->IsZombie()) { printf("  missing closure %s %s\n", cfg, obs[o]); continue; }
    TH1D* hunf = (TH1D*)f->Get("h_unfolded_nuwro");
    TH1D* htru = (TH1D*)f->Get("h_fakedata_truth");
    TH1D* htun = (TH1D*)f->Get("h_genie_tune");
    if (!hunf) { f->Close(); continue; }
    hunf = (TH1D*)hunf->Clone(); hunf->SetDirectory(0); keep.push_back(hunf);
    if (htru) { htru=(TH1D*)htru->Clone(); htru->SetDirectory(0); keep.push_back(htru); }
    if (htun) { htun=(TH1D*)htun->Clone(); htun->SetDirectory(0); keep.push_back(htun); }
    hunf->SetMarkerStyle(20); hunf->SetMarkerSize(1.1);
    hunf->SetLineColor(kBlack); hunf->SetMarkerColor(kBlack);
    // Generator overlays: use the A_C-SMEARED predictions dumped by the unfolder
    // (h_gen_<Label>), NOT the raw truth-level FTE files. The data and the uB tune both
    // live in the A_C-smeared measurement space, so every model must be smeared the same
    // way; drawing raw generators here made them look systematically high.
    const char* gens[4]={"h_gen_GENIE","h_gen_GiBUU","h_gen_NEUT","h_gen_NuWro"};
    int gcol[4]={TColor::GetColor("#0072B2"),TColor::GetColor("#009E73"),TColor::GetColor("#CC79A7"),TColor::GetColor("#D55E00")};
    int gsty[4]={1,2,7,9};
    std::vector<TH1D*> gh; std::vector<int> gi;
    for (int g=0;g<4;++g){
      TH1D* hg0=(TH1D*)f->Get(gens[g]);
      if(hg0 && hg0->GetNbinsX()==hunf->GetNbinsX()){
        TH1D* hg=(TH1D*)hg0->Clone(Form("g_%s_%d_%d",cfg,o,g)); hg->SetDirectory(0);
        hg->SetLineColor(gcol[g]); hg->SetLineStyle(gsty[g]); hg->SetLineWidth(2); hg->SetMarkerSize(0);
        gh.push_back(hg); gi.push_back(g); keep.push_back(hg);
      }
    }
    f->Close();
    // A panel without generator curves says why: the RHC and combined W_had/TKI sidecars carry none
    // (the FHC ones do), and the inclusive-kinematic observables of this sample have none at all.
    if ( gh.empty() ) {
      bool fhc = false;
      if ( std::string(cfg) != "FHC5" ) {
        TFile* ff = TFile::Open(Form("%sclosure_hists_xsec_ccpi1p_FHC5_%s.root", PROC, obs[o]));
        if ( ff && !ff->IsZombie() ) { for ( auto gname : gens ) if ( ff->Get(gname) ) fhc = true; ff->Close(); }
      }
      P[o].gen_note = fhc ? "generator predictions: FHC only" : "no generator predictions";
    }
    // y-axis max over data(+error), truth, tune and every generator curve
    double ymax = 0.;
    for (int b=1;b<=hunf->GetNbinsX();++b) ymax = std::max(ymax, hunf->GetBinContent(b)+hunf->GetBinError(b));
    if (htru) for (int b=1;b<=htru->GetNbinsX();++b) ymax = std::max(ymax, htru->GetBinContent(b));
    if (htun) for (int b=1;b<=htun->GetNbinsX();++b) ymax = std::max(ymax, htun->GetBinContent(b));
    for (auto hg:gh) for (int b=1;b<=hg->GetNbinsX();++b) ymax = std::max(ymax, hg->GetBinContent(b));
    hunf->SetMinimum(0); hunf->SetMaximum(1.35*ymax);
    if (htru) { htru->SetLineColor(kGray+2); htru->SetLineWidth(2); htru->SetLineStyle(2); }
    if (htun) { htun->SetLineColor(kBlack); htun->SetLineWidth(2); }
    P[o].data = hunf; P[o].truth = htru; P[o].tune = htun; P[o].gen = gh; P[o].gidx = gi;
    // open-ended top bin: its physical lower edge, from the result file
    if ( ds_obs(obs[o])->open_top ) P[o].open_lo = hunf->GetXaxis()->GetBinLowEdge(hunf->GetNbinsX());
  }
  std::vector<DsPanel*> all;
  for (auto& p : P) if (p.data) all.push_back(&p);

  // text in pixels for a 1800 px wide figure printed at the text width (~9 pt titles); 700 px tall
  // panels leave room for the full vertical-axis title with its units
  const DsStyle S = { /*title*/34, /*label*/31, /*binlabel*/31, /*legend*/34, /*ptitle*/38, /*note*/30,
                      /*lm*/0.215, /*rm*/0.045, /*tm*/0.095, /*bm*/0.16, /*xoff*/1.3, /*yoff*/1.85 };
  TString out = Form("unfold_output/dsigma_ccpi1p%s_%s.pdf", incl ? "_incl" : ( noW ? "_noW" : "" ), cfg);
  ds_figure(out.Data(), all, 3, 2, 1800, 700, 150, S, cfg, keep);
  if ( noW ) {
    // the same six panels as two three-panel figures, large enough to read at the text width
    auto pick = [&](std::vector<const char*> keys) {
      std::vector<DsPanel*> v;
      for (auto k : keys) for (auto& p : P) if (p.data && p.key == k) v.push_back(&p);
      return v;
    };
    TString a = out, b = out; a.ReplaceAll(".pdf", "_a.pdf"); b.ReplaceAll(".pdf", "_b.pdf");
    ds_figure(a.Data(), pick({"Whad","pp","pn2bin"}),               3, 1, 1800, 700, 150, S, cfg, keep);
    ds_figure(b.Data(), pick({"dpt2bin","dalphat2bin","dphit3bin"}), 3, 1, 1800, 700, 150, S, cfg, keep);
  }
}
