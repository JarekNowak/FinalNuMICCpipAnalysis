// cutflow_perrun.C — per-run cut-flow figures for all 5 reco observables.
// Reads the per-cut RECO cut-flow histograms (h_cf_<obs>_<stage>) that the extended
// CC1mu1piXp selection now writes into each processed run file. For each observable
// it makes ONE figure with 10 panels (one per cut stage); each panel stacks the
// per-run contributions. Run AFTER the framework is rebuilt and the run files are
// reprocessed.
//   usage: root -l -b -q 'macros/cutflow_perrun.C("fhc")'   // or "rhc"
void cutflow_perrun(const char* mode = "fhc") {
  const char* PROC = "/data/uboone/processed/";
  // (label, file) per run for each mode
  std::vector<std::pair<TString,TString>> runs;
  if (TString(mode) == "fhc") {
    runs = {
      {"Run1", "xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root"},
      {"Run2", "xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root"},
      {"Run4", "xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"},
      {"Run5", "xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root"} };
  } else {
    runs = {
      {"Run1",  "xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root"},
      {"Run2",  "xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root"},
      {"Run3",  "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root"}, // aa-ae summed below
      {"Run4a", "xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root"},
      {"Run4b", "xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root"},
      {"Run4c", "xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root"} };
  }
  const char* obs[5]   = {"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* obsTtl[5]= {"p_{#mu}","p_{#pi}","cos#theta_{#mu}","cos#theta_{#pi}","#theta_{#mu#pi}"};
  const char* stage[10]= {"cut0_none","cut1_vertex","cut2_topology","cut3_tracklike",
                          "cut4_pioncontained","cut5_muongap","cut6_piongap",
                          "cut7_shower","cut8_opening","cut9_final"};
  int col[8] = {kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kGray+2};
  gStyle->SetOptStat(0);

  for (int o = 0; o < 5; ++o) {
    TCanvas c(Form("cf_%s_%s",mode,obs[o]), "", 1600, 900);
    c.Divide(5, 2);
    std::vector<TH1D*> keep;
    for (int s = 0; s < 10; ++s) {
      c.cd(s+1);
      THStack* hs = new THStack(Form("hs_%s_%d",obs[o],s), Form("%s;%s;Events",stage[s],obsTtl[o]));
      TLegend* lg = new TLegend(0.55,0.60,0.88,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
      for (size_t r = 0; r < runs.size(); ++r) {
        TFile* f = TFile::Open(PROC + runs[r].second);
        if (!f || f->IsZombie()) continue;
        TH1D* h = (TH1D*)f->Get(Form("h_cf_%s_%s", obs[o], stage[s]));
        if (!h) { f->Close(); continue; }
        TH1D* hc = (TH1D*)h->Clone(Form("c_%s_%d_%zu",obs[o],s,r)); hc->SetDirectory(0);
        // Run3 rhc is split aa-ae: add ab-ae onto the aa clone
        if (runs[r].first == "Run3") {
          for (const char* p : {"ab","ac","ad","ae"}) {
            TFile* f2 = TFile::Open(Form("%sxsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",PROC,p));
            if (f2 && !f2->IsZombie()) { TH1D* h2=(TH1D*)f2->Get(Form("h_cf_%s_%s",obs[o],stage[s])); if(h2) hc->Add(h2); f2->Close(); }
          }
        }
        hc->SetFillColor(col[r % 8]); hc->SetLineColor(kBlack); hc->SetLineWidth(1);
        hs->Add(hc); keep.push_back(hc);
        lg->AddEntry(hc, Form("%s (%.0f)", runs[r].first.Data(), hc->Integral()), "f");
        f->Close();
      }
      hs->Draw("hist");
      lg->Draw();
    }
    TString out = Form("unfold_output/cutflow_%s_%s.pdf", mode, obs[o]);
    c.SaveAs(out);
    printf("wrote %s\n", out.Data());
  }
}
