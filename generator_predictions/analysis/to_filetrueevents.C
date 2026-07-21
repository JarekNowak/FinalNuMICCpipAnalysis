// to_filetrueevents.C — convert the combined per-observable differential NuWro
// prediction (physical x-axis, dsigma/dx in 10^-38 cm^2/Ar) into the format
// FileTrueEvents expects: one histogram per observable whose x-axis is the flat
// true-signal-bin index and whose content is the bin-INTEGRATED cross section
// (dsigma/dx * bin_width) in 10^-38 cm^2/Ar.
//
// FileTrueEvents does pred = content * conv_factor to get event counts, and
// UnfolderNuMI converts back with 1/(conv_factor*width); the round trip
// requires content = dsigma/dx * width. The number of bins must equal the
// number of SIGNAL true bins for that observable (the background true bin is
// not part of the prediction).
//
// Run:
//   root -l -b -q 'to_filetrueevents.C("out/nuwro_cc1pi_final.root","out/nuwro_pred_fte.root")'

void to_filetrueevents(const char* in_file, const char* out_file) {
  TFile fin(in_file);
  if (fin.IsZombie()) { printf("ERROR opening %s\n", in_file); return; }
  const char* names[5] = {"pmu","ppi","costhmu","costhpi","thmupi"};
  TFile fout(out_file, "recreate");
  for (auto n : names) {
    TH1D* h = (TH1D*)fin.Get(n);
    if (!h) { printf("missing %s\n", n); continue; }
    int nb = h->GetNbinsX();
    // Flat index histogram: nb bins from 0..nb.
    TH1D* hf = new TH1D(Form("%s_fte", n), n, nb, 0, nb);
    for (int b = 1; b <= nb; ++b) {
      double integ = h->GetBinContent(b) * h->GetBinWidth(b);
      double err   = h->GetBinError(b)   * h->GetBinWidth(b);
      hf->SetBinContent(b, integ);
      hf->SetBinError(b, err);
    }
    hf->Write();
    printf("  %-8s -> %s (%d bins, total sigma = %.4e [1e-38 cm^2/Ar])\n",
      n, hf->GetName(), nb, hf->Integral());
  }
  fout.Close();
  printf("wrote %s\n", out_file);
}
