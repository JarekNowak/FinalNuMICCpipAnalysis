// export_flux.C — convert the framework NuMI flux histograms into the flux
// formats the generators need.
//
//   NuWro : a `beam_energy = Elow_MeV Ehigh_MeV c1 c2 ... cN` line, contents
//           being the per-bin flux heights (arbitrary norm; NuWro samples the
//           shape and reports an absolute cross section).
//
// We combine numu + numubar, because the CC1mu1piXp signal is |pdg|==14, i.e.
// both contribute. The generators are run for numu (pdg 14); the numubar
// contribution to CC1pi+ is small in FHC and is folded in later if needed. For
// the proving run we use the numu shape.
//
// Run:  root -l -b -q export_flux.C

void export_flux() {
  TFile f("uboone_numi_flux_histograms.root");
  TH1D* h_numu    = (TH1D*)f.Get("h_numu");
  TH1D* h_numubar = (TH1D*)f.Get("h_numubar");
  if (!h_numu) { printf("ERROR: h_numu not found\n"); return; }

  int nb = h_numu->GetNbinsX();
  double elo = h_numu->GetXaxis()->GetXmin(); // GeV
  double ehi = h_numu->GetXaxis()->GetXmax();
  printf("h_numu: %d bins, [%.1f, %.1f] GeV, integral(width) numu=%.4e",
    nb, elo, ehi, h_numu->Integral("width"));
  if (h_numubar) printf(", numubar=%.4e", h_numubar->Integral("width"));
  printf("\n");

  // NuWro beam_energy line: range in MeV, then nb bin contents (numu only).
  FILE* out = fopen("nuwro_numu_beam.txt", "w");
  fprintf(out, "beam_energy = %.0f %.0f", elo * 1000., ehi * 1000.);
  for (int i = 1; i <= nb; ++i) fprintf(out, " %.6g", h_numu->GetBinContent(i));
  fprintf(out, "\n");
  fclose(out);
  printf("wrote nuwro_numu_beam.txt (%d bins over %.0f-%.0f MeV)\n",
    nb, elo * 1000., ehi * 1000.);

  // Also a plain two-column (E_GeV  flux) table for GENIE/GiBUU/NEUT use later.
  FILE* t = fopen("numi_numu_flux.txt", "w");
  for (int i = 1; i <= nb; ++i)
    fprintf(t, "%.5f %.6e\n", h_numu->GetBinCenter(i), h_numu->GetBinContent(i));
  fclose(t);
  printf("wrote numi_numu_flux.txt (E_GeV flux, %d rows)\n", nb);
}
