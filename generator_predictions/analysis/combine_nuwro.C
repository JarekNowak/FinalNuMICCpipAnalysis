// combine_nuwro.C — combine the numu and numubar per-nucleon NuWro predictions
// into the framework's per-nucleus, combined-flux, 10^-38 cm^2/Ar units.
//
// The framework measures  dsigma/dx = N_signal / (T * Phi_total * dx),  where
// T is the number of argon NUCLEI and Phi_total is the integrated numu+numubar
// flux. A generator prediction that matches must be the flux-combined,
// per-nucleus differential cross section:
//
//   dsigma/dx[bin] = A * ( ds_numu[bin]*Phi_numu + ds_numubar[bin]*Phi_numubar )
//                        / Phi_total / 1e-38
//
// where ds_X[bin] is the per-nucleon dsigma/dx from nuwro_cc1pi (cm^2/nucleon),
// A = 40 nucleons per Ar, and the Phi_X are the per-POT integrated fluxes.
//
// Run:
//   root -l -b -q 'combine_nuwro.C("out/nuwro_numu_pred.root","out/nuwro_numubar_pred.root","out/nuwro_cc1pi_final.root")'

void combine_nuwro(const char* numu_file, const char* numubar_file,
                   const char* out_file)
{
  // Per-POT integrated fluxes (Integral("width") of h_numu / h_numubar from
  // uboone_numi_flux_histograms.root).
  const double PHI_NUMU    = 1.476120e-9;
  const double PHI_NUMUBAR = 8.963927e-10;
  const double PHI_TOTAL   = PHI_NUMU + PHI_NUMUBAR;
  const double A_AR        = 40.0;   // nucleons per argon nucleus

  TFile fnu(numu_file);
  TFile fnb(numubar_file);
  if (fnu.IsZombie() || fnb.IsZombie()) { printf("ERROR opening inputs\n"); return; }

  const char* names[5] = {"pmu","ppi","costhmu","costhpi","thmupi"};
  TFile fout(out_file, "recreate");

  for (auto n : names) {
    TH1D* h_nu = (TH1D*)fnu.Get(n);
    TH1D* h_nb = (TH1D*)fnb.Get(n);
    if (!h_nu || !h_nb) { printf("missing %s\n", n); continue; }

    TH1D* h = (TH1D*)h_nu->Clone(n);
    h->Reset();
    for (int b = 1; b <= h_nu->GetNbinsX(); ++b) {
      double dnu = h_nu->GetBinContent(b);
      double dnb = h_nb->GetBinContent(b);
      double enu = h_nu->GetBinError(b);
      double enb = h_nb->GetBinError(b);
      double val = A_AR * (dnu*PHI_NUMU + dnb*PHI_NUMUBAR) / PHI_TOTAL / 1e-38;
      double err = A_AR * std::sqrt(std::pow(enu*PHI_NUMU,2)
                       + std::pow(enb*PHI_NUMUBAR,2)) / PHI_TOTAL / 1e-38;
      h->SetBinContent(b, val);
      h->SetBinError(b, err);
    }
    h->Write();
    printf("  %-8s integral(width) = %.4e  [1e-38 cm^2/Ar]\n", n, h->Integral("width"));
  }
  fout.Close();
  printf("wrote %s\n", out_file);
}
