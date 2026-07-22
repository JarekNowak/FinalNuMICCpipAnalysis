// Prepare newg4 flux for the generators: zero the sub-60-MeV bins (removes the
// 25-30 MeV artifact) and write cleaned histograms + NuWro beam_energy lines.
void prep_flux(){
  TFile fin("flux_newg4.root");
  const char* names[2]={"numu_cv_fhc","numubar_cv_fhc"};
  const char* tags[2]={"numu","numubar"};
  TFile fout("flux_newg4_clean.root","recreate");
  for(int k=0;k<2;++k){
    TH1D* h=(TH1D*)fin.Get(names[k]);
    int nb=h->GetNbinsX();
    // zero bins with upper edge <= 0.060 GeV (bins 1..12, 0-60 MeV)
    for(int b=1;b<=nb;++b) if(h->GetBinLowEdge(b)+h->GetBinWidth(b) <= 0.0600001) h->SetBinContent(b,0.0);
    fout.cd(); TH1D* hc=(TH1D*)h->Clone(Form("flux_%s",tags[k])); hc->Write();
    // NuWro beam_energy line: range in MeV, 2000 bin heights
    FILE* o=fopen(Form("nuwro_%s_beam.txt",tags[k]),"w");
    fprintf(o,"beam_energy = %.0f %.0f", h->GetXaxis()->GetXmin()*1000., h->GetXaxis()->GetXmax()*1000.);
    for(int b=1;b<=nb;++b) fprintf(o," %.6g", h->GetBinContent(b));
    fprintf(o,"\n"); fclose(o);
    printf("%s: cleaned integral(>60MeV)=%.4e  peak bin now=%d (%.3f GeV)\n",
      tags[k], h->Integral(), h->GetMaximumBin(), h->GetBinCenter(h->GetMaximumBin()));
  }
  fout.Close();
}
