// Combine numu+numubar per-nucleon predictions -> per-nucleus 10^-38 cm2/Ar,
// using the AUTHORITATIVE newg4 FHC flux integrals (>60 MeV):
//   Phi_numu = 4.43515e-10, Phi_numubar = 2.37644e-10 /POT/cm2
void combine_newg4(const char* numu_f,const char* numubar_f,const char* out_f){
  const double PN=4.43515e-10, PNB=2.37644e-10, PT=PN+PNB, A=40.0;
  TFile fn(numu_f), fnb(numubar_f), fo(out_f,"recreate");
  for(const char* n:{"pmu","ppi","costhmu","costhpi","thmupi"}){
    TH1D* a=(TH1D*)fn.Get(n); TH1D* b=(TH1D*)fnb.Get(n);
    if(!a||!b){printf("missing %s\n",n);continue;}
    TH1D* h=(TH1D*)a->Clone(n); h->Reset();
    for(int k=1;k<=a->GetNbinsX();++k){
      double v=A*(a->GetBinContent(k)*PN+b->GetBinContent(k)*PNB)/PT/1e-38;
      double e=A*sqrt(pow(a->GetBinError(k)*PN,2)+pow(b->GetBinError(k)*PNB,2))/PT/1e-38;
      h->SetBinContent(k,v); h->SetBinError(k,e);
    }
    h->Write(); printf("  %-8s integral(width)=%.4e [1e-38 cm2/Ar]\n",n,h->Integral("width"));
  }
  fo.Close();
}
