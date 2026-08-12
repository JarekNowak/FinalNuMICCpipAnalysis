// combine_wtki.C — combine per-nucleon numu+numubar W/TKI predictions into the FHC
// per-nucleus 10^-38 cm^2/Ar prediction and write it directly as flat-index FTE
// histograms (content = dsigma/dx * physical binwidth), the format dsigma_ccpi1p.C
// overlays. FHC flux integrals (>60 MeV): Phi_numu=4.43515e-10, Phi_numubar=2.37644e-10.
//   root -l -b -q 'combine_wtki.C("newg4/genie_1p_numu.root","newg4/genie_1p_numubar.root","newg4/genie_wtki_fte.root")'
void combine_wtki(const char* numu_f,const char* numubar_f,const char* out_f){
  const double PN=4.43515e-10, PNB=2.37644e-10, PT=PN+PNB, A=40.0;
  const char* OBS[6]={"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  TFile fn(numu_f), fnb(numubar_f), fo(out_f,"recreate");
  for(auto n:OBS){
    TH1D* a=(TH1D*)fn.Get(n); TH1D* b=(TH1D*)fnb.Get(n);
    if(!a||!b){printf("missing %s\n",n);continue;}
    int nb=a->GetNbinsX();
    TH1D* fte=new TH1D(Form("%s_fte",n),"",nb,0,nb);
    for(int k=1;k<=nb;++k){
      double dsdx=A*(a->GetBinContent(k)*PN+b->GetBinContent(k)*PNB)/PT/1e-38;      // per-Ar 1e-38
      double e  =A*sqrt(pow(a->GetBinError(k)*PN,2)+pow(b->GetBinError(k)*PNB,2))/PT/1e-38;
      double w=a->GetBinWidth(k);                                                   // physical binwidth
      fte->SetBinContent(k,dsdx*w); fte->SetBinError(k,e*w);
    }
    fte->SetDirectory(&fo); fte->Write();
    printf("  %-8s sigma_int=%.4f [1e-38 cm2/Ar]\n",n,fte->Integral());
  }
  fo.Close(); printf("wrote %s\n",out_f);
}
