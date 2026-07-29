void combine_rhc(const char* numu_f,const char* numubar_f,const char* out_f){
  const double PN=2.92348e-10, PNB=3.52298e-10, PT=PN+PNB, A=40.0;
  TFile fn(numu_f), fnb(numubar_f), fo(out_f,"recreate");
  for(const char* n:{"pmu","ppi","costhmu","costhpi","thmupi"}){
    TH1D* a=(TH1D*)fn.Get(n); TH1D* b=(TH1D*)fnb.Get(n); if(!a||!b){printf("missing %s\n",n);continue;}
    TH1D* h=(TH1D*)a->Clone(n); h->Reset();
    for(int k=1;k<=a->GetNbinsX();++k){
      h->SetBinContent(k, A*(a->GetBinContent(k)*PN+b->GetBinContent(k)*PNB)/PT/1e-38);
      h->SetBinError(k,   A*sqrt(pow(a->GetBinError(k)*PN,2)+pow(b->GetBinError(k)*PNB,2))/PT/1e-38); }
    h->Write();
  }
  fo.Close(); printf("wrote %s\n",out_f);
}
