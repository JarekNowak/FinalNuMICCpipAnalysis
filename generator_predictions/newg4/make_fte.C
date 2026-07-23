// Convert a *_newg4_final.root (per-observable dsigma/dx histos) to the
// FileTrueEvents (FTE) format the xsec configs consume: a flat-index histogram
// <obs>_fte (nbins = observable bins, x-range 0..N) whose bin content is the
// per-bin cross section dsigma/dx * binwidth [1e-38 cm^2/Ar]. Matches how the
// NuWro/GENIE FTE files were built (verified against nuwro_newg4_fte.root).
void make_fte(const char* infile, const char* outfile){
  const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  TFile fi(infile), fo(outfile,"recreate");
  for(int o=0;o<5;o++){
    TH1D* h=(TH1D*)fi.Get(obs[o]);
    if(!h){printf("missing %s in %s\n",obs[o],infile);continue;}
    int nb=h->GetNbinsX();
    TH1D* fte=new TH1D(Form("%s_fte",obs[o]),"",nb,0,nb);
    for(int b=1;b<=nb;b++){
      fte->SetBinContent(b, h->GetBinContent(b)*h->GetBinWidth(b));
      fte->SetBinError(b,   h->GetBinError(b)*h->GetBinWidth(b));
    }
    fte->SetDirectory(&fo); fte->Write();
  }
  fo.Close();
  printf("wrote %s\n",outfile);
}
