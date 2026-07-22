// flux-avg GENIE tot_cc per Ar over the cleaned newg4 flux, write to a text file
void genie_fluxavg(){
  const char* sp[2]={"../genie/splines_numu_ar.root","../genie/splines_numubar_ar.root"};
  const char* dir[2]={"nu_mu_Ar40","nu_mu_bar_Ar40"};
  const char* fx[2]={"flux_numu","flux_numubar"};
  TFile ff("flux_newg4_clean.root");
  FILE* o=fopen("genie_sigmatot.txt","w");
  for(int k=0;k<2;++k){
    TFile fs(sp[k]); TDirectory* d=(TDirectory*)fs.Get(dir[k]);
    if(!d){ printf("no dir %s in %s\n",dir[k],sp[k]); fs.ls(); continue; }
    TGraph* g=(TGraph*)d->Get("tot_cc");
    TH1D* h=(TH1D*)ff.Get(fx[k]);
    double num=0,den=0; for(int b=1;b<=h->GetNbinsX();++b){double E=h->GetBinCenter(b); if(E<0.06)continue; double phi=h->GetBinContent(b); num+=g->Eval(E)*phi; den+=phi; }
    double favg=num/den; printf("GENIE %s flux-avg tot_cc (per Ar) = %.4e [1e-38]\n",fx[k],favg);
    fprintf(o,"%s %.6e\n",fx[k],favg);
  }
  fclose(o);
}
