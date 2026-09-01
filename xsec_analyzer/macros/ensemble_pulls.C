// ensemble_pulls.C -- pull and coverage from the fake-data ensemble.
//
// The closure test in the note uses ONE throw. It shows the extraction reproduces its
// input; it cannot show whether the quoted uncertainty is the right size, because a
// single residual has no width. This measures that.
//
// Reference: each throw's own h_genie_tune, i.e. the FIXED central-value truth smeared
// by THAT throw's A_C. Using h_fakedata_truth instead would be wrong -- it follows the
// Poisson draw, so numerator and reference fluctuate together and the pull width comes
// out too small. (Measured: tune varies 0.04% between throws, purely through A_C;
// fakedata_truth varies 0.6%.)
//
//   root -l -b -q 'macros/ensemble_pulls.C(32)'
void ensemble_pulls(int nthrow=32){
  std::vector<std::vector<double>> pulls;   // [bin][throw]
  int nb=0, used=0;
  for(int t=1;t<=nthrow;++t){
    TFile* f=TFile::Open(Form("/data/uboone/processed/ens/closure_hists_xsec_t%d.root",t));
    if(!f||f->IsZombie()) continue;
    TH1* u=(TH1*)f->Get("h_unfolded_nuwro");
    TH1* r=(TH1*)f->Get("h_genie_tune");
    if(!u||!r){ f->Close(); continue; }
    if(nb==0){ nb=u->GetNbinsX(); pulls.assign(nb,{}); }
    for(int b=1;b<=nb;++b){
      double e=u->GetBinError(b);
      if(e>0) pulls[b-1].push_back((u->GetBinContent(b)-r->GetBinContent(b))/e);
    }
    ++used; f->Close();
  }
  printf("\n  ensemble: %d members, %d bins\n\n", used, nb);
  if(used<3){ printf("  too few members\n"); return; }

  printf("  %4s %9s %9s %9s %9s\n","bin","mean","width","cov68","cov95");
  printf("  %s\n", std::string(46,'-').c_str());
  double gm=0,gw=0,g68=0,g95=0; int gn=0;
  for(int b=0;b<nb;++b){
    auto& p=pulls[b]; int n=p.size(); if(n<3) continue;
    double m=0; for(double x:p) m+=x; m/=n;
    double v=0; for(double x:p) v+=(x-m)*(x-m); double w=std::sqrt(v/(n-1));
    int c68=0,c95=0; for(double x:p){ if(std::fabs(x)<1.0) ++c68; if(std::fabs(x)<1.96) ++c95; }
    printf("  %4d %9.3f %9.3f %8.1f%% %8.1f%%\n", b, m, w, 100.*c68/n, 100.*c95/n);
    gm+=m; gw+=w; g68+=100.*c68/n; g95+=100.*c95/n; ++gn;
  }
  printf("  %s\n", std::string(46,'-').c_str());
  printf("  %4s %9.3f %9.3f %8.1f%% %8.1f%%\n","all",gm/gn,gw/gn,g68/gn,g95/gn);
  printf("\n  expected for correctly-sized uncertainties: mean 0, width 1, 68%%, 95%%\n");
  double w=gw/gn;
  if(w<0.8)      printf("  -> width %.2f < 1: uncertainties are CONSERVATIVE (over-covering)\n",w);
  else if(w>1.2) printf("  -> width %.2f > 1: uncertainties are UNDER-estimated\n",w);
  else           printf("  -> width %.2f: consistent with correctly-sized uncertainties\n",w);
  printf("  width is determined to about +-%.0f%% with %d members.\n", 100./std::sqrt(2.*(used-1)), used);
}
