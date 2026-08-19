// w_binning_scan.C -- how many W bins the resolution supports. Same method and criterion
// as macros/ppi_binning_scan.C. Result: W_pipr supports TWO bins (83.2/77.9%), against the
// six currently used (87.6/22.3/11.7/12.0/10.0/2.3%). See the note, Sec. w_binning.
//   usage: root -l -b -q "macros/w_binning_scan.C(\"W_pipr\",1.08,2.60)"
// How many W_pipr bins can the resolution actually support?
// Same method as macros/ppi_binning_scan.C: fine 2D (true,reco) histogram, cumulative sums,
// exhaustive edge scan, criterion = migration diagonal > 0.68 with >=100 events per bin.
#include <functional>
#include <vector>
#include <algorithm>
void w_scan(const char* var="W_pipr", double LO=1.08, double HI=2.60, double STEP=0.01,
            double MINW=0.05, int MINEVT=100) {
  const double CRIT=0.68;
  TChain ch("stv_tree");
  for (auto r : {"Run1","Run2","Run4"})
    ch.Add(Form("/data/uboone/processed/w/xsec-ana-%s_fhc_new_numi_flux_fhc_pandora_ntuple.root", r));
  const char* S="CC1mu1pi1p";
  int NT=(int)std::lround((HI-LO)/STEP), NR=(int)std::lround((HI-LO)/STEP);
  ch.Draw(Form("%s_%s_reco : %s_%s_true >> h2(%d,%g,%g,%d,%g,%g)",S,var,S,var,NT,LO,HI,NR,LO,HI),
          Form("%s_Selected && %s_MC_Signal",S,S), "goff");
  TH2D* h2=(TH2D*)gDirectory->Get("h2");
  if(!h2||!h2->GetEntries()){printf("  no entries\n");return;}
  std::vector<std::vector<double>> C(NT+1, std::vector<double>(NR+1,0.));
  for(int i=1;i<=NT;++i) for(int j=1;j<=NR;++j)
    C[i][j]=h2->GetBinContent(i,j)+C[i-1][j]+C[i][j-1]-C[i-1][j-1];
  auto box=[&](int t0,int t1,int r0,int r1){return C[t1][r1]-C[t0][r1]-C[t1][r0]+C[t0][r0];};
  auto ix=[&](double e){return (int)std::lround((e-LO)/STEP);};
  printf("\n  %s: %.0f selected signal events, range %.2f-%.2f\n", var, box(0,NT,0,NR), LO, HI);
  auto eval=[&](const std::vector<double>& e, std::vector<double>& d, std::vector<double>& n){
    int nb=e.size()-1; d.assign(nb,0.); n.assign(nb,0.);
    for(int b=0;b<nb;++b){
      int t0=ix(e[b]), t1=(b==nb-1?NT:ix(e[b+1]));
      int r0=(b==0?0:ix(e[b])), r1=(b==nb-1?NR:ix(e[b+1]));
      double tot=box(t0,t1,0,NR); n[b]=tot;
      d[b]= tot>0? box(t0,t1,r0,r1)/tot : 0.;
    }
  };
  std::vector<double> cand;
  for(double e=LO+MINW;e<=HI-MINW+1e-9;e+=STEP) cand.push_back(e);
  for(int nb=2;nb<=4;++nb){
    std::vector<double> best,bd,bn,d,n; double bm=-1;
    std::vector<int> idx(nb-1,0);
    std::function<void(int,int)> rec=[&](int pos,int start){
      if(pos==nb-1){
        std::vector<double> e{LO};
        for(int k=0;k<nb-1;++k) e.push_back(cand[idx[k]]);
        e.push_back(HI);
        for(size_t k=1;k<e.size();++k) if(e[k]-e[k-1]<MINW-1e-9) return;
        eval(e,d,n);
        if(*std::min_element(n.begin(),n.end())<MINEVT) return;
        double m=*std::min_element(d.begin(),d.end());
        if(m>bm){bm=m;best=e;bd=d;bn=n;}
        return;
      }
      for(int i=start;i<(int)cand.size();++i){idx[pos]=i;rec(pos+1,i+1);}
    };
    rec(0,0);
    printf("\n  --- %d bins ---\n", nb);
    if(bm<0){printf("    none satisfies the >=%d events/bin floor\n",MINEVT);continue;}
    printf("    edges:"); for(size_t k=0;k+1<best.size();++k) printf(" %.2f",best[k]); printf(" -> open\n");
    printf("    diagonals:"); for(double x:bd) printf(" %5.1f%%",100*x);
    printf("   worst %5.1f%% %s\n",100*bm, bm>CRIT?"PASSES":"FAILS");
    printf("    events:   "); for(double x:bn) printf(" %6.0f",x); printf("\n");
  }
  // the current scheme, for comparison
  std::vector<double> cur{1.08,1.19,1.23,1.27,1.34,1.47,HI}, d,n;
  eval(cur,d,n);
  printf("\n  current 6-bin scheme diagonals:");
  for(double x:d) printf(" %5.1f%%",100*x);
  printf("\n");
}
