// Evaluate the pion-momentum regression on the HELD-OUT run (run_id==1, i.e. Run 2).
// Reports bias and RMS per true bin AND the migration diagonals, because a regression can
// flatter its own residuals while sculpting the response -- the trap documented for the
// tight golden-pion working point.
void pion_regression_eval() {
  TFile f("/data/uboone/processed/validate/pion_regr_sample.root");
  TTree* t=(TTree*)f.Get("reg"); if(!t){printf("  no sample\n");return;}
  float p_true,p_range,p_mcs,calo,len,wstd,wmean,wsep;
  float bragg_pion,bragg_mip,bragg_mu,bragg_p,tscore,ndau; int run_id,golden;
  t->SetBranchAddress("p_true",&p_true);   t->SetBranchAddress("p_range",&p_range);
  t->SetBranchAddress("p_mcs",&p_mcs);     t->SetBranchAddress("calo",&calo);
  t->SetBranchAddress("len",&len);         t->SetBranchAddress("wstd",&wstd);
  t->SetBranchAddress("wmean",&wmean);     t->SetBranchAddress("wsep",&wsep);
  t->SetBranchAddress("bragg_pion",&bragg_pion); t->SetBranchAddress("bragg_mip",&bragg_mip);
  t->SetBranchAddress("bragg_mu",&bragg_mu);     t->SetBranchAddress("bragg_p",&bragg_p);
  t->SetBranchAddress("tscore",&tscore);   t->SetBranchAddress("ndau",&ndau);
  t->SetBranchAddress("run_id",&run_id);   t->SetBranchAddress("golden",&golden);

  TMVA::Reader r("!Color:!Silent");
  r.AddVariable("p_range",&p_range); r.AddVariable("p_mcs",&p_mcs); r.AddVariable("calo",&calo);
  r.AddVariable("len",&len);         r.AddVariable("wstd",&wstd);   r.AddVariable("wmean",&wmean);
  r.AddVariable("wsep",&wsep);       r.AddVariable("bragg_pion",&bragg_pion);
  r.AddVariable("bragg_mip",&bragg_mip); r.AddVariable("bragg_mu",&bragg_mu);
  r.AddVariable("bragg_p",&bragg_p); r.AddVariable("tscore",&tscore); r.AddVariable("ndau",&ndau);
  r.BookMVA("BDTG","pion_regr/weights/PionMom_BDTG.weights.xml");

  const int NB=6; double e[NB+1]={0.175,0.250,0.320,0.420,0.550,0.800,1.500};
  std::vector<double> n(NB,0),sr(NB,0),sr2(NB,0),sp(NB,0),sp2(NB,0);
  std::vector<double> dr(NB,0),dp(NB,0);
  // two-bin scheme diagonals as well
  double n2[2]={0,0}, d2r[2]={0,0}, d2p[2]={0,0};
  auto bin=[&](double v){ for(int k=0;k<NB;++k) if(v>=e[k]&&v<e[k+1]) return k; return v<e[0]?0:-1; };

  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
    if(run_id!=1) continue;                       // held-out run
    double pred=r.EvaluateRegression("BDTG")[0];
    int tb=bin(p_true); if(tb<0||p_true>=e[NB]) continue;
    int rb=bin(p_range), pb=bin(pred);
    n[tb]++; sr[tb]+=p_range/p_true; sr2[tb]+=pow(p_range/p_true,2);
    sp[tb]+=pred/p_true; sp2[tb]+=pow(pred/p_true,2);
    if(rb==tb) dr[tb]++;
    if(pb==tb) dp[tb]++;
    int t2 = p_true<0.205?0:1, r2 = p_range<0.205?0:1, p2 = pred<0.205?0:1;
    n2[t2]++; if(r2==t2) d2r[t2]++; if(p2==t2) d2p[t2]++;
  }
  printf("\n  held-out run (Run 2), %.0f pion tracks\n", n[0]+n[1]+n[2]+n[3]+n[4]+n[5]);
  printf("\n  true p_pi        N    RANGE <r/t> RMS  diag    REGRESSION <p/t> RMS  diag\n");
  for(int b=0;b<NB;++b){ if(!n[b])continue;
    double mr=sr[b]/n[b], rr=sqrt(std::max(0.,sr2[b]/n[b]-mr*mr));
    double mp=sp[b]/n[b], rp=sqrt(std::max(0.,sp2[b]/n[b]-mp*mp));
    printf("  %5.3f-%5.3f %6.0f   %5.3f %5.3f %5.1f%%      %5.3f %5.3f %5.1f%%\n",
      e[b],e[b+1],n[b],mr,rr,100*dr[b]/n[b],mp,rp,100*dp[b]/n[b]); }
  printf("\n  two-bin scheme (split 0.205) diagonals:\n");
  printf("    range      %5.1f%% / %5.1f%%\n", 100*d2r[0]/n2[0], 100*d2r[1]/n2[1]);
  printf("    regression %5.1f%% / %5.1f%%\n", 100*d2p[0]/n2[0], 100*d2p[1]/n2[1]);
}
