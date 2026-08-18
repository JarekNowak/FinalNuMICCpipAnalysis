// golden_diagonals_check.C -- migration diagonals with and without a PERFECT truth-level
// golden selection, computed on the LABELLED TRACK SAMPLE (golden_train_nolen.root), which
// is a looser sample than the analysis one: any truth-matched pion track passing purity>0.5,
// trk_score>0.5, len>5 cm, from events of any topology. Its diagonals are therefore lower
// than the analysis response and the two must not be compared. Written to settle where the
// note's 99/24/13/6/6 and 97/62 figures came from: they are this sample, not the analysis.
// Migration diagonals with and without a PERFECT truth-level golden-pion selection,
// computed from the labelled training sample (p_true, p_corr, golden).
void golden_check() {
  TFile f("/data/uboone/processed/validate/golden_train_nolen.root");
  TTree* t=(TTree*)f.Get("train");
  if(!t){printf("  no tree\n");return;}
  float p_true,p_corr; int golden;
  t->SetBranchAddress("p_true",&p_true); t->SetBranchAddress("p_corr",&p_corr);
  t->SetBranchAddress("golden",&golden);

  std::vector<std::vector<double>> schemes = {
    {0.175,0.250,0.320,0.420,0.550,1e9},   // original five-bin
    {0.175,0.205,1e9}                      // adopted two-bin
  };
  const char* nm[2]={"five-bin","two-bin"};

  for(int s=0;s<2;++s){
    const auto& e=schemes[s]; int nb=e.size()-1;
    std::vector<double> tot_all(nb,0),dia_all(nb,0),tot_g(nb,0),dia_g(nb,0);
    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
      int tb=-1,rb=-1;
      for(int b=0;b<nb;++b){ if(p_true>=e[b]&&p_true<e[b+1]) tb=b;
                             if(p_corr>=e[b]&&p_corr<e[b+1]) rb=b; }
      if(p_corr<e[0]) rb=0;              // lowest reco bin open below
      if(tb<0) continue;
      tot_all[tb]++; if(rb==tb) dia_all[tb]++;
      if(golden){ tot_g[tb]++; if(rb==tb) dia_g[tb]++; } }
    printf("\n  %s scheme\n    all pions   :",nm[s]);
    for(int b=0;b<nb;++b) printf(" %5.1f%%", tot_all[b]?100*dia_all[b]/tot_all[b]:0.);
    printf("\n    golden only :");
    for(int b=0;b<nb;++b) printf(" %5.1f%%", tot_g[b]?100*dia_g[b]/tot_g[b]:0.);
    printf("\n    golden frac :");
    for(int b=0;b<nb;++b) printf(" %5.1f%%", tot_all[b]?100*tot_g[b]/tot_all[b]:0.);
    printf("\n    events(all) :");
    for(int b=0;b<nb;++b) printf(" %6.0f", tot_all[b]);
    printf("\n");
  }
}
