// golden_pion_compare.C — the with/without comparison for the golden-pion BDT.
//
// Evaluates the trained classifier on the HELD-OUT run (run_id==1, i.e. Run 2, which was
// not used for training) and reports, for a range of BDT cuts:
//   - selection efficiency and the golden fraction it achieves
//   - the p_pi migration matrix diagonal in the five analysis bins
//
// The migration diagonal is the metric that matters: the BNB CC1pi note requires it to
// exceed 0.68 for a bin to be usable, and it is what decides whether the top p_pi bins are
// measured or inferred from the response model.
//
//   usage: root -l -b -q 'macros/golden_pion_compare.C("nolen")'

void golden_pion_compare( const char* varset = "nolen" ) {

  TFile f( Form("/data/uboone/processed/validate/golden_train_%s.root", varset) );
  TTree* t = (TTree*) f.Get( "train" );
  if ( !t ) { printf( "  no training tree\n" ); return; }

  float wstd,wmean,wsep,bragg_pion,bragg_mip,bragg_mu,bragg_p,tscore,ndau,len;
  float p_true,p_range,p_corr; int golden,run_id;
  t->SetBranchAddress("wstd",&wstd);             t->SetBranchAddress("wmean",&wmean);
  t->SetBranchAddress("wsep",&wsep);             t->SetBranchAddress("bragg_pion",&bragg_pion);
  t->SetBranchAddress("bragg_mip",&bragg_mip);   t->SetBranchAddress("bragg_mu",&bragg_mu);
  t->SetBranchAddress("bragg_p",&bragg_p);       t->SetBranchAddress("tscore",&tscore);
  t->SetBranchAddress("ndau",&ndau);             t->SetBranchAddress("len",&len);
  t->SetBranchAddress("p_true",&p_true);         t->SetBranchAddress("p_range",&p_range);
  t->SetBranchAddress("p_corr",&p_corr);
  t->SetBranchAddress("golden",&golden);         t->SetBranchAddress("run_id",&run_id);

  TMVA::Reader reader( "!Color:!Silent" );
  reader.AddVariable( "wstd", &wstd );
  reader.AddVariable( "wmean", &wmean );
  reader.AddVariable( "wsep", &wsep );
  reader.AddVariable( "bragg_pion", &bragg_pion );
  reader.AddVariable( "bragg_mip", &bragg_mip );
  reader.AddVariable( "bragg_mu", &bragg_mu );
  reader.AddVariable( "bragg_p", &bragg_p );
  reader.AddVariable( "tscore", &tscore );
  reader.AddVariable( "ndau", &ndau );
  if ( std::string(varset) == "withlen" ) reader.AddVariable( "len", &len );
  reader.BookMVA( "BDT", Form("golden_%s/weights/GoldenPion_BDT.weights.xml", varset) );

  // analysis bins (true and reco share edges; reco bin 1 is everything below 0.250)
  const int NB = 5;
  double edge[NB+1] = { 0.175, 0.250, 0.320, 0.420, 0.550, 1e9 };
  auto bin_of = []( double p, double* e, int n ) {
    for ( int i = 0; i < n; ++i ) if ( p >= e[i] && p < e[i+1] ) return i;
    return -1;
  };

  const int NC = 7;
  double cuts[NC] = { -1e9, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15 };
  long   sel[NC] = {0}, selg[NC] = {0};
  long   diag[NC][NB] = {{0}}, tot[NC][NB] = {{0}};
  long   n_test = 0, n_test_gold = 0;

  Long64_t N = t->GetEntries();
  for ( Long64_t i = 0; i < N; ++i ) {
    t->GetEntry(i);
    if ( run_id != 1 ) continue;                 // held-out run only
    ++n_test; if ( golden ) ++n_test_gold;
    double resp = reader.EvaluateMVA( "BDT" );
    int tb = bin_of( p_true, edge, NB );
    int rb = bin_of( p_corr, edge, NB );
    if ( p_corr < edge[0] ) rb = 0;              // reco bin 1 is open below
    for ( int c = 0; c < NC; ++c ) {
      if ( resp < cuts[c] ) continue;
      ++sel[c]; if ( golden ) ++selg[c];
      if ( tb >= 0 ) { ++tot[c][tb]; if ( rb == tb ) ++diag[c][tb]; }
    }
  }

  printf( "\n  held-out sample (Run 2, not used in training): %ld tracks, golden fraction %.1f%%\n",
          n_test, 100.*n_test_gold/n_test );
  printf( "\n  BDT cut    eff     golden      migration diagonal per true bin [%%]\n" );
  printf( "                     fraction    b1     b2     b3     b4     b5\n" );
  for ( int c = 0; c < NC; ++c ) {
    printf( "  %s  %5.1f%%   %5.1f%%   ",
            c==0 ? "  none " : Form(" >%+.2f", cuts[c]),
            100.*sel[c]/n_test, 100.*selg[c]/std::max(1L,sel[c]) );
    for ( int b = 0; b < NB; ++b )
      printf( "%5.1f  ", tot[c][b] ? 100.*diag[c][b]/tot[c][b] : 0. );
    int npass = 0;
    for ( int b = 0; b < NB; ++b ) if ( tot[c][b] && 100.*diag[c][b]/tot[c][b] > 68. ) ++npass;
    printf( "  (%d/5 bins > 0.68)\n", npass );
  }
  printf( "\n  BNB reference: golden fraction 36.3%% -> 67.7%% at ~47%% of generic efficiency.\n" );
}
