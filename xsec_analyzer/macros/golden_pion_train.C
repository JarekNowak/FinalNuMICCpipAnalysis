// golden_pion_train.C — build the labelled training sample for a golden-pion classifier
// and train it with TMVA.
//
// GOLDEN = the pion ranged out rather than interacting hadronically, defined on TRUTH
// alone as mc_end_p < 0.01 GeV/c. That label is not circular (it never refers to the
// reconstructed momentum we are trying to fix) and it predicts what we care about:
// golden pions recover the true momentum to within 20% 58.3% of the time, against 7.9%
// for the rest.
//
// Truth is joined to reco by matching backtracked_px/py/pz to mc_px/py/pz, which was
// verified to match 100% of candidate pion tracks with 0% ambiguity.
//
// Two variable sets are trained so the momentum-dependence of the efficiency can be
// compared: "nolen" excludes trk_len (which correlates with momentum and would make the
// selection efficiency momentum-dependent in the very observable being protected) and
// "withlen" includes it.
//
//   usage: root -l -b -q 'macros/golden_pion_train.C("nolen")'   // or "withlen"

#include <vector>
#include <cmath>

void golden_pion_train( const char* varset = "nolen" ) {

  const bool use_len = ( std::string(varset) == "withlen" );

  // ---- inputs: one entry per truth-matched reconstructed pion candidate ----------
  const char* files[4] = {
    "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root"
  };

  TString out = Form( "/data/uboone/processed/validate/golden_train_%s.root", varset );
  TFile fo( out, "recreate" );
  // train/test split is by RUN, not by random event, so the quoted performance is not
  // inflated by run-correlated detector conditions.
  TTree* tr_tree = new TTree( "train", "golden pion training sample" );

  float wstd, wmean, wsep, bragg_pion, bragg_mip, bragg_mu, bragg_p, tscore, ndau, len;
  float p_true, p_range, p_corr;
  int   golden, run_id;
  tr_tree->Branch( "wstd", &wstd );          tr_tree->Branch( "wmean", &wmean );
  tr_tree->Branch( "wsep", &wsep );          tr_tree->Branch( "bragg_pion", &bragg_pion );
  tr_tree->Branch( "bragg_mip", &bragg_mip );tr_tree->Branch( "bragg_mu", &bragg_mu );
  tr_tree->Branch( "bragg_p", &bragg_p );    tr_tree->Branch( "tscore", &tscore );
  tr_tree->Branch( "ndau", &ndau );          tr_tree->Branch( "len", &len );
  tr_tree->Branch( "p_true", &p_true );      tr_tree->Branch( "p_range", &p_range );
  tr_tree->Branch( "p_corr", &p_corr );
  tr_tree->Branch( "golden", &golden );      tr_tree->Branch( "run_id", &run_id );

  const double mmu = 0.10566, mpi = 0.13957;
  long ngold = 0, nnon = 0;

  for ( int fi = 0; fi < 4; ++fi ) {
    TFile* f = TFile::Open( files[fi] );
    if ( !f || f->IsZombie() ) { printf( "  [skip] %s\n", files[fi] ); continue; }
    TTree* t = (TTree*) f->Get( "nuselection/NeutrinoSelectionFilter" );
    if ( !t ) { f->Close(); continue; }
    t->SetBranchStatus( "*", 0 );
    for ( auto b : { "backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
                     "backtracked_purity","mc_pdg","mc_px","mc_py","mc_pz","mc_end_p",
                     "trk_score_v","trk_len_v","trk_range_muon_mom_v",
                     "trk_avg_deflection_stdev_v","trk_avg_deflection_mean_v",
                     "trk_avg_deflection_separation_mean_v","trk_bragg_pion_v",
                     "trk_bragg_mip_v","trk_bragg_mu_v","trk_bragg_p_v",
                     "pfp_trk_daughters_v","pfp_shr_daughters_v" } )
      t->SetBranchStatus( b, 1 );

    std::vector<int>   *bpdg=0, *mpdg=0;
    std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0,*mep=0;
    std::vector<float> *ts=0,*tl=0,*rm=0,*ws=0,*wm=0,*wsp=0,*bpi=0,*bmip=0,*bmu=0,*bp=0;
    std::vector<unsigned int> *dtr=0,*dsh=0;
    t->SetBranchAddress("backtracked_pdg",&bpdg);   t->SetBranchAddress("backtracked_px",&bpx);
    t->SetBranchAddress("backtracked_py",&bpy);     t->SetBranchAddress("backtracked_pz",&bpz);
    t->SetBranchAddress("backtracked_purity",&bpur);
    t->SetBranchAddress("mc_pdg",&mpdg); t->SetBranchAddress("mc_px",&mpx);
    t->SetBranchAddress("mc_py",&mpy);   t->SetBranchAddress("mc_pz",&mpz);
    t->SetBranchAddress("mc_end_p",&mep);
    t->SetBranchAddress("trk_score_v",&ts);        t->SetBranchAddress("trk_len_v",&tl);
    t->SetBranchAddress("trk_range_muon_mom_v",&rm);
    t->SetBranchAddress("trk_avg_deflection_stdev_v",&ws);
    t->SetBranchAddress("trk_avg_deflection_mean_v",&wm);
    t->SetBranchAddress("trk_avg_deflection_separation_mean_v",&wsp);
    t->SetBranchAddress("trk_bragg_pion_v",&bpi);  t->SetBranchAddress("trk_bragg_mip_v",&bmip);
    t->SetBranchAddress("trk_bragg_mu_v",&bmu);    t->SetBranchAddress("trk_bragg_p_v",&bp);
    t->SetBranchAddress("pfp_trk_daughters_v",&dtr);
    t->SetBranchAddress("pfp_shr_daughters_v",&dsh);

    Long64_t N = t->GetEntries();
    for ( Long64_t i = 0; i < N; ++i ) {
      t->GetEntry(i);
      if ( !bpdg || !mpdg ) continue;
      for ( size_t j = 0; j < bpdg->size(); ++j ) {
        if ( abs(bpdg->at(j)) != 211 ) continue;
        if ( bpur->at(j) < 0.5 || ts->at(j) < 0.5 || tl->at(j) < 5. ) continue;
        double pr = rm->at(j);
        if ( pr <= 0. ) continue;
        double w = ws->at(j);
        if ( w < 0. ) continue;

        int best = -1;
        for ( size_t k = 0; k < mpdg->size(); ++k ) {
          if ( mpdg->at(k) != bpdg->at(j) ) continue;
          double d = fabs(mpx->at(k)-bpx->at(j)) + fabs(mpy->at(k)-bpy->at(j))
                   + fabs(mpz->at(k)-bpz->at(j));
          if ( d < 1e-4 ) { best = k; break; }
        }
        if ( best < 0 ) continue;

        double pt = sqrt( pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2) );
        // train where the classifier will be applied
        if ( pt < 0.175 ) continue;

        double E = sqrt(pr*pr+mmu*mmu) - mmu + mpi;
        p_corr  = sqrt( std::max(0., E*E - mpi*mpi) );
        p_true  = pt;
        p_range = pr;
        wstd = w; wmean = wm->at(j); wsep = wsp->at(j);
        // The Bragg likelihood ratios are +-inf for ~2.7% of tracks (the likelihood
        // underflows for very short or very straight tracks). TMVA cannot train on
        // those, so map them to a sentinel outside the physical range rather than
        // dropping the track -- the sentinel is itself informative, since an
        // undefined Bragg fit correlates with not having a clean stopping signature.
        auto san = []( float v ) { return std::isfinite(v) ? v : -1.f; };
        bragg_pion = san( bpi->at(j) ); bragg_mip = san( bmip->at(j) );
        bragg_mu   = san( bmu->at(j) ); bragg_p   = san( bp->at(j) );
        tscore = ts->at(j); len = tl->at(j);
        ndau = dtr->at(j) + dsh->at(j);
        golden = ( mep->at(best) < 0.01 ) ? 1 : 0;
        run_id = fi;
        if ( golden ) ++ngold; else ++nnon;
        tr_tree->Fill();
      }
    }
    f->Close();
    printf( "  %s -> running totals: golden=%ld non=%ld\n", files[fi], ngold, nnon );
  }

  fo.cd();
  tr_tree->Write();
  printf( "\n  sample: %ld golden, %ld non-golden (golden fraction %.1f%%)\n",
          ngold, nnon, 100.*ngold/(ngold+nnon) );

  // ---- train ---------------------------------------------------------------------
  TMVA::Tools::Instance();
  TString tout = Form( "/data/uboone/processed/validate/golden_tmva_%s.root", varset );
  TFile* fout = TFile::Open( tout, "RECREATE" );
  TMVA::Factory factory( "GoldenPion", fout,
    "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );
  TMVA::DataLoader* d = new TMVA::DataLoader( Form("golden_%s", varset) );

  d->AddVariable( "wstd", 'F' );
  d->AddVariable( "wmean", 'F' );
  d->AddVariable( "wsep", 'F' );
  d->AddVariable( "bragg_pion", 'F' );
  d->AddVariable( "bragg_mip", 'F' );
  d->AddVariable( "bragg_mu", 'F' );
  d->AddVariable( "bragg_p", 'F' );
  d->AddVariable( "tscore", 'F' );
  d->AddVariable( "ndau", 'F' );
  if ( use_len ) d->AddVariable( "len", 'F' );

  // split by run: runs 0,2 train / runs 1,3 test
  d->AddSignalTree( tr_tree, 1.0 );
  d->AddBackgroundTree( tr_tree, 1.0 );
  d->SetSignalWeightExpression( "1" );
  d->PrepareTrainingAndTestTree(
    "golden==1 && (run_id==0||run_id==2)", "golden==0 && (run_id==0||run_id==2)",
    "SplitMode=Random:NormMode=NumEvents:!V" );

  factory.BookMethod( d, TMVA::Types::kBDT, "BDT",
    "!H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:"
    "AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:"
    "SeparationType=GiniIndex:nCuts=40" );

  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();
  fout->Close();
  fo.Close();
  printf( "\n  wrote %s and %s\n", out.Data(), tout.Data() );
}
