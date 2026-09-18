// throw_altmodel_fhc.C -- FHC fake data thrown from an ALTERNATIVE interaction model.
//
// The closure tests reported with every extraction throw the pseudo-data from the same simulation
// that builds the response, so they test the implementation and cannot test robustness to the model
// being wrong. This throws the pseudo-data from a reweighted model instead, while the response stays
// nominal: the bias that then appears in the unfolded result is the thing the nominal closure cannot
// see -- how much of the answer is the data and how much is the model assumed in unfolding it.
//
// Identical to throw_perrun_fhc.C except for one factor: the Poisson mean of each event is
// multiplied by its weight in the chosen universe. The weight applies to signal AND background, as
// a different interaction model would, so the test also covers a mismodelled background. Truth
// branches are carried through unchanged, so the framework's FakeData universe -- which is filled
// from these files in truth and reco (SystematicsCalculator.cxx:893) -- becomes the ALTERNATIVE
// model's truth, and the closure it reports is exactly the model-mismatch bias.
//
//   root -l -b -q 'macros/throw_altmodel_fhc.C("weight_All_UBGenie",350,"alt350",1)'
//   root -l -b -q 'macros/throw_altmodel_fhc.C("weight_Theta_Delta2Npi_UBGenie",0,"altdelta",1)'
#include <vector>
#include <string>

static double safe_w( double w ) {
  return ( !std::isfinite(w) || w < 0. || w > 30. ) ? 1.0 : w;
}

static void throw_alt_one( const char* infile, const char* outfile, double dpot, double mcpot,
                           const char* wbranch, int univ, int seed ) {
  const char* P = "/data/uboone/processed/";
  double potscale = dpot / mcpot;
  gRandom->SetSeed( seed );
  TChain cin( "stv_tree" );
  cin.Add( Form("%s%s", P, infile) );
  cin.SetBranchStatus( "*", 1 );
  // the multisim blocks are large; keep only the one being used
  for ( auto b : { "weight_All_UBGenie", "weight_ppfx_all", "weight_reint_all" } )
    if ( std::string(b) != std::string(wbranch) ) cin.SetBranchStatus( b, 0 );

  float tcv, pcv, nw;
  cin.SetBranchAddress( "tuned_cv_weight", &tcv );
  cin.SetBranchAddress( "ppfx_cv_weight", &pcv );
  cin.SetBranchAddress( "normalisation_weight", &nw );
  std::vector<double>* wv = nullptr;
  cin.SetBranchAddress( wbranch, &wv );

  TFile* out = new TFile( Form("%s%s", P, outfile), "recreate" );
  out->SetCompressionLevel( 1 );
  TTree* ot = cin.CloneTree( 0 );
  float one = 1.0f;
  ot->SetBranchAddress( "tuned_cv_weight", &one );
  ot->SetBranchAddress( "ppfx_cv_weight", &one );
  ot->SetBranchAddress( "normalisation_weight", &one );

  // UniverseMaker::apply_cv_correction_weights: a UBGenie universe weight REPLACES the tune, so the
  // event weight is w x ppfx x norm and NOT w x tune x ppfx x norm. Multiplying the tune back in
  // inflates the rate by ~68% and would have made this test a normalisation artefact rather than a
  // model variation. reint and SCC keep the tune; ppfx_all replaces the flux CV instead.
  const bool ubgenie = std::string( wbranch ).find( "UBGenie" ) != std::string::npos;
  const bool ppfx    = std::string( wbranch ).find( "ppfx" ) != std::string::npos;

  Long64_t N = cin.GetEntries();
  long kept = 0; double sum_cv = 0., sum_alt = 0.;
  for ( Long64_t i = 0; i < N; ++i ) {
    cin.GetEntry( i );
    double cv = safe_w( tcv * pcv * nw );
    double wu = ( wv && univ < (int)wv->size() ) ? safe_w( wv->at(univ) ) : 1.0;
    double alt = ubgenie ? wu * pcv * nw : ppfx ? wu * tcv * nw : wu * tcv * pcv * nw;
    sum_cv += cv; sum_alt += alt;
    int nc = gRandom->Poisson( alt * potscale );
    for ( int c = 0; c < nc; ++c ) { one = 1.0f; ot->Fill(); ++kept; }
  }
  ot->Write( "", TObject::kOverwrite );
  TParameter<float> sp( "summed_pot", (float)dpot );
  sp.Write( "summed_pot", TObject::kOverwrite );
  out->Close();
  printf( "  %-52s POTSCALE=%.5f  model rate x %.4f  kept=%ld\n",
          outfile, potscale, sum_alt / sum_cv, kept );
}

void throw_altmodel_fhc( const char* wbranch = "weight_All_UBGenie", int univ = 350,
                         const char* tag = "alt350", int seed = 1 ) {
  printf( "FHC alternative-model fake data: %s universe %d (tag %s, seed %d)\n",
          wbranch, univ, tag, seed );
  auto nm = [&]( const char* run ) { return Form("xsec-ana-fakedata_%s_fhc_%s.root", tag, run); };
  throw_alt_one( "xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root", nm("run1"),
                 3.283e20, 2.3282e21, wbranch, univ, seed );
  throw_alt_one( "xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root", nm("run2"),
                 1.268e20, 2.4934e21, wbranch, univ, seed + 1 );
  throw_alt_one( "xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root", nm("run4"),
                 2.075e20, 2.8335e21, wbranch, univ, seed + 2 );
  throw_alt_one( "xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root", nm("run5"),
                 2.231e20, 1.9300e21, wbranch, univ, seed + 3 );
}
