// throw_altmodel_rhc.C -- RHC fake data thrown from an ALTERNATIVE interaction model.
//
// The RHC counterpart of throw_altmodel_fhc.C: pseudo-data thrown from a reweighted model while the
// response, efficiency and background stay nominal, so the closure the framework reports is the
// model-mismatch bias. RHC runs 3 and 4 are multi-file groups, chained and thrown at the group's
// data POT exactly as throw_perrun_rhc.C does.
//
// The UBGenie weight rule applies here too: a UBGenie universe weight REPLACES the MicroBooNE tune,
// so the event weight is w x ppfx x norm and not w x tune x ppfx x norm.
//
//   root -l -b -q 'macros/throw_altmodel_rhc.C("weight_All_UBGenie",545,"altgenie",1)'
//   root -l -b -q 'macros/throw_altmodel_rhc.C("weight_Theta_Delta2Npi_UBGenie",0,"altdelta",1)'
#include <vector>
#include <string>

static double safe_w_r( double w ) {
  return ( !std::isfinite(w) || w < 0. || w > 30. ) ? 1.0 : w;
}

static void throw_alt_group( std::vector<const char*> infiles, const char* outfile,
                             double dpot, double mcpot, const char* wbranch, int univ, int seed ) {
  const char* P = "/data/uboone/processed/";
  double potscale = dpot / mcpot;
  gRandom->SetSeed( seed );
  TChain cin( "stv_tree" );
  for ( auto f : infiles ) cin.Add( Form("%s%s", P, f) );
  cin.SetBranchStatus( "*", 1 );
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

  const bool ubgenie = std::string( wbranch ).find( "UBGenie" ) != std::string::npos;
  const bool isppfx  = std::string( wbranch ).find( "ppfx" ) != std::string::npos;

  Long64_t N = cin.GetEntries();
  long kept = 0; double sum_cv = 0., sum_alt = 0.;
  for ( Long64_t i = 0; i < N; ++i ) {
    cin.GetEntry( i );
    double cv = safe_w_r( tcv * pcv * nw );
    double wu = ( wv && univ < (int)wv->size() ) ? safe_w_r( wv->at(univ) ) : 1.0;
    double alt = ubgenie ? wu * pcv * nw : isppfx ? wu * tcv * nw : wu * tcv * pcv * nw;
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

void throw_altmodel_rhc( const char* wbranch = "weight_All_UBGenie", int univ = 545,
                         const char* tag = "altgenie", int seed = 1 ) {
  printf( "RHC alternative-model fake data: %s universe %d (tag %s, seed %d)\n",
          wbranch, univ, tag, seed );
  auto nm = [&]( const char* run ) { return Form("xsec-ana-fakedata_%s_rhc_%s.root", tag, run); };
  throw_alt_group( { "xsec-ana-Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root" }, nm("run1"),
                   0.6053e20, 8.9972e20, wbranch, univ, seed );
  throw_alt_group( { "xsec-ana-Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root" }, nm("run2"),
                   2.591e20, 5.7865e21, wbranch, univ, seed + 1 );
  throw_alt_group( { "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root",
                     "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root",
                     "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root",
                     "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root",
                     "xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root" }, nm("run3"),
                   5.003e20, 5.5185e21, wbranch, univ, seed + 2 );
  throw_alt_group( { "xsec-ana-Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root",
                     "xsec-ana-Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root",
                     "xsec-ana-Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root" }, nm("run4"),
                   2.883e20, 3.2586e21, wbranch, univ, seed + 3 );
}
