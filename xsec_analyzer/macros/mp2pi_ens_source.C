// mp2pi_ens_source.C -- sources of the pseudo-experiment ensembles of the two-pion one-bin total: for each pseudo-data
// file of a throw plan (scripts/mp2pi_setup.py, macros/throw_mp2pi.C), the events of its inputs (overlays and dirt of
// that run period) that enter a reco bin of configs/mp2pi_total_bin_config.txt (CC1mu2piBDT_evt_region >= 0) or the
// true signal bin (CC1mu2piBDT_MC_Signal), which are the only events through which a thrown pseudo-data file enters the
// framework: the reco bins (signal region and sidebands) and the truth of the pseudo-data. Each event carries ens_mu =
// CV weight x scale, the Poisson mean of its thrown count (same safe-weight rule as throw_mp2pi.C); the universe weight
// vectors are dropped, as in the throws, after the weights of two alternative models relative to the CV are stored
// (ens_w545: GENIE multisim universe 545, ens_wdelta: the Delta->N pi angular unisim; both replace the tune, as in
// scripts/mp_evt.py, so the ratio is the universe weight over tuned_cv_weight; 1 when the event carries no vector).
// The sources are thrown by macros/mp2pi_ens_throw.C.
//   root -l -b -q 'macros/mp2pi_ens_source.C("/data/uboone/processed/mp2pi_mcscal/throw_plan.txt","/data/uboone/processed/mp2pi_mcscal/ens/source")'
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static void source_one( const std::string& outfile, double dpot, const std::vector<std::string>& in,
                        const std::vector<double>& scale ) {
  TChain cin( "stv_tree" );
  for ( const auto& f : in ) cin.Add( f.c_str() );
  cin.SetBranchStatus( "*", 1 );
  for ( auto b : { "weight_All_UBGenie", "weight_ppfx_all", "weight_reint_all" } ) cin.SetBranchStatus( b, 0 );
  float tcv, pcv, nw; int region; bool sig;
  cin.SetBranchAddress( "tuned_cv_weight", &tcv ); cin.SetBranchAddress( "ppfx_cv_weight", &pcv );
  cin.SetBranchAddress( "normalisation_weight", &nw );
  cin.SetBranchAddress( "CC1mu2piBDT_evt_region", &region ); cin.SetBranchAddress( "CC1mu2piBDT_MC_Signal", &sig );
  // the alternative-model weights, read in lockstep from a second chain of the same files
  TChain cw( "stv_tree" );
  for ( const auto& f : in ) cw.Add( f.c_str() );
  cw.SetBranchStatus( "*", 0 );
  std::vector<double>* vg = nullptr; std::vector<double>* vd = nullptr;
  for ( auto b : { "weight_All_UBGenie", "weight_Theta_Delta2Npi_UBGenie" } ) cw.SetBranchStatus( b, 1 );
  cw.SetBranchAddress( "weight_All_UBGenie", &vg ); cw.SetBranchAddress( "weight_Theta_Delta2Npi_UBGenie", &vd );
  TFile* out = new TFile( ( outfile + ".tmp" ).c_str(), "recreate" ); out->SetCompressionLevel( 1 );
  TTree* ot = cin.CloneTree( 0 );
  double mu = 0.; float w545 = 1.f, wdelta = 1.f;
  ot->Branch( "ens_mu", &mu, "ens_mu/D" ); ot->Branch( "ens_w545", &w545, "ens_w545/F" ); ot->Branch( "ens_wdelta", &wdelta, "ens_wdelta/F" );
  auto ratio = []( const std::vector<double>* v, size_t u, float t ) -> float {
    if ( !v || v->size() <= u || !( t > 0.f ) ) return 1.f;
    const double r = v->at( u ) / t; return ( std::isfinite( r ) && r >= 0. && r <= 30. ) ? float( r ) : 1.f; };
  const Long64_t N = cin.GetEntries(); long kept = 0; double sum_all = 0., sum_kept = 0.;
  for ( Long64_t i = 0; i < N; ++i ) {
    cin.GetEntry( i );
    double cv = double( tcv ) * pcv * nw;
    if ( !std::isfinite( cv ) || cv < 0 || cv > 30 ) cv = 1.;
    mu = cv * scale.at( cin.GetTreeNumber() ); sum_all += mu;
    if ( region < 0 && !sig ) continue;
    // a file without the branch (the Delta unisim is absent from most dirt samples) keeps no stale value: weight 1
    cw.GetEntry( i ); TTree* ct = cw.GetTree();
    w545 = ( ct && ct->GetBranch( "weight_All_UBGenie" ) ) ? ratio( vg, 545, tcv ) : 1.f;
    wdelta = ( ct && ct->GetBranch( "weight_Theta_Delta2Npi_UBGenie" ) ) ? ratio( vd, 0, tcv ) : 1.f;
    sum_kept += mu; ot->Fill(); ++kept;
  }
  ot->Write( "", TObject::kOverwrite );
  TParameter<float> sp( "summed_pot", float( dpot ) ); sp.Write( "summed_pot", TObject::kOverwrite );
  out->Close();
  std::rename( ( outfile + ".tmp" ).c_str(), outfile.c_str() );
  printf( "  %s: %zu inputs, %lld entries, %ld kept (expected thrown count %.1f of %.1f)\n", outfile.c_str(), in.size(),
          N, kept, sum_kept, sum_all );
}

void mp2pi_ens_source( const char* plan, const char* outdir ) {
  std::ifstream fin( plan );
  std::string line, out; double dpot = 0.; int seed = 0;
  std::vector<std::string> in; std::vector<double> sc;
  gSystem->mkdir( outdir, true );
  auto flush = [&]() { if ( !out.empty() ) source_one( out, dpot, in, sc ); in.clear(); sc.clear(); };
  while ( std::getline( fin, line ) ) {
    std::istringstream iss( line ); std::string key; iss >> key;
    if ( key == "OUT" ) { flush(); std::string name; iss >> name >> dpot >> seed; out = std::string( outdir ) + "/" + name; }
    else if ( key == "IN" ) { std::string f; double s; iss >> f >> s; in.push_back( f ); sc.push_back( s ); }
  }
  flush();
}
