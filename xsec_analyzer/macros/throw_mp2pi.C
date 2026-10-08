// throw_mp2pi.C -- per-run Poisson CV pseudo-data of the two-pion one-bin total (Phase 4, item 3, of
// report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md), from the plan written by scripts/mp2pi_setup.py: for each beam-on
// entry the overlays AND the dirt of that run period, each event thrown Poisson(CV weight x scale) times with the
// scale data POT / distinct summed_pot of its type, weights set to one, summed_pot = data POT. The framework adds the
// beam-off to fake data itself. Same safe-weight rule as throw_perrun_*.C (non-finite, negative or > 30 -> 1).
//   root -l -b -q 'macros/throw_mp2pi.C("/data/uboone/processed/mp2pi/throw_plan.txt")'
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static void throw_one_mp2pi( const std::string& outfile, double dpot, int seed,
                             const std::vector<std::string>& in, const std::vector<double>& scale ) {
  gRandom->SetSeed( seed );
  TChain cin( "stv_tree" );
  for ( const auto& f : in ) cin.Add( f.c_str() );
  cin.SetBranchStatus( "*", 1 );
  for ( auto b : { "weight_All_UBGenie", "weight_ppfx_all", "weight_reint_all" } ) cin.SetBranchStatus( b, 0 );
  float tcv, pcv, nw;
  cin.SetBranchAddress( "tuned_cv_weight", &tcv ); cin.SetBranchAddress( "ppfx_cv_weight", &pcv );
  cin.SetBranchAddress( "normalisation_weight", &nw );
  TFile* out = new TFile( outfile.c_str(), "recreate" ); out->SetCompressionLevel( 1 );
  TTree* ot = cin.CloneTree( 0 );
  float one = 1.f;
  ot->SetBranchAddress( "tuned_cv_weight", &one ); ot->SetBranchAddress( "ppfx_cv_weight", &one );
  ot->SetBranchAddress( "normalisation_weight", &one );
  Long64_t N = cin.GetEntries(); long kept = 0;
  for ( Long64_t i = 0; i < N; ++i ) {
    cin.GetEntry( i );
    const double s = scale.at( cin.GetTreeNumber() );
    double cv = double( tcv ) * pcv * nw;
    if ( !std::isfinite( cv ) || cv < 0 || cv > 30 ) cv = 1.;
    const int nc = gRandom->Poisson( cv * s );
    for ( int c = 0; c < nc; ++c ) { one = 1.f; ot->Fill(); ++kept; }
  }
  ot->Write( "", TObject::kOverwrite );
  TParameter<float> sp( "summed_pot", float( dpot ) ); sp.Write( "summed_pot", TObject::kOverwrite );
  out->Close();
  printf( "  %s: %zu inputs, kept %ld events\n", outfile.c_str(), in.size(), kept );
}

void throw_mp2pi( const char* plan ) {
  std::ifstream fin( plan );
  std::string line, out; double dpot = 0.; int seed = 0;
  std::vector<std::string> in; std::vector<double> sc;
  auto flush = [&]() { if ( !out.empty() ) throw_one_mp2pi( out, dpot, seed, in, sc ); in.clear(); sc.clear(); };
  const std::string dir = std::string( plan ).substr( 0, std::string( plan ).rfind( '/' ) + 1 );
  while ( std::getline( fin, line ) ) {
    std::istringstream iss( line ); std::string key; iss >> key;
    if ( key == "OUT" ) { flush(); std::string name; iss >> name >> dpot >> seed; out = dir + name; }
    else if ( key == "IN" ) { std::string f; double s; iss >> f >> s; in.push_back( f ); sc.push_back( s ); }
  }
  flush();
}
