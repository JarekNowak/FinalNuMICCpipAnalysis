// mp2pi_ens_throw.C -- one pseudo-experiment of the two-pion ensembles: every source file of macros/mp2pi_ens_source.C
// thrown event by event Poisson(ens_mu) times, with the CV weights set to one and summed_pot kept, into the member's own
// directory (one directory per member, so that no two jobs write the same file). The seed of file k of member m is
// base + 1000 m + k. For an alternative model, ens_mu is multiplied by the per-event weight branch given (e.g. a
// unisim weight stored in the source), so that the member is a throw of that model.
//   root -l -b -q 'macros/mp2pi_ens_throw.C("<source dir>","<member dir>",<member>,<base seed>)'
#include <string>
#include <vector>
#include "TSystemDirectory.h"

void mp2pi_ens_throw( const char* srcdir, const char* outdir, int member, int base = 20261016, const char* altw = "" ) {
  gSystem->mkdir( outdir, true );
  TSystemDirectory dir( "src", srcdir );
  std::vector<std::string> files;
  for ( TObject* o : *dir.GetListOfFiles() ) {
    const std::string n = o->GetName();
    if ( n.size() > 5 && n.substr( n.size() - 5 ) == ".root" ) files.push_back( n );
  }
  std::sort( files.begin(), files.end() );
  for ( size_t k = 0; k < files.size(); ++k ) {
    TFile fin( ( std::string( srcdir ) + "/" + files[k] ).c_str() );
    TTree* t = (TTree*)fin.Get( "stv_tree" );
    auto* sp = (TParameter<float>*)fin.Get( "summed_pot" );
    double mu = 0.; float aw = 1.f; float one = 1.f;
    t->SetBranchAddress( "ens_mu", &mu );
    const bool alt = std::string( altw ).size() > 0;
    if ( alt ) t->SetBranchAddress( altw, &aw );
    TFile* out = new TFile( ( std::string( outdir ) + "/" + files[k] + ".tmp" ).c_str(), "recreate" );
    out->SetCompressionLevel( 1 );
    TTree* ot = t->CloneTree( 0 );
    ot->SetBranchAddress( "tuned_cv_weight", &one ); ot->SetBranchAddress( "ppfx_cv_weight", &one );
    ot->SetBranchAddress( "normalisation_weight", &one );
    TRandom3 rng( base + 1000 * member + int( k ) );
    long kept = 0;
    for ( Long64_t i = 0; i < t->GetEntries(); ++i ) {
      t->GetEntry( i );
      double m = mu * ( alt ? ( ( std::isfinite( aw ) && aw >= 0 && aw <= 30 ) ? aw : 1. ) : 1. );
      const int nc = rng.Poisson( m );
      for ( int c = 0; c < nc; ++c ) { one = 1.f; ot->Fill(); ++kept; }
    }
    ot->Write( "", TObject::kOverwrite );
    if ( sp ) { TParameter<float> p( "summed_pot", sp->GetVal() ); p.Write( "summed_pot", TObject::kOverwrite ); }
    out->Close();
    std::rename( ( std::string( outdir ) + "/" + files[k] + ".tmp" ).c_str(), ( std::string( outdir ) + "/" + files[k] ).c_str() );
    printf( "  member %d %s: %ld events\n", member, files[k].c_str(), kept );
  }
}
