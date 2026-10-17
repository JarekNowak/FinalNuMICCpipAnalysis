// mp2pi_ens_throw.C -- one pseudo-experiment of the two-pion ensembles: every source file of macros/mp2pi_ens_source.C
// thrown event by event Poisson(ens_mu) times, with the CV weights set to one and summed_pot kept, into the member's own
// directory (one directory per member, so that no two jobs write the same file). The seed of file k of member m is a
// hash of (base, m, k) (splitmix64): consecutive seeds base + 1000 m + k, used first, gave Mersenne-Twister streams
// whose signal-region counts were sub-Poisson across members (variance / mean 0.86 over 1000 members); hashed seeds give
// 1.0. For an alternative model, ens_mu is multiplied by a per-event weight: a branch of the source (e.g. a unisim
// weight), or, when the argument starts with '=', a TTreeFormula expression of the source branches (stress models).
// mode 1 keeps the consecutive seeds of the first ensembles (for their reproduction).
//   root -l -b -q 'macros/mp2pi_ens_throw.C("<source dir>","<member dir>",<member>,<base seed>)'
#include <string>
#include <vector>
#include "TSystemDirectory.h"
#include "TTreeFormula.h"

static UInt_t ens_seed( ULong64_t x ) {
  x += 0x9e3779b97f4a7c15ULL; x = ( x ^ ( x >> 30 ) ) * 0xbf58476d1ce4e5b9ULL;
  x = ( x ^ ( x >> 27 ) ) * 0x94d049bb133111ebULL; x ^= x >> 31;
  const UInt_t s = UInt_t( x ^ ( x >> 32 ) ); return s ? s : 1u;
}

void mp2pi_ens_throw( const char* srcdir, const char* outdir, int member, int base = 20261016, const char* altw = "",
                      int mode = 0 ) {
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
    const std::string as( altw );
    const bool alt = as.size() > 0, formula = alt && as[0] == '=';
    TTreeFormula* tf = formula ? new TTreeFormula( "altw", as.substr( 1 ).c_str(), t ) : nullptr;
    if ( alt && !formula ) t->SetBranchAddress( altw, &aw );
    TFile* out = new TFile( ( std::string( outdir ) + "/" + files[k] + ".tmp" ).c_str(), "recreate" );
    out->SetCompressionLevel( 1 );
    TTree* ot = t->CloneTree( 0 );
    ot->SetBranchAddress( "tuned_cv_weight", &one ); ot->SetBranchAddress( "ppfx_cv_weight", &one );
    ot->SetBranchAddress( "normalisation_weight", &one );
    TRandom3 rng( mode == 1 ? UInt_t( base + 1000 * member + int( k ) )
                            : ens_seed( ( ULong64_t( UInt_t( base ) ) << 32 ) ^ ( ULong64_t( member ) << 8 ) ^ ULong64_t( k ) ) );
    long kept = 0;
    for ( Long64_t i = 0; i < t->GetEntries(); ++i ) {
      t->GetEntry( i );
      if ( tf ) aw = ( tf->GetNdata() > 0 ) ? float( tf->EvalInstance() ) : 1.f;   // an empty array: weight 1
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
