// compare_reprocess.C -- gate for promoting a reprocessed ntuple over the live one.
//
// Every numeric leaf that exists in the OLD stv_tree is compared entry by entry with the NEW
// file (scalars and vector elements, bit-exact; NaN == NaN), plus the entry count and the
// summed_pot parameter (set_summed_pot*.C rewrites that in place, so a fresh reprocess could
// silently change the normalisation). Branches only in NEW are listed. Optionally the new
// pion-proton opening angle is checked against the friend tree of macros/proton_angle_friend.C
// on selected events (reco) and on events with a leading true pion and proton (truth).
//   root -l -b -q 'macros/compare_reprocess.C+("old.root","new.root","friend.pa.root")'
// Last line: [CMPREPRO] <new> OK|FAIL ...
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTreeFormula.h"
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

void compare_reprocess( const char* oldf, const char* newf, const char* friendf = "" ) {
  TFile fo( oldf ), fn( newf );
  TTree* to = (TTree*)fo.Get( "stv_tree" ), *tn = (TTree*)fn.Get( "stv_tree" );
  if ( !to || !tn ) { printf( "[CMPREPRO] %s FAIL missing tree\n", newf ); return; }
  bool ok = true;
  Long64_t N = to->GetEntries();
  if ( N != tn->GetEntries() ) { printf( "  entries differ %lld vs %lld\n", N, tn->GetEntries() ); ok = false; }

  auto pot = []( TFile& f ) -> double { auto* p = dynamic_cast<TParameter<float>*>( f.Get( "summed_pot" ) );
    if ( p ) return p->GetVal(); auto* q = dynamic_cast<TParameter<double>*>( f.Get( "summed_pot" ) ); return q ? q->GetVal() : -1.; };
  double po = pot( fo ), pn = pot( fn );
  if ( po != pn ) { printf( "  summed_pot differs: old %.6g new %.6g\n", po, pn ); ok = false; }

  // leaves: only numeric leaves are compared (vectors of numbers via TTreeFormula instances)
  std::vector<std::string> common; std::set<std::string> newonly;
  for ( auto* b : *tn->GetListOfBranches() ) newonly.insert( b->GetName() );
  for ( auto* b : *to->GetListOfBranches() ) {
    std::string n = b->GetName(); newonly.erase( n );
    if ( !tn->GetBranch( n.c_str() ) ) { printf( "  branch missing in new: %s\n", n.c_str() ); ok = false; continue; }
    common.push_back( n );
  }
  for ( auto& n : newonly ) printf( "  new branch: %s\n", n.c_str() );

  std::vector<TTreeFormula*> Fo, Fn; std::vector<std::string> names;
  for ( auto& n : common ) {
    TBranch* b = to->GetBranch( n.c_str() );
    std::string cls = b->GetClassName();
    // skip non-numeric objects (maps, vectors of vectors); vectors of arithmetic types are fine
    if ( !cls.empty() && cls.find( "vector<" ) != 0 ) continue;
    if ( cls.find( "vector<vector" ) == 0 || cls.find( "string" ) != std::string::npos ) continue;
    auto* a = new TTreeFormula( ( "o" + n ).c_str(), n.c_str(), to );
    auto* c = new TTreeFormula( ( "n" + n ).c_str(), n.c_str(), tn );
    if ( a->GetNdim() == 0 || c->GetNdim() == 0 ) { delete a; delete c; continue; }
    a->SetQuickLoad( false ); c->SetQuickLoad( false );
    Fo.push_back( a ); Fn.push_back( c ); names.push_back( n );
  }
  std::map<std::string, long> ndiff;
  for ( Long64_t i = 0; i < N; ++i ) {
    to->LoadTree( i ); tn->LoadTree( i );
    for ( size_t k = 0; k < Fo.size(); ++k ) {
      int na = Fo[k]->GetNdata(), nb = Fn[k]->GetNdata();
      if ( na != nb ) { ++ndiff[names[k]]; continue; }
      for ( int j = 0; j < na; ++j ) {
        double x = Fo[k]->EvalInstance( j ), y = Fn[k]->EvalInstance( j );
        if ( !( x == y || ( std::isnan( x ) && std::isnan( y ) ) ) ) { ++ndiff[names[k]]; break; }
      }
    }
  }
  for ( auto& d : ndiff ) { printf( "  DIFF %-45s %ld entries\n", d.first.c_str(), d.second ); ok = false; }
  printf( "  compared %zu branches over %lld entries\n", names.size(), N );

  if ( friendf && friendf[0] && tn->GetBranch( "CC1mu1pi1p_pi_pr_opening_angle_reco" ) ) {
    TFile ff( friendf ); TTree* pa = (TTree*)ff.Get( "pa" );
    if ( !pa || pa->GetEntries() != N ) { printf( "  friend unusable\n" ); ok = false; }
    else {
      tn->AddFriend( pa );
      TTreeFormula sel( "s", "CC1mu1pi1p_Selected", tn ), sig( "g", "CC1mu1pi1p_MC_Signal", tn ),
        rn( "rn", "CC1mu1pi1p_pi_pr_opening_angle_reco", tn ), rf( "rf", "pa.thpipr_reco", tn ),
        tt( "tt", "CC1mu1pi1p_pi_pr_opening_angle_true", tn ), tf( "tf", "pa.thpipr_true", tn );
      long nr = 0, br = 0, nt = 0, bt = 0;
      for ( Long64_t i = 0; i < N; ++i ) {
        tn->LoadTree( i );
        for ( auto* f : { &sel, &sig, &rn, &rf, &tt, &tf } ) f->GetNdata();
        if ( sel.EvalInstance() ) { ++nr; if ( std::abs( rn.EvalInstance() - rf.EvalInstance() ) > 1e-5 ) ++br; }
        if ( sig.EvalInstance() && tf.EvalInstance() > -9000. ) { ++nt; if ( std::abs( tt.EvalInstance() - tf.EvalInstance() ) > 1e-5 ) ++bt; }
      }
      printf( "  theta_pipr vs friend: reco %ld/%ld differ, truth(signal) %ld/%ld differ\n", br, nr, bt, nt );
      if ( br || bt ) ok = false;
    }
  }
  printf( "[CMPREPRO] %s %s\n", newf, ok ? "OK" : "FAIL" );
}
