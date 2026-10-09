// merge_univ_dirs.C -- add the per-file universe directories of one univmake output (src) to another (dst), for
// replacing a subset of the input files without re-running univmake on all of them (the reweighted detector-variation
// samples of the two-pion analysis, scripts/detvar_fluxw.py). Both files must come from the same bin configuration,
// so that they share the root directory name. Directories of dst with the same name are replaced; the cached totals
// ("total_*" subdirectories) of dst are deleted, so the next SystematicsCalculator rebuilds them from the per-file
// histograms. dst is modified in place: run it on a copy.
//   root -l -b -q 'macros/merge_univ_dirs.C("small_univmake.root","copy_of_univmake.root")'
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <iostream>
#include <set>
#include <string>

static void copy_dir( TDirectory* src, TDirectory* dst )
{
  std::set<std::string> done;
  for ( TObject* ko : *src->GetListOfKeys() ) {
    TKey* key = static_cast<TKey*>( ko );
    const std::string name = key->GetName();
    if ( !done.insert( name ).second ) continue;      // highest cycle only
    if ( std::string( key->GetClassName() ).find( "TDirectory" ) != std::string::npos ) {
      TDirectory* s = src->GetDirectory( name.c_str() );
      TDirectory* d = dst->mkdir( name.c_str(), s->GetTitle() );
      copy_dir( s, d );
    }
    else {
      TObject* o = key->ReadObj();
      dst->cd(); o->Write( name.c_str() );
      delete o;
    }
  }
}

int merge_univ_dirs( const char* src_name, const char* dst_name )
{
  TFile src( src_name, "read" ), dst( dst_name, "update" );
  if ( src.IsZombie() || dst.IsZombie() ) return 1;
  int n_root = 0, n_copied = 0;
  for ( TObject* ko : *src.GetListOfKeys() ) {
    TKey* key = static_cast<TKey*>( ko );
    if ( std::string( key->GetClassName() ).find( "TDirectory" ) == std::string::npos ) continue;
    TDirectory* sroot = src.GetDirectory( key->GetName() );
    TDirectory* droot = dst.GetDirectory( key->GetName() );
    if ( !droot ) { std::cerr << "no directory " << key->GetName() << " in " << dst_name << '\n'; return 1; }
    ++n_root;
    std::set<std::string> done;
    for ( TObject* ko2 : *droot->GetListOfKeys() ) {
      const std::string name = ko2->GetName();
      if ( name.rfind( "total_", 0 ) == 0 && done.insert( name ).second ) {
        droot->Delete( ( name + ";*" ).c_str() );
        std::cout << "deleted cached totals " << name << '\n';
      }
    }
    for ( TObject* ko2 : *sroot->GetListOfKeys() ) {
      TKey* k2 = static_cast<TKey*>( ko2 );
      const std::string name = k2->GetName();
      if ( std::string( k2->GetClassName() ).find( "TDirectory" ) == std::string::npos ) continue;
      if ( name.rfind( "total_", 0 ) == 0 || name.rfind( "+", 0 ) != 0 ) continue;   // per-file directories only
      if ( droot->GetDirectory( name.c_str() ) ) droot->Delete( ( name + ";*" ).c_str() );
      TDirectory* d = droot->mkdir( name.c_str() );
      copy_dir( sroot->GetDirectory( name.c_str() ), d );
      ++n_copied;
      std::cout << "copied " << name << '\n';
    }
  }
  dst.Write( "", TObject::kOverwrite );
  dst.Close();
  std::cout << n_copied << " per-file directories copied under " << n_root << " root directories\n";
  return n_root == 1 ? 0 : 1;
}
