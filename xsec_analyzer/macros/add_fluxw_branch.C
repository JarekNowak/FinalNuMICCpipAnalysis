// add_fluxw_branch.C -- copy of a processed detector-variation ntuple with the per-event new / old NuMI flux weight
// (scripts/detvar_fluxw.py writes it as raw float32, one value per stv_tree entry) as the branch numi_flux_rw_weight,
// which UniverseMaker applies in every universe. The weight_* universe branches are left out: the detector variations
// enter only through the unweighted universe, and without weight_TunedCentralValue_UBGenie univmake builds only that.
// Every other object of the file (summed_pot, cut-flow histograms) is copied at its highest cycle.
//   root -l -b -q 'macros/add_fluxw_branch.C("in.root","in.fluxw.f32","out.root")'
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

int add_fluxw_branch( const char* in, const char* wfile, const char* out )
{
  std::ifstream wf( wfile, std::ios::binary | std::ios::ate );
  if ( !wf ) { std::cerr << "cannot read " << wfile << '\n'; return 1; }
  const std::streamsize nbytes = wf.tellg(); wf.seekg( 0 );
  std::vector<float> w( nbytes / sizeof( float ) );
  wf.read( reinterpret_cast<char*>( w.data() ), nbytes );

  TFile fin( in, "read" );
  TTree* t = nullptr; fin.GetObject( "stv_tree", t );
  if ( !t ) { std::cerr << "no stv_tree in " << in << '\n'; return 1; }
  const Long64_t n = t->GetEntries();
  if ( n != static_cast<Long64_t>( w.size() ) ) {
    std::cerr << "entries " << n << " but " << w.size() << " weights\n"; return 1;
  }

  const std::string tmp = std::string( out ) + ".tmp";
  TFile fout( tmp.c_str(), "recreate" );
  std::set<std::string> done;
  for ( TObject* ko : *fin.GetListOfKeys() ) {
    const std::string name = static_cast<TKey*>( ko )->GetName();
    if ( name == "stv_tree" || !done.insert( name ).second ) continue;
    TObject* o = fin.Get( name.c_str() );
    fout.cd(); o->Write( name.c_str() );
  }

  t->SetBranchStatus( "*", 1 );
  t->SetBranchStatus( "weight_*", 0 );
  fout.cd();
  TTree* o = t->CloneTree( 0 );
  float fw = 1.f;
  o->Branch( "numi_flux_rw_weight", &fw, "numi_flux_rw_weight/F" );
  double sum = 0.;
  for ( Long64_t i = 0; i < n; ++i ) {
    t->GetEntry( i ); fw = w[i]; sum += fw; o->Fill();
  }
  o->Write( "", TObject::kOverwrite );
  fout.Close();
  if ( std::rename( tmp.c_str(), out ) != 0 ) { std::cerr << "rename failed\n"; return 1; }
  std::cout << "wrote " << out << ": " << n << " entries, mean weight " << sum / n << '\n';
  return 0;
}
