// dirt_perrun.h -- per-run dirt normalisation shared by the count-level macros (2026-10-01).
//
// Each run's dirt sample(s) are scaled by that run's data POT over the summed_pot of the
// sample(s) (the sum of the DISTINCT per-file values), exactly as SystematicsCalculator scales
// dirtMC. Runs, data POT and dirt samples are read from the release file-properties list, so
// nothing is hard-coded here. The returned scale does NOT contain the 0.65 dirt normalisation:
// ProcessNTuples stores it in normalisation_weight, so it enters with the central-value weight
// tune * ppfx_cv * normalisation, once.
//
// Until 2026-10-01 the macros scaled the Run-1 dirt sample alone by the data POT over the summed
// OVERLAY POT (0.081020 FHC, 0.071666 RHC, or 0.092402 before the Run-1 exposure correction),
// several of them also applying 0.65 a second time.
//
//   mode : "fhc" (runs 1, 2, 4, 5), "rhc" (runs 11-14) or "comb"
//   dir  : directory holding the dirt trees, under the release alias names (xsec-ana-dirt_*.root) or
//          the processed-file names (/data/uboone/processed/ carries the inclusive control-region
//          flags, .../w/ the proton-tagged); dirt_perrun_in takes one directory per run instead
//   skip : runs to leave out (e.g. {2, 12}: no Run-2 beam-on data)
// The POT is read from the file that is read. bkgfit/dirt.py is the Python counterpart.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TParameter.h"
#include "TSystem.h"

struct DirtSrc { std::string file; double scale; int run; };

// release alias -> name of the processed file it stands for (a fresh reprocess holds only the latter)
inline std::string dirt_orig_name( const std::string& alias )
{
  static const std::map< std::string, std::string > orig = {
    { "dirt_fhc_run1",  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot" },
    { "dirt_fhc_run2",  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot" },
    { "dirt_fhc_run4c", "numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana" },
    { "dirt_fhc_run4d", "numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana" },
    { "dirt_fhc_run5",  "run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns" },
    { "dirt_rhc_run1",  "neutrinoselection_filt_run3b_dirt_overlay" },
    { "dirt_rhc_run2",  "neutrinoselection_filt_run3b_dirt_overlay" },
    { "dirt_rhc_run3",  "neutrinoselection_filt_run3b_dirt_overlay" },
    { "dirt_rhc_run4a", "run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana" },
    { "dirt_rhc_run4b", "numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana" } };
  const std::string stem = alias.substr( 9, alias.size() - 9 - 5 );   // strip "xsec-ana-" and ".root"
  const auto it = orig.find( stem );
  return "xsec-ana-" + ( it != orig.end() ? it->second : stem ) + ".root";
}

inline std::vector<DirtSrc> dirt_perrun_in( const std::string& mode,
  const std::function< std::string( int ) >& dir_of, const std::vector<int>& skip = {},
  const std::string& fp_list = "" )
{
  const char* xa = gSystem->Getenv( "XSEC_ANALYZER_DIR" );
  const std::string conf = !fp_list.empty() ? fp_list
    : std::string( xa ? xa : "." ) + "/configs/file_properties_numi_comb_w.txt";
  std::ifstream in( conf );
  if ( !in ) { printf( "dirt_perrun: cannot read %s\n", conf.c_str() ); return {}; }
  std::map< int, double > data_pot;
  std::map< int, std::vector< std::string > > files;
  std::string line;
  while ( std::getline( in, line ) ) {
    const auto hash = line.find( '#' );
    if ( hash != std::string::npos ) line = line.substr( 0, hash );
    std::istringstream ss( line );
    std::string path, type; int run;
    if ( !( ss >> path >> run >> type ) ) continue;
    if ( type == "onBNB" ) { double trig = 0., pot = 0.; ss >> trig >> pot; data_pot[ run ] = pot; }
    else if ( type == "dirtMC" ) files[ run ].push_back( path.substr( path.find_last_of( '/' ) + 1 ) );
  }
  std::vector< DirtSrc > out;
  for ( const auto& rf : files ) {
    const int run = rf.first;
    const bool fhc = run < 10;
    if ( ( mode == "fhc" && !fhc ) || ( mode == "rhc" && fhc ) ) continue;
    if ( std::find( skip.begin(), skip.end(), run ) != skip.end() ) continue;
    const std::string dir = dir_of( run );
    std::vector< std::string > paths;
    for ( const auto& alias : rf.second ) {
      std::string f = dir + alias;
      if ( gSystem->AccessPathName( f.c_str() ) ) f = dir + dirt_orig_name( alias );   // true if NOT accessible
      paths.push_back( f );
    }
    std::vector< double > seen; double pot = 0.;
    for ( const auto& f : paths ) {
      TFile x( f.c_str() );
      auto* p = dynamic_cast< TParameter<float>* >( x.Get( "summed_pot" ) );
      if ( !p ) { printf( "dirt_perrun: no summed_pot in %s\n", f.c_str() ); continue; }
      const double v = p->GetVal();
      bool dup = false;
      for ( double s : seen ) if ( std::abs( s - v ) <= 1e-6 * std::abs( s ) ) dup = true;
      if ( !dup ) { seen.push_back( v ); pot += v; }
    }
    if ( pot <= 0. || !data_pot.count( run ) ) { printf( "dirt_perrun: run %d skipped\n", run ); continue; }
    for ( const auto& f : paths ) out.push_back( { f, data_pot.at( run ) / pot, run } );
  }
  return out;
}

inline std::vector<DirtSrc> dirt_perrun( const std::string& mode,
  const std::string& dir = "/data/uboone/processed/", const std::vector<int>& skip = {},
  const std::string& fp_list = "" )
{
  return dirt_perrun_in( mode, [dir]( int ) { return dir; }, skip, fp_list );
}
