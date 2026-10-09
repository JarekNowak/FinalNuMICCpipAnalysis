// mcs_cal_apply.C -- copy of a processed ntuple with the MCS momentum of the simulated exiting muon candidate calibrated
// to beam-on data (scripts/mcs_cal_derive.py: per cell of |dir x| and visible length of the candidate track, a map g of
// the MCS momentum). Mapped: candidates measured by MCS (reco > 0 and reco != range, as the impact study) and matched to
// a simulated particle (backtracked_pdg != 0); overlaid cosmic tracks are data already and stay as they are. The
// inclusive family (CC1mu1piXp_) changes candidate_muon_mom_reco only; the proton-tagged family (CC1mu1pi1p_, the w/
// files) also recomputes W_had and the TKI from the calibrated muon, as macros/mcs_map_fakedata_1p.C (delta phi_T reads
// only the muon direction and keeps its stored value). Every other object of the file is copied at its highest cycle.
//   root -l -b -e 'gSystem->AddIncludePath("-I$PWD/include")' -q 'macros/mcs_cal_apply.C+("in.root","out.root","map.txt")'
//   check = 1: no output; the proton-tagged recomputation is run on the uncalibrated muon and compared with the stored
//   values (they must agree)
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TVector3.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "XSecAnalyzer/NuMIBeamFrame.hh"
#include "../src/utils/STVTools.cxx"   // needs -I include (see usage)

namespace {
  inline double ssqrt( double x ) { return x > 0. ? std::sqrt( x ) : 0.; }
  struct Cell { double x0, x1, l0, l1; std::vector<double> g; };
  std::vector<Cell> cells;
  const double GSTEP = 0.01;
  bool load_map( const char* f ) {
    std::ifstream in( f ); std::string l;
    while ( std::getline( in, l ) ) {
      if ( l.empty() || l[0] == '#' ) continue;
      std::istringstream ss( l ); Cell c; ss >> c.x0 >> c.x1 >> c.l0 >> c.l1; double v;
      while ( ss >> v ) c.g.push_back( v );
      if ( c.g.size() < 2 ) return false;
      cells.push_back( c );
    }
    return !cells.empty();
  }
  const Cell* find_cell( double adx, double len ) {
    for ( const auto& c : cells ) if ( adx >= c.x0 && adx < c.x1 && len >= c.l0 && len < c.l1 ) return &c;
    return nullptr;
  }
  double gmap( const Cell& c, double p ) {
    const double pmax = GSTEP * ( c.g.size() - 1 );
    if ( p >= pmax ) return p * c.g.back() / pmax;
    const double u = p / GSTEP; const size_t i = static_cast<size_t>( u ); const double t = u - i;
    return c.g[i] + t * ( c.g[i + 1] - c.g[i] );
  }
}

int mcs_cal_apply( const char* in, const char* out, const char* mapfile, int check = 0 )
{
  if ( !load_map( mapfile ) ) { std::cerr << "cannot read the map " << mapfile << '\n'; return 1; }
  TFile fin( in, "read" );
  TTree* t = nullptr; fin.GetObject( "stv_tree", t );
  if ( !t ) { std::cerr << "no stv_tree in " << in << '\n'; return 1; }
  const bool oneP = t->GetBranch( "CC1mu1pi1p_Selected" ) != nullptr;
  const std::string S = oneP ? "CC1mu1pi1p_" : "CC1mu1piXp_";
  if ( !t->GetBranch( ( S + "CandidateMuonIndex" ).c_str() ) || !t->GetBranch( "backtracked_pdg" ) ) {
    std::cerr << "missing candidate or backtracked_pdg branches in " << in << '\n'; return 1; }
  t->SetBranchStatus( "*", check ? 0 : 1 );
  if ( check ) for ( const std::string& b : { S + "*", std::string( "trk_*" ), std::string( "backtracked_pdg" ),
                                              std::string( "reco_nu_vtx_sce_*" ) } ) t->SetBranchStatus( b.c_str(), 1 );

  int imu = -1, ipi = -1, ipr = -1;
  double reco = 0, range = 0, ppi = 0, Whad = 0, dat = 0, dpht = 0, dpt = 0, pn = 0;
  float vx = 0, vy = 0, vz = 0;
  std::vector<float> *dx = nullptr, *dy = nullptr, *dz = nullptr, *len = nullptr, *kep = nullptr;
  std::vector<int>* bt = nullptr;
  t->SetBranchAddress( ( S + "CandidateMuonIndex" ).c_str(), &imu );
  t->SetBranchAddress( ( S + "candidate_muon_mom_reco" ).c_str(), &reco );
  t->SetBranchAddress( ( S + "candidate_muon_mom_range" ).c_str(), &range );
  t->SetBranchAddress( "trk_dir_x_v", &dx ); t->SetBranchAddress( "trk_len_v", &len );
  t->SetBranchAddress( "backtracked_pdg", &bt );
  if ( oneP ) {
    t->SetBranchAddress( ( S + "CandidatePionIndex" ).c_str(), &ipi );
    t->SetBranchAddress( ( S + "CandidateProtonIndex" ).c_str(), &ipr );
    t->SetBranchAddress( ( S + "candidate_pion_mom_reco" ).c_str(), &ppi );
    t->SetBranchAddress( ( S + "W_had_reco" ).c_str(), &Whad );
    t->SetBranchAddress( ( S + "deltaAlphaT_reco" ).c_str(), &dat );
    t->SetBranchAddress( ( S + "deltaPhiT_reco" ).c_str(), &dpht );
    t->SetBranchAddress( ( S + "deltaPt_reco" ).c_str(), &dpt );
    t->SetBranchAddress( ( S + "pn_reco" ).c_str(), &pn );
    t->SetBranchAddress( "reco_nu_vtx_sce_x", &vx ); t->SetBranchAddress( "reco_nu_vtx_sce_y", &vy );
    t->SetBranchAddress( "reco_nu_vtx_sce_z", &vz );
    t->SetBranchAddress( "trk_dir_y_v", &dy ); t->SetBranchAddress( "trk_dir_z_v", &dz );
    t->SetBranchAddress( "trk_energy_proton_v", &kep );
  }

  const std::string tmp = std::string( out ) + ".tmp";
  TFile* fout = nullptr; TTree* o = nullptr;
  double o_reco = 0, o_Whad = 0, o_dat = 0, o_dpht = 0, o_dpt = 0, o_pn = 0;
  if ( !check ) {
    fout = new TFile( tmp.c_str(), "recreate" );
    std::set<std::string> done;
    for ( TObject* ko : *fin.GetListOfKeys() ) {
      const std::string name = static_cast<TKey*>( ko )->GetName();
      if ( name == "stv_tree" || !done.insert( name ).second ) continue;
      TObject* ob = fin.Get( name.c_str() ); fout->cd(); ob->Write( name.c_str() );
    }
    fout->cd();
    o = t->CloneTree( 0 );
    o->SetBranchAddress( ( S + "candidate_muon_mom_reco" ).c_str(), &o_reco );
    if ( oneP ) {
      o->SetBranchAddress( ( S + "W_had_reco" ).c_str(), &o_Whad );
      o->SetBranchAddress( ( S + "deltaAlphaT_reco" ).c_str(), &o_dat );
      o->SetBranchAddress( ( S + "deltaPhiT_reco" ).c_str(), &o_dpht );
      o->SetBranchAddress( ( S + "deltaPt_reco" ).c_str(), &o_dpt );
      o->SetBranchAddress( ( S + "pn_reco" ).c_str(), &o_pn );
    }
  }
  const Long64_t n = t->GetEntries();
  long nmcs = 0, nmap = 0, ncos = 0, ntki = 0, nbad = 0, nnocell = 0; double sr = 0., maxdev = 0.;
  for ( Long64_t i = 0; i < n; ++i ) {
    t->GetEntry( i );
    o_reco = reco; o_Whad = Whad; o_dat = dat; o_dpht = dpht; o_dpt = dpt; o_pn = pn;
    bool mapped = false;
    if ( imu >= 0 && reco > 0 && reco != range && imu < static_cast<int>( dx->size() ) ) {
      ++nmcs;
      if ( bt->at( imu ) != 0 ) {
        const Cell* c = find_cell( std::fabs( dx->at( imu ) ), len->at( imu ) );
        if ( c ) { o_reco = gmap( *c, reco ); sr += o_reco / reco; ++nmap; mapped = true; } else ++nnocell;
      } else ++ncos;
    }
    if ( oneP && ipr >= 0 && imu >= 0 && ipi >= 0 && ( mapped || check ) ) {
      const double pmu = check ? reco : o_reco;
      TVector3 mu( dx->at( imu ), dy->at( imu ), dz->at( imu ) ); mu.SetMag( pmu );
      const double Emu = ssqrt( pmu * pmu + MUON_MASS * MUON_MASS );
      TVector3 pi( dx->at( ipi ), dy->at( ipi ), dz->at( ipi ) ); pi.SetMag( ppi );
      const double Epi = ssqrt( ppi * ppi + PI_PLUS_MASS * PI_PLUS_MASS );
      const double KEp = kep->at( ipr ); const double pmom = ssqrt( KEp * KEp + 2. * PROTON_MASS * KEp );
      TVector3 pr( dx->at( ipr ), dy->at( ipr ), dz->at( ipr ) ); pr.SetMag( pmom );
      const double Epr = KEp + PROTON_MASS;
      const TVector3 nu = NuMIBeam::nu_dir_from_vertex( vx, vy, vz );
      const TVector3 mu_b = NuMIBeam::to_beam_frame( mu, nu ), had_b = NuMIBeam::to_beam_frame( pi + pr, nu );
      STVTools stv; stv.CalculateSTVs( mu_b, had_b, Emu, Epi + Epr, kOpt1 );
      o_dpt = stv.ReturnPt(); o_pn = stv.ReturnPn(); o_dat = stv.ReturnDeltaAlphaT();
      o_dpht = check ? stv.ReturnDeltaPhiT() : dpht;
      o_Whad = ssqrt( PROTON_MASS * PROTON_MASS + 2. * PROTON_MASS * ( stv.ReturnECal() - Emu ) - stv.ReturnQ2() );
      ++ntki;
      if ( check ) {
        const double d = std::max( { std::fabs( o_dpt - dpt ), std::fabs( o_pn - pn ), std::fabs( o_dat - dat ),
                                     std::fabs( o_dpht - dpht ), std::fabs( o_Whad - Whad ) } );
        if ( d > maxdev ) maxdev = d;
        if ( d > 1e-6 ) ++nbad;
      }
    }
    if ( o ) o->Fill();
  }
  std::printf( "%s %s: %lld entries, %ld MCS-measured candidates, %ld calibrated (mean ratio %.4f), %ld overlaid cosmic, "
               "%ld outside the cells, %ld W_had/TKI recomputed", oneP ? "1p" : "incl", in, n, nmcs, nmap,
               nmap ? sr / nmap : 0., ncos, nnocell, ntki );
  if ( check ) { std::printf( ", %ld deviate > 1e-6 (max %.3g)\n", nbad, maxdev ); return nbad ? 2 : 0; }
  std::printf( "\n" );
  o->Write( "", TObject::kOverwrite );
  fout->Close();
  if ( std::rename( tmp.c_str(), out ) != 0 ) { std::cerr << "rename failed\n"; return 1; }
  std::printf( "wrote %s\n", out );
  return nnocell ? 3 : 0;
}
