// dump_pid_tracks.C -- training sample for the multi-pion particle classifier (Phase 1 of
// report/MULTIPION_BNB_ADAPTATION_PLAN.md).
//
// Reads a raw NuMI PeLEE ntuple (nuselection/NeutrinoSelectionFilter) and writes one row per primary
// track-like PFParticle of the CC1mu1piXp candidate pool, with the reconstructed PID features, the
// backtracked truth and the scores of the BDTs the selection uses today, evaluated on the same track
// with the same inputs as CC1mu1piXp.cxx. Reading the raw ntuple directly avoids a reprocessing: the
// processed files do not carry the chi2, PIDA, trunk dE/dx, calorimetry or end-spacepoint vectors.
//
// Event preselection (as the selection): one neutrino slice and the SCE-corrected reco vertex in the
// reco fiducial volume outside the dead z region. The topological score is stored, not cut on.
// Track pool (as the pion-candidate loop): generation 2, track score >= 0.3, LLR score in (-1, 2),
// Bragg p/mu/mip in (0, 500), length > 0. Training cuts (purity, track score, length) are applied later.
//
//   root -l -b -q 'macros/mp_pid/dump_pid_tracks.C+("in.root","out.root",1,0,0)'
//   arguments: input, output, run period (1-5), horn (0 FHC, 1 RHC), max events (0 = all)
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TMVA/Reader.h"

namespace {
template <class T> struct Vec { std::vector<T>* p = nullptr; bool ok = false; };

template <class T> void bind( TTree* t, const char* name, Vec<T>& v ) {
  if ( !t->GetBranch( name ) ) { printf( "  (branch %s absent)\n", name ); return; }
  t->SetBranchStatus( name, 1 ); t->SetBranchAddress( name, &v.p ); v.ok = true;
}
template <class T> void bind_scalar( TTree* t, const char* name, T& x ) {
  if ( !t->GetBranch( name ) ) { printf( "  (branch %s absent)\n", name ); return; }
  t->SetBranchStatus( name, 1 ); t->SetBranchAddress( name, &x );
}
template <class T> float at( const Vec<T>& v, size_t i, float dflt = -9999.f ) {
  return ( v.ok && v.p && i < v.p->size() ) ? float( v.p->at( i ) ) : dflt;
}
bool inside( double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1 ) {
  return x > x0 && x < x1 && y > y0 && y < y1 && z > z0 && z < z1;
}
// reco fiducial volume of CC1mu1piXp (define_reco_FV) and the z dead region of the vertex FV
bool in_vertex_fv( double x, double y, double z ) {
  if ( !inside( x, y, z, 10., 246., -101., 101., 10., 986. ) ) return false;
  return !( z > 675.1 && z < 775.1 );
}
// containment volume of CC1mu1piXp (define_containment_FV)
bool contained( double x, double y, double z ) { return inside( x, y, z, 2.0, 254.35, -113.53, 107.47, 2.1, 1034.9 ); }
float safe_weight( double w ) { return ( std::isfinite( w ) && w >= 0. && w <= 30. ) ? float( w ) : 1.f; }
}

void dump_pid_tracks( const char* in_path, const char* out_path, int period, int horn, long max_events = 0 ) {
  TFile fin( in_path );
  TTree* t = (TTree*)fin.Get( "nuselection/NeutrinoSelectionFilter" );
  if ( !t ) { printf( "no nuselection/NeutrinoSelectionFilter in %s\n", in_path ); return; }
  t->SetBranchStatus( "*", 0 );

  int run = 0, sub = 0, evt = 0, nslice = 0, nu_pdg = 0, ccnc = -1;
  float topo = 0.f, vx = 0.f, vy = 0.f, vz = 0.f, w_tune = 1.f, w_ppfx = 1.f;
  bind_scalar( t, "run", run ); bind_scalar( t, "sub", sub ); bind_scalar( t, "evt", evt );
  bind_scalar( t, "nslice", nslice ); bind_scalar( t, "topological_score", topo );
  bind_scalar( t, "reco_nu_vtx_sce_x", vx ); bind_scalar( t, "reco_nu_vtx_sce_y", vy ); bind_scalar( t, "reco_nu_vtx_sce_z", vz );
  bind_scalar( t, "weightSplineTimesTune", w_tune ); bind_scalar( t, "ppfx_cv", w_ppfx );
  bind_scalar( t, "nu_pdg", nu_pdg ); bind_scalar( t, "ccnc", ccnc );

  Vec<unsigned int> gen, ntrkd, nshrd, ndesc;
  Vec<float> ts, len, sx, sy, sz, ex, ey, ez, llr, llr_u, llr_v, llr_y;
  Vec<float> bp, bmu, bmip, bpi, bp_u, bp_v, bmu_u, bmu_v, bpi_u, bpi_v, bmip_u, bmip_v;
  Vec<float> pida, pida_u, pida_v, chipr, chimu, chipi, chika, chipr_u, chipr_v, chimu_u, chimu_v, chipi_u, chipi_v;
  Vec<float> trunk_y, trunk_u, trunk_v, trunkrr_y, trunkrr_u, trunkrr_v, calo_y, calo_u, calo_v;
  Vec<float> defl_mean, defl_std, defl_sep, mcs, range_mu, e_p;
  Vec<int> endsp, nh_u, nh_v, nh_y, ph_u, ph_v, ph_y, bt_pdg;
  Vec<bool> fwd_p, fwd_mu, fwd_pi;
  Vec<float> bt_pur, bt_comp, bt_px, bt_py, bt_pz;
  bind( t, "pfp_generation_v", gen ); bind( t, "pfp_trk_daughters_v", ntrkd ); bind( t, "pfp_shr_daughters_v", nshrd );
  bind( t, "pfp_n_descendents_v", ndesc );
  bind( t, "trk_score_v", ts ); bind( t, "trk_len_v", len );
  bind( t, "trk_sce_start_x_v", sx ); bind( t, "trk_sce_start_y_v", sy ); bind( t, "trk_sce_start_z_v", sz );
  bind( t, "trk_sce_end_x_v", ex ); bind( t, "trk_sce_end_y_v", ey ); bind( t, "trk_sce_end_z_v", ez );
  bind( t, "trk_llr_pid_score_v", llr ); bind( t, "trk_llr_pid_u_v", llr_u ); bind( t, "trk_llr_pid_v_v", llr_v );
  bind( t, "trk_llr_pid_y_v", llr_y );
  bind( t, "trk_bragg_p_v", bp ); bind( t, "trk_bragg_mu_v", bmu ); bind( t, "trk_bragg_mip_v", bmip );
  bind( t, "trk_bragg_pion_v", bpi );
  bind( t, "trk_bragg_p_u_v", bp_u ); bind( t, "trk_bragg_p_v_v", bp_v ); bind( t, "trk_bragg_mu_u_v", bmu_u );
  bind( t, "trk_bragg_mu_v_v", bmu_v ); bind( t, "trk_bragg_pion_u_v", bpi_u ); bind( t, "trk_bragg_pion_v_v", bpi_v );
  bind( t, "trk_bragg_mip_u_v", bmip_u ); bind( t, "trk_bragg_mip_v_v", bmip_v );
  bind( t, "trk_bragg_p_fwd_preferred_v", fwd_p ); bind( t, "trk_bragg_mu_fwd_preferred_v", fwd_mu );
  bind( t, "trk_bragg_pion_fwd_preferred_v", fwd_pi );
  bind( t, "trk_pida_v", pida ); bind( t, "trk_pida_u_v", pida_u ); bind( t, "trk_pida_v_v", pida_v );
  bind( t, "trk_pid_chipr_v", chipr ); bind( t, "trk_pid_chimu_v", chimu ); bind( t, "trk_pid_chipi_v", chipi );
  bind( t, "trk_pid_chika_v", chika ); bind( t, "trk_pid_chipr_u_v", chipr_u ); bind( t, "trk_pid_chipr_v_v", chipr_v );
  bind( t, "trk_pid_chimu_u_v", chimu_u ); bind( t, "trk_pid_chimu_v_v", chimu_v ); bind( t, "trk_pid_chipi_u_v", chipi_u );
  bind( t, "trk_pid_chipi_v_v", chipi_v );
  bind( t, "trk_trunk_dEdx_y_v", trunk_y ); bind( t, "trk_trunk_dEdx_u_v", trunk_u ); bind( t, "trk_trunk_dEdx_v_v", trunk_v );
  bind( t, "trk_trunk_rr_dEdx_y_v", trunkrr_y ); bind( t, "trk_trunk_rr_dEdx_u_v", trunkrr_u );
  bind( t, "trk_trunk_rr_dEdx_v_v", trunkrr_v );
  bind( t, "trk_calo_energy_y_v", calo_y ); bind( t, "trk_calo_energy_u_v", calo_u ); bind( t, "trk_calo_energy_v_v", calo_v );
  bind( t, "trk_avg_deflection_mean_v", defl_mean ); bind( t, "trk_avg_deflection_stdev_v", defl_std );
  bind( t, "trk_avg_deflection_separation_mean_v", defl_sep );
  bind( t, "trk_mcs_muon_mom_v", mcs ); bind( t, "trk_range_muon_mom_v", range_mu ); bind( t, "trk_energy_proton_v", e_p );
  bind( t, "trk_end_spacepoints_v", endsp );
  bind( t, "trk_nhits_u_v", nh_u ); bind( t, "trk_nhits_v_v", nh_v ); bind( t, "trk_nhits_y_v", nh_y );
  bind( t, "pfnplanehits_U", ph_u ); bind( t, "pfnplanehits_V", ph_v ); bind( t, "pfnplanehits_Y", ph_y );
  bind( t, "backtracked_pdg", bt_pdg ); bind( t, "backtracked_purity", bt_pur ); bind( t, "backtracked_completeness", bt_comp );
  bind( t, "backtracked_px", bt_px ); bind( t, "backtracked_py", bt_py ); bind( t, "backtracked_pz", bt_pz );

  // ---- the BDTs of CC1mu1piXp, with the same inputs and order ----
  const char* xa = gSystem->Getenv( "XSEC_ANALYZER_DIR" );
  const std::string bdt_dir = std::string( xa ? xa : "." ) + "/../booster_decision_tree/";
  float r_bp, r_bmu, r_bmip, r_llr, r_ts, r_len, r_ex, r_ey, r_ez, m_llr, m_bp, m_bmu, m_bmip, m_bpi, m_ts, m_len, m_dist;
  TMVA::Reader rd_mip( "!Color:Silent" ), rd_mu( "!Color:Silent" ), rd_pi( "!Color:Silent" );
  for ( TMVA::Reader* r : { &rd_mip, &rd_mu, &rd_pi } ) {
    r->AddVariable( "trk_bragg_p_v", &r_bp ); r->AddVariable( "trk_bragg_mu_v", &r_bmu );
    r->AddVariable( "trk_bragg_mip_v", &r_bmip ); r->AddVariable( "trk_llr_pid_score_v", &r_llr );
    r->AddVariable( "trk_score_v", &r_ts );
    if ( r == &rd_mu ) r->AddVariable( "trk_len_v", &r_len );
    r->AddVariable( "trk_sce_end_x_v", &r_ex ); r->AddVariable( "trk_sce_end_y_v", &r_ey ); r->AddVariable( "trk_sce_end_z_v", &r_ez );
  }
  TMVA::Reader rd_mp( "!Color:Silent" ), rd_mp_soft( "!Color:Silent" ), rd_mp_hard( "!Color:Silent" );
  for ( TMVA::Reader* r : { &rd_mp, &rd_mp_soft, &rd_mp_hard } ) {
    r->AddVariable( "llr", &m_llr ); r->AddVariable( "bragg_p", &m_bp ); r->AddVariable( "bragg_mu", &m_bmu );
    r->AddVariable( "bragg_mip", &m_bmip ); r->AddVariable( "bragg_pion", &m_bpi ); r->AddVariable( "trk_score", &m_ts );
    r->AddVariable( "length", &m_len ); r->AddVariable( "dist", &m_dist );
  }
  auto book = [&]( TMVA::Reader& r, const char* ds ) {
    r.BookMVA( "BDT", ( bdt_dir + ds + "/weights/TMVAClassification_BDT.weights.xml" ).c_str() ); };
  book( rd_mip, "dataset_MIP_BDT_no_len" ); book( rd_mu, "dataset_muon_BDT" ); book( rd_pi, "dataset_pion_BDT_no_len" );
  book( rd_mp, "mp_pion_bdt" ); book( rd_mp_soft, "mp_pion_bdt_soft" ); book( rd_mp_hard, "mp_pion_bdt_hard" );

  // ---- output ----
  TFile fout( out_path, "RECREATE" );
  TTree o( "trk", "primary track-like PFParticles of the CC1mu1piXp candidate pool" );
  int o_period = period, o_horn = horn, o_run, o_sub, o_evt, o_pfp, o_pdg, o_label, o_charge, o_mucand, o_ccnc, o_nupdg;
  Long64_t o_entry; int o_nprimtrk, o_nprimshr;
  float o_w, o_topo, o_pur, o_comp, o_truep;
  o.Branch( "period", &o_period ); o.Branch( "horn", &o_horn ); o.Branch( "run", &o_run ); o.Branch( "sub", &o_sub );
  o.Branch( "evt", &o_evt ); o.Branch( "entry", &o_entry ); o.Branch( "pfp", &o_pfp ); o.Branch( "w_cv", &o_w ); o.Branch( "topo", &o_topo );
  o.Branch( "nu_pdg", &o_nupdg ); o.Branch( "ccnc", &o_ccnc );
  o.Branch( "pdg", &o_pdg ); o.Branch( "label", &o_label ); o.Branch( "charge", &o_charge );
  o.Branch( "purity", &o_pur ); o.Branch( "completeness", &o_comp ); o.Branch( "true_p", &o_truep );
  o.Branch( "is_mu_cand", &o_mucand ); o.Branch( "n_prim_trk", &o_nprimtrk ); o.Branch( "n_prim_shr", &o_nprimshr );
  // features: name -> value, filled per track in this order
  std::vector<std::string> fnames = {
    "len", "ts", "dist", "contained", "start_contained", "llr", "llr_u", "llr_v", "llr_y",
    "bragg_p", "bragg_mu", "bragg_mip", "bragg_pion", "bragg_p_u", "bragg_p_v", "bragg_mu_u", "bragg_mu_v",
    "bragg_pion_u", "bragg_pion_v", "bragg_mip_u", "bragg_mip_v", "fwd_p", "fwd_mu", "fwd_pion",
    "pida", "pida_u", "pida_v", "chipr", "chimu", "chipi", "chika", "chipr_u", "chipr_v", "chimu_u", "chimu_v",
    "chipi_u", "chipi_v", "trunk_dedx_y", "trunk_dedx_u", "trunk_dedx_v", "trunk_rr_dedx_y", "trunk_rr_dedx_u",
    "trunk_rr_dedx_v", "calo_e_y", "calo_e_u", "calo_e_v", "defl_mean", "defl_stdev", "defl_sep_mean",
    "mcs_mom", "range_mom_mu", "mcs_over_range", "e_proton", "end_sp", "nhits_u", "nhits_v", "nhits_y",
    "planehits_u", "planehits_v", "planehits_y", "n_trk_daughters", "n_shr_daughters", "n_descendents",
    "s_mip", "s_muon", "s_pion", "s_mppi", "s_mppi_soft", "s_mppi_hard" };
  std::vector<float> f( fnames.size() );
  for ( size_t k = 0; k < fnames.size(); ++k ) o.Branch( fnames[k].c_str(), &f[k] );

  const long n = ( max_events > 0 && max_events < t->GetEntries() ) ? max_events : t->GetEntries();
  long n_ev_pass = 0, n_trk = 0;
  for ( long ie = 0; ie < n; ++ie ) {
    t->GetEntry( ie );
    if ( nslice != 1 || !in_vertex_fv( vx, vy, vz ) || !gen.ok || !gen.p ) continue;
    ++n_ev_pass;
    const size_t np = gen.p->size();
    // pool of the pion-candidate loop, with the per-track BDT scores
    struct T { size_t i; float s[6]; float dist; };
    std::vector<T> pool;
    for ( size_t i = 0; i < np; ++i ) {
      if ( gen.p->at( i ) != 2u ) continue;
      const float tsc = at( ts, i ), l = at( llr, i ), b1 = at( bp, i ), b2 = at( bmu, i ), b3 = at( bmip, i ), ln = at( len, i );
      if ( !( tsc >= 0.3f && l > -1.f && l < 2.f && b1 > 0.f && b1 < 500.f && b2 > 0.f && b2 < 500.f && b3 > 0.f && b3 < 500.f
              && ln > 0.f && ln < 1e6f ) ) continue;
      T tr; tr.i = i;
      const double dx = at( sx, i ) - vx, dy = at( sy, i ) - vy, dz = at( sz, i ) - vz;
      tr.dist = float( std::sqrt( dx * dx + dy * dy + dz * dz ) );
      r_bp = b1; r_bmu = b2; r_bmip = b3; r_llr = l; r_ts = tsc; r_len = ln; r_ex = at( ex, i ); r_ey = at( ey, i ); r_ez = at( ez, i );
      tr.s[0] = rd_mip.EvaluateMVA( "BDT" ); tr.s[1] = rd_mu.EvaluateMVA( "BDT" ); tr.s[2] = rd_pi.EvaluateMVA( "BDT" );
      m_llr = l; m_bp = b1; m_bmu = b2; m_bmip = b3; m_bpi = at( bpi, i, 1.f ); m_ts = tsc; m_len = ln; m_dist = tr.dist;
      tr.s[3] = rd_mp.EvaluateMVA( "BDT" ); tr.s[4] = rd_mp_soft.EvaluateMVA( "BDT" ); tr.s[5] = rd_mp_hard.EvaluateMVA( "BDT" );
      pool.push_back( tr );
    }
    // muon candidate exactly as CC1mu1piXp: only in events with >= 2 primaries of track score > 0.5,
    // the highest muon-BDT score among the qualifying tracks
    int n_prim_trk = 0, n_prim_shr = 0;
    for ( size_t i = 0; i < np; ++i ) if ( gen.p->at( i ) == 2u ) ( at( ts, i ) > 0.5f ? n_prim_trk : n_prim_shr )++;
    int mu_idx = -1; float best = -1e30f;
    if ( n_prim_trk >= 2 ) for ( const auto& tr : pool ) {
      const size_t i = tr.i;
      if ( at( ts, i ) >= 0.5f && at( llr, i ) > 0.2f && at( len, i ) > 10.f && at( ts, i ) > 0.8f && tr.dist <= 4.f && tr.s[0] >= -0.1f )
        if ( mu_idx < 0 || tr.s[1] > best ) { mu_idx = int( i ); best = tr.s[1]; }
    }
    for ( const auto& tr : pool ) {
      const size_t i = tr.i;
      o_run = run; o_sub = sub; o_evt = evt; o_entry = ie; o_pfp = int( i ); o_topo = topo; o_nupdg = nu_pdg; o_ccnc = ccnc;
      o_w = safe_weight( double( w_tune ) * double( w_ppfx ) );
      o_pdg = int( at( bt_pdg, i, 0.f ) ); o_pur = at( bt_pur, i, 0.f ); o_comp = at( bt_comp, i, 0.f );
      const double px = at( bt_px, i, 0.f ), py = at( bt_py, i, 0.f ), pz = at( bt_pz, i, 0.f );
      o_truep = float( std::sqrt( px * px + py * py + pz * pz ) );
      const int apdg = std::abs( o_pdg );
      o_label = apdg == 13 ? 0 : apdg == 211 ? 1 : o_pdg == 2212 ? 2 : 3;
      o_charge = ( apdg == 13 || apdg == 211 ) ? ( o_pdg > 0 ? ( apdg == 13 ? -1 : 1 ) : ( apdg == 13 ? 1 : -1 ) ) : ( o_pdg == 2212 ? 1 : 0 );
      o_mucand = ( int( i ) == mu_idx ); o_nprimtrk = n_prim_trk; o_nprimshr = n_prim_shr;
      const float mcs_v = at( mcs, i ), rng_v = at( range_mu, i );
      const float vals[] = {
        at( len, i ), at( ts, i ), tr.dist, float( contained( at( ex, i ), at( ey, i ), at( ez, i ) ) ),
        float( contained( at( sx, i ), at( sy, i ), at( sz, i ) ) ), at( llr, i ), at( llr_u, i ), at( llr_v, i ), at( llr_y, i ),
        at( bp, i ), at( bmu, i ), at( bmip, i ), at( bpi, i ), at( bp_u, i ), at( bp_v, i ), at( bmu_u, i ), at( bmu_v, i ),
        at( bpi_u, i ), at( bpi_v, i ), at( bmip_u, i ), at( bmip_v, i ), at( fwd_p, i ), at( fwd_mu, i ), at( fwd_pi, i ),
        at( pida, i ), at( pida_u, i ), at( pida_v, i ), at( chipr, i ), at( chimu, i ), at( chipi, i ), at( chika, i ),
        at( chipr_u, i ), at( chipr_v, i ), at( chimu_u, i ), at( chimu_v, i ), at( chipi_u, i ), at( chipi_v, i ),
        at( trunk_y, i ), at( trunk_u, i ), at( trunk_v, i ), at( trunkrr_y, i ), at( trunkrr_u, i ), at( trunkrr_v, i ),
        at( calo_y, i ), at( calo_u, i ), at( calo_v, i ), at( defl_mean, i ), at( defl_std, i ), at( defl_sep, i ),
        mcs_v, rng_v, ( rng_v > 0.f && mcs_v > 0.f ) ? mcs_v / rng_v : -9999.f, at( e_p, i ), at( endsp, i ),
        at( nh_u, i ), at( nh_v, i ), at( nh_y, i ), at( ph_u, i ), at( ph_v, i ), at( ph_y, i ),
        at( ntrkd, i ), at( nshrd, i ), at( ndesc, i ),
        tr.s[0], tr.s[1], tr.s[2], tr.s[3], tr.s[4], tr.s[5] };
      static_assert( sizeof( vals ) / sizeof( float ) == 69, "feature list and names out of step" );
      for ( size_t k = 0; k < f.size(); ++k ) f[k] = vals[k];
      o.Fill(); ++n_trk;
    }
  }
  fout.cd(); o.Write();
  printf( "DUMP_DONE %s: %ld events read, %ld pass the preselection, %ld tracks\n", in_path, n, n_ev_pass, n_trk );
}
