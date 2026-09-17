// bkgfit_templates.C -- templates for the joint FHC+RHC background fit (BACKGROUND_FIT_PLAN.md).
//
// Everything the fit needs, per RUN PERIOD, in one ROOT file:
//   * neutrino-MC central value, sum of w^2, and every systematic universe, split by
//     background CLASS x (nu | nubar);
//   * detector-variation central value and the eight knobs (Run-4 samples, per beam mode), by class;
//   * per-run beam-off (gate-scaled), dirt (exposure-scaled), and beam-on data (control-region skim).
// Analysis bins (flattened, see the "bins" TNamed):
//   0      SR     selected signal region, total                       (MC only -- blind)
//   1      SIGGEN generated signal in the fiducial volume, total       (MC only)
//   2..22  CC0pi  reco cos(theta_mu) [7] x p_mu [3]
//   23..43 PI0    reco cos(theta_mu) [7] x p_mu [3]
//   44     MULTI  multi-pi total (overlaps PI0: diagnostic only)
//   45     COSMIC cosmic-enriched total
// Classes (truth, same counters as sb_transfer.C): 0 signal, 1 numuCC 0pi, 2 CC+NC with a pi0,
// 3 numuCC >=2pi+- no pi0, 4 everything else (OOFV, nueCC, NC without pi0, CC 1pi+- outside the
// phase space, other CC); index = 2*class + (nu_pdg<0).
// Periods: 0 FHC-R1 1 FHC-R2 2 FHC-R4 3 FHC-R5 4 RHC-R1 5 RHC-R2 6 RHC-R3 7 RHC-R4. Run 2 has MC only
// (no beam-on file exists) and no beam-off; it is used for full-exposure sensitivity studies.
// Weights: framework safe-weight rule; universes follow UniverseMaker::apply_cv_correction_weights
// (UBGenie * ppfx * norm; ppfx_all * tune * norm; reint and SCC * tune * ppfx * norm).
//   root -l -b -q macros/bkgfit_templates.C+    -> /data/uboone/processed/bkgfit/templates.root
#include "sb_guard.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TNamed.h"
#include "TEntryList.h"
#include "TLeaf.h"
#include "TDirectory.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <map>
#include <cmath>

namespace bkt {
  const int NB = 46, NC = 10, NP = 8;
  const std::vector<double> CTH = { -1, 0.0, 0.45, 0.65, 0.8, 0.9, 0.95, 1.0 };
  const std::vector<double> PMU = { 0.15, 0.35, 0.75 };   // last bin open
  inline int kin_bin( double c, double p ) {
    int ic = -1; for ( int i = 0; i + 1 < (int)CTH.size(); ++i ) if ( c >= CTH[i] && c < CTH[i+1] ) ic = i;
    if ( c == 1.0 ) ic = (int)CTH.size() - 2;
    int ip = -1; if ( p >= PMU[0] ) { ip = 0; if ( p >= PMU[1] ) ip = 1; if ( p >= PMU[2] ) ip = 2; }
    if ( ic < 0 || ip < 0 ) return -1;
    return ic * 3 + ip;
  }
  struct Syst { std::string branch; std::string name; int rule; };  // rule 0 UBGenie, 1 ppfx, 2 tune*ppfx
  const std::vector<Syst> SYST = {
    { "weight_ppfx_all", "flux", 1 }, { "weight_All_UBGenie", "xsec_multi", 0 }, { "weight_reint_all", "reint", 2 },
    { "weight_AxFFCCQEshape_UBGenie", "xsec_AxFFCCQEshape", 0 }, { "weight_DecayAngMEC_UBGenie", "xsec_DecayAngMEC", 0 },
    { "weight_NormCCCOH_UBGenie", "xsec_NormCCCOH", 0 }, { "weight_NormNCCOH_UBGenie", "xsec_NormNCCOH", 0 },
    { "weight_RPA_CCQE_UBGenie", "xsec_RPA_CCQE", 0 }, { "weight_ThetaDelta2NRad_UBGenie", "xsec_ThetaDelta2NRad", 0 },
    { "weight_Theta_Delta2Npi_UBGenie", "xsec_Theta_Delta2Npi", 0 }, { "weight_VecFFCCQEshape_UBGenie", "xsec_VecFFCCQEshape", 0 },
    { "weight_XSecShape_CCMEC_UBGenie", "xsec_XSecShape_CCMEC", 0 },
    { "weight_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fa3_SCC", 2 }, { "weight_xsr_scc_Fv3_SCC", "xsec_xsr_scc_Fv3_SCC", 2 } };
  inline double safe( double w ) { return ( std::isfinite( w ) && w >= 0 && w <= 30 ) ? w : 1.0; }
}

// Fill CV / w2 (and optionally universes) for one file into period p. mode: 0 nu-MC, 1 unweighted (EXT/data),
// 2 weighted without classes (dirt)
static void fill_file( const std::string& file, double scale, int kind, bool universes,
                       std::vector<double>& cv, std::vector<double>& w2,
                       std::vector<std::vector<std::vector<double>>>* U ) {
  using namespace bkt;
  TFile f( file.c_str() ); TTree* t = (TTree*)f.Get( "stv_tree" );
  if ( !t ) { printf( "  MISSING %s\n", file.c_str() ); return; }
  const char* S = "CC1mu1piXp";
  TString sel = Form( "%s_sb_cc0pi||%s_sb_pi0||%s_sb_multipi||%s_sb_cosmic", S, S, S, S );
  sel += Form( "||%s_Selected", S );            // beam-off and dirt in the signal region belong to its background
  if ( kind != 1 ) sel += Form( "||%s_MC_Signal", S );
  t->Draw( ">>el", sel, "entrylist" ); TEntryList* el = (TEntryList*)gDirectory->Get( "el" );
  t->SetBranchStatus( "*", 0 );
  auto on = [&]( const char* b ) { if ( t->GetBranch( b ) ) t->SetBranchStatus( b, 1 ); };
  const char* flags[6] = { "Selected", "MC_Signal", "sb_cc0pi", "sb_pi0", "sb_multipi", "sb_cosmic" };
  for ( auto fl : flags ) on( Form( "%s_%s", S, fl ) );
  on( Form( "%s_candidate_muon_costh_reco", S ) ); on( Form( "%s_candidate_muon_mom_reco", S ) );
  bool b[6] = { false }; for ( int i = 0; i < 6; ++i ) if ( t->GetBranch( Form( "%s_%s", S, flags[i] ) ) ) t->SetBranchAddress( Form( "%s_%s", S, flags[i] ), &b[i] );
  float cth = 0, pmu = 0;
  TLeaf* lc = t->GetLeaf( Form( "%s_candidate_muon_costh_reco", S ) ), *lp = t->GetLeaf( Form( "%s_candidate_muon_mom_reco", S ) );
  float tune = 1, ppfx = 1, norm = 1; int cat = -1, npi = 0, npi0 = 0, nupdg = 14;
  TLeaf* lcat = nullptr, *lnpi = nullptr, *lnpi0 = nullptr, *lpdg = nullptr;
  if ( kind != 1 ) {
    for ( auto w : { "tuned_cv_weight", "ppfx_cv_weight", "normalisation_weight" } ) on( w );
    t->SetBranchAddress( "tuned_cv_weight", &tune ); t->SetBranchAddress( "ppfx_cv_weight", &ppfx ); t->SetBranchAddress( "normalisation_weight", &norm );
    on( Form( "%s_EventCategory", S ) ); on( Form( "%s_mc_n_threshold_pionpm", S ) ); on( Form( "%s_mc_n_threshold_pion0", S ) ); on( "mc_nu_pdg" );
    lcat = t->GetLeaf( Form( "%s_EventCategory", S ) ); lnpi = t->GetLeaf( Form( "%s_mc_n_threshold_pionpm", S ) );
    lnpi0 = t->GetLeaf( Form( "%s_mc_n_threshold_pion0", S ) ); lpdg = t->GetLeaf( "mc_nu_pdg" );
  }
  std::vector<std::vector<double>*> wv( SYST.size(), nullptr );
  if ( universes ) for ( size_t s = 0; s < SYST.size(); ++s ) if ( t->GetBranch( SYST[s].branch.c_str() ) ) { on( SYST[s].branch.c_str() ); t->SetBranchAddress( SYST[s].branch.c_str(), &wv[s] ); }

  std::vector<int> bins; bins.reserve( 4 );
  for ( Long64_t e = 0; e < el->GetN(); ++e ) {
    t->GetEntry( el->GetEntry( e ) );
    int cls = 0;
    if ( kind == 0 ) {
      cat = (int)lcat->GetValue(); npi = (int)lnpi->GetValue(); npi0 = (int)lnpi0->GetValue(); nupdg = (int)lpdg->GetValue();
      int c;
      if ( cat == 0 ) c = 0;
      else if ( cat == 3 && npi == 0 && npi0 == 0 ) c = 1;
      else if ( ( cat == 3 || cat == 5 ) && npi0 >= 1 ) c = 2;
      else if ( cat == 3 && npi >= 2 && npi0 == 0 ) c = 3;
      else c = 4;
      cls = 2 * c + ( nupdg < 0 ? 1 : 0 );
    }
    double w0 = ( kind == 1 ) ? 1.0 : safe( tune * ppfx * norm ), w = w0 * scale;
    bins.clear();
    if ( b[0] ) bins.push_back( 0 );
    if ( b[1] && kind == 0 && cat == 0 ) bins.push_back( 1 );
    double c = lc ? lc->GetValue() : -9, p = lp ? lp->GetValue() : -9;
    int kb = kin_bin( c, p );
    if ( b[2] && kb >= 0 ) bins.push_back( 2 + kb );
    if ( b[3] && kb >= 0 ) bins.push_back( 23 + kb );
    if ( b[4] ) bins.push_back( 44 );
    if ( b[5] ) bins.push_back( 45 );
    if ( bins.empty() ) continue;
    for ( int bb : bins ) { cv[cls * NB + bb] += w; w2[cls * NB + bb] += w * w; }
    if ( universes && U ) for ( size_t s = 0; s < SYST.size(); ++s ) {
      if ( !wv[s] ) continue;
      auto& us = ( *U )[s]; size_t nu = wv[s]->size();
      if ( us.empty() ) us.assign( nu, std::vector<double>( NC * NB, 0. ) );
      for ( size_t u = 0; u < nu && u < us.size(); ++u ) {
        double wu = ( *wv[s] )[u];
        double ww = SYST[s].rule == 0 ? wu * ppfx * norm : SYST[s].rule == 1 ? wu * tune * norm : wu * tune * ppfx * norm;
        ww = safe( ww ) * scale;
        for ( int bb : bins ) us[u][cls * NB + bb] += ww;
      }
    }
  }
  delete el;
}

void bkgfit_templates( const char* outdir = "/data/uboone/processed/bkgfit" ) {
  using namespace bkt;
  // sb_pi0/ (2026-09-06): the sideband-instrumented reprocess with the CURRENT selection. The older sb/ (2026-08-11)
  // predates the beam-frame correction and stores cos(theta_mu) about detector z, while the beam-on skim stores the
  // beam-frame angle -- comparing the two produced the control-region muon-angle "shape failure".
  const std::string P = "/data/uboone/processed/sb_pi0/", DV = "/data/uboone/processed/sb_pi0/", BO = "/data/uboone/processed/beamon_skim/",
                    E = "/data/uboone/processed/sb_pi0/xsec-ana-";
  const char* pname[NP] = { "FHC_R1", "FHC_R2", "FHC_R4", "FHC_R5", "RHC_R1", "RHC_R2", "RHC_R3", "RHC_R4" };
  const double POT[NP] = { 3.283, 1.268, 2.075, 2.231, 0.6053, 2.591, 5.003, 2.883 };
  std::vector<std::vector<std::pair<std::string, double>>> mc( NP ), ext( NP );
  std::vector<std::vector<std::string>> data( NP );
  auto fmc = [&]( const char* n ) { return P + "xsec-ana-" + n + ".root"; };
  mc[0] = { { fmc( "Run1_fhc_new_numi_flux_fhc_pandora_ntuple" ), 0.14101 } };
  mc[1] = { { fmc( "Run2_fhc_new_numi_flux_fhc_pandora_ntuple" ), 0.05085 } };
  mc[2] = { { fmc( "Run4_fhc_new_numi_flux_fhc_pandora_ntuple" ), 0.07323 } };
  mc[3] = { { fmc( "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc" ), 0.11560 } };
  mc[4] = { { fmc( "Run1_rhc_new_numi_flux_rhc_pandora_ntuple" ), 0.06728 } };
  mc[5] = { { fmc( "Run2_rhc_new_numi_flux_rhc_pandora_ntuple" ), 0.04478 } };
  for ( auto s : { "aa", "ab", "ac", "ad", "ae" } ) mc[6].push_back( { fmc( ( std::string( "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_" ) + s ).c_str() ), 0.09066 } );
  for ( auto s : { "Run4a", "Run4b", "Run4c" } ) mc[7].push_back( { fmc( ( std::string( s ) + "_rhc_new_numi_flux_rhc_pandora_ntuple" ).c_str() ), 0.08847 } );
  // beam-off: exactly the gate ratios of macros/sb_protocol_data.C
  const double G1 = 4582248.27, G3 = 32649128.65, G4 = 34831148.625, G5 = 19256341.475, OCCX = 0.98;
  const char* R4[4] = { "numi_pelee_ntuple_beam_off_run4a_rhc_ana.root", "numi_pelee_ntuple_beam_off_run4b_rhc_ana.root",
                        "numi_pelee_ntuple_beam_off_run4c_fhc_ana.root", "numi_pelee_ntuple_beam_off_run4d_fhc_ana.root" };
  ext[0] = { { E + "neutrinoselection_filt_run1_beamoff.root", OCCX * 9846635. / G1 } };
  for ( auto f : R4 ) ext[2].push_back( { E + f, OCCX * 4131149. / G4 } );
  ext[3] = { { E + "numi_pelee_ntuple_beam_off_run5_fhc_ana.root", OCCX * 5154196. / G5 } };
  ext[4] = { { E + "neutrinoselection_filt_run1_beamoff.root", OCCX * 1458253. / G1 } };
  ext[6] = { { E + "neutrinoselection_filt_run3b_beamoff.root", OCCX * 10349610. / G3 } };
  for ( auto f : R4 ) ext[7].push_back( { E + f, OCCX * 6304167. / G4 } );
  data[0] = { BO + "xsec-ana-beamon_fhc_run1.root" }; data[2] = { BO + "xsec-ana-beamon_fhc_run4c.root", BO + "xsec-ana-beamon_fhc_run4d.root" };
  data[3] = { BO + "xsec-ana-beamon_fhc_run5.root" }; data[4] = { BO + "xsec-ana-beamon_rhc_run1.root" };
  data[6] = { BO + "xsec-ana-beamon_rhc_run3b.root" }; data[7] = { BO + "xsec-ana-beamon_rhc_run4a.root", BO + "xsec-ana-beamon_rhc_run4b.root" };
  const std::string DIRT = P + "xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root";

  gSystem->mkdir( outdir, true );
  TFile fo( Form( "%s/templates.root", outdir ), "RECREATE" );
  auto write1 = [&]( const std::string& name, const std::vector<double>& v, int n ) {
    TH1D h( name.c_str(), "", n, 0, n ); for ( int i = 0; i < n; ++i ) h.SetBinContent( i + 1, v[i] ); fo.cd(); h.Write(); };
  std::map<std::string, int> nuniv;
  for ( int p = 0; p < NP; ++p ) {
    printf( "== period %s\n", pname[p] ); fflush( stdout );
    std::vector<double> cv( NC * NB, 0. ), w2( NC * NB, 0. );
    std::vector<std::vector<std::vector<double>>> U( SYST.size() );
    for ( auto& x : mc[p] ) { fill_file( x.first, x.second, 0, true, cv, w2, &U ); printf( "   mc %s\n", x.first.c_str() ); fflush( stdout ); }
    write1( Form( "cv_%s", pname[p] ), cv, NC * NB ); write1( Form( "w2_%s", pname[p] ), w2, NC * NB );
    for ( size_t s = 0; s < SYST.size(); ++s ) {
      int nu = U[s].size(); nuniv[SYST[s].name] = nu; if ( !nu ) continue;
      TH2F h( Form( "u_%s_%s", SYST[s].name.c_str(), pname[p] ), "", nu, 0, nu, NC * NB, 0, NC * NB );
      for ( int u = 0; u < nu; ++u ) for ( int i = 0; i < NC * NB; ++i ) h.SetBinContent( u + 1, i + 1, U[s][u][i] );
      fo.cd(); h.Write();
    }
    std::vector<double> e( NC * NB, 0. ), e2( NC * NB, 0. );
    for ( auto& x : ext[p] ) fill_file( x.first, x.second, 1, false, e, e2, nullptr );
    write1( Form( "ext_%s", pname[p] ), e, NB ); write1( Form( "ext2_%s", pname[p] ), e2, NB );
    std::vector<double> d( NC * NB, 0. ), d2( NC * NB, 0. );
    // dirt: the single run-1 snapshot scaled per mode exactly as sb_protocol_data.C, split by period exposure
    double sdirt = p < 4 ? 0.092402 * 0.65 * POT[p] / 8.857 : 0.071666 * 0.65 * POT[p] / 11.082;
    fill_file( DIRT, sdirt, 2, false, d, d2, nullptr );
    std::vector<double> dsum( NB, 0. ), dsum2( NB, 0. );
    for ( int c = 0; c < NC; ++c ) for ( int i = 0; i < NB; ++i ) { dsum[i] += d[c * NB + i]; dsum2[i] += d2[c * NB + i]; }
    dsum[1] = 0; dsum2[1] = 0;   // generated signal is neutrino MC only
    write1( Form( "dirt_%s", pname[p] ), dsum, NB ); write1( Form( "dirt2_%s", pname[p] ), dsum2, NB );
    std::vector<double> dat( NC * NB, 0. ), dat2( NC * NB, 0. );
    for ( auto& x : data[p] ) { sb_guard_data( x ); fill_file( x, 1.0, 1, false, dat, dat2, nullptr ); }
    dat[0] = dat[1] = 0;   // the skim holds no signal-region event by construction
    write1( Form( "data_%s", pname[p] ), dat, NB );
    TNamed hasdata( Form( "hasdata_%s", pname[p] ), data[p].empty() ? "0" : "1" ); fo.cd(); hasdata.Write();
  }
  // detector variations, per mode, by class (unscaled; only fractional shifts are used)
  const char* knobs[8] = { "LYdown", "LYrayl", "Recomb2", "SCE", "WMAngleXZ", "WMAngleYZ", "WMX", "WMYZ" };
  for ( auto m : { "run4fhc", "run4rhc" } ) {
    std::vector<std::string> names = { "CV" }; for ( auto k : knobs ) names.push_back( k );
    for ( auto& k : names ) {
      std::vector<double> cv( NC * NB, 0. ), w2( NC * NB, 0. );
      fill_file( DV + "xsec-ana-detvar_" + m + "_" + k + ".root", 1.0, 0, false, cv, w2, nullptr );
      write1( Form( "dv_%s_%s", m, k.c_str() ), cv, NC * NB ); write1( Form( "dv2_%s_%s", m, k.c_str() ), w2, NC * NB );
      printf( "   detvar %s %s\n", m, k.c_str() ); fflush( stdout );
    }
  }
  std::string meta = "NB=46 NC=10 NP=8; bins: 0 SR, 1 SIGGEN, 2-22 CC0pi cth[-1,0,0.45,0.65,0.8,0.9,0.95,1]xpmu[0.15,0.35,0.75,inf] (index 2+ic*3+ip), "
                     "23-43 PI0 same, 44 MULTI, 45 COSMIC; classes idx=2*c+nubar, c: 0 signal 1 CC0pi 2 pi0(CC+NC) 3 multipi 4 other; "
                     "periods FHC_R1 FHC_R2 FHC_R4 FHC_R5 RHC_R1 RHC_R2 RHC_R3 RHC_R4; POT(1e20) 3.283 1.268 2.075 2.231 0.6053 2.591 5.003 2.883";
  TNamed tn( "bins", meta.c_str() ); fo.cd(); tn.Write();
  std::string nus; for ( auto& kv : nuniv ) nus += kv.first + "=" + std::to_string( kv.second ) + " ";
  TNamed tu( "nuniv", nus.c_str() ); tu.Write();
  fo.Close();
  printf( "wrote %s/templates.root  universes: %s\n", outdir, nus.c_str() );
}
