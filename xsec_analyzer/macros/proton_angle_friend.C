// proton_angle_friend.C -- proton-tagged angle observables WITHOUT reprocessing processed/w.
//
// CC1mu1pi1p writes proton_costh_{reco,true} but neither the proton candidate index nor the
// pion-proton opening angle. This macro repeats the proton pick of
// src/selections/CC1mu1pi1p.cxx::selection() on the pass-through PFP vectors (LLR < 0.05,
// track score >= 0.5, p >= 0.3 GeV/c from the proton-hypothesis KE, track end inside the
// reco FV, not the muon or pion candidate; lowest LLR wins) and writes an entry-aligned
// friend tree "pa" with
//   thpipr_reco / thpipr_true   pion-proton opening angle [rad] (frame independent)
//   costhp_reco / costhp_true   proton cos(theta) about the neutrino direction
//                               (reco: target->vertex beta convention, recomputed; true: copied
//                               from CC1mu1pi1p_proton_costh_true, since the output tree does
//                               not carry the true neutrino momentum)
//   costhp_chk                  |costhp_reco - CC1mu1pi1p_proton_costh_reco|, the check that
//                               the pick is replicated exactly (must be ~1e-6 when selected)
//   prmatch                     reco proton candidate truth: -1 none, 0 not a proton,
//                               1 a proton but not the leading true proton, 2 the leading one
// Only events that are selected or signal are filled; everything else is -9999.
//   root -l -b -q 'macros/proton_angle_friend.C+("in.root","out_friend.root")'
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TVector3.h"
#include <vector>
#include <cmath>
#include <string>
#include "XSecAnalyzer/NuMIBeamFrame.hh"

void proton_angle_friend( const char* in, const char* out ) {
  const double BOGUS = -9999., MP = 0.93827208;
  const float LLRCUT = 0.05, TSCUT = 0.5;
  const double FV[6] = { 10., 246., -101., 101., 10., 986. };  // CC1mu1piXp reco FV
  const char* S = "CC1mu1pi1p";

  TFile fin( in ); TTree* t = (TTree*)fin.Get( "stv_tree" );
  if ( !t ) { printf( "no stv_tree in %s\n", in ); return; }

  std::vector<float>* llr=nullptr,*ts=nullptr,*kep=nullptr,*ex=nullptr,*ey=nullptr,*ez=nullptr,
    *dx=nullptr,*dy=nullptr,*dz=nullptr,*mpx=nullptr,*mpy=nullptr,*mpz=nullptr,
    *bpx=nullptr,*bpy=nullptr,*bpz=nullptr;
  std::vector<int>* mpdg=nullptr,*bpdg=nullptr;
  t->SetBranchStatus( "*", 0 );
  auto on = [&]( const char* n ){ t->SetBranchStatus( n, 1 ); };
  for ( auto n : { "trk_llr_pid_score_v","trk_score_v","trk_energy_proton_v","trk_sce_end_x_v",
    "trk_sce_end_y_v","trk_sce_end_z_v","trk_dir_x_v","trk_dir_y_v","trk_dir_z_v","mc_pdg","mc_px",
    "mc_py","mc_pz","backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
    "reco_nu_vtx_sce_x","reco_nu_vtx_sce_y","reco_nu_vtx_sce_z" } ) on( n );
  for ( auto n : { "Selected","MC_Signal","CandidateMuonIndex","CandidatePionIndex","proton_costh_reco","proton_costh_true" } )
    on( Form( "%s_%s", S, n ) );
  t->SetBranchAddress( "trk_llr_pid_score_v", &llr ); t->SetBranchAddress( "trk_score_v", &ts );
  t->SetBranchAddress( "trk_energy_proton_v", &kep );
  t->SetBranchAddress( "trk_sce_end_x_v", &ex ); t->SetBranchAddress( "trk_sce_end_y_v", &ey );
  t->SetBranchAddress( "trk_sce_end_z_v", &ez );
  t->SetBranchAddress( "trk_dir_x_v", &dx ); t->SetBranchAddress( "trk_dir_y_v", &dy );
  t->SetBranchAddress( "trk_dir_z_v", &dz );
  t->SetBranchAddress( "mc_pdg", &mpdg ); t->SetBranchAddress( "mc_px", &mpx );
  t->SetBranchAddress( "mc_py", &mpy ); t->SetBranchAddress( "mc_pz", &mpz );
  t->SetBranchAddress( "backtracked_pdg", &bpdg ); t->SetBranchAddress( "backtracked_px", &bpx );
  t->SetBranchAddress( "backtracked_py", &bpy ); t->SetBranchAddress( "backtracked_pz", &bpz );
  // scalars through formulas: independent of the stored leaf types
  TTreeFormula fsel( "sel", Form( "%s_Selected", S ), t ), fsig( "sig", Form( "%s_MC_Signal", S ), t ),
    fmu( "mu", Form( "%s_CandidateMuonIndex", S ), t ), fpi( "pi", Form( "%s_CandidatePionIndex", S ), t ),
    fcos( "cos", Form( "%s_proton_costh_reco", S ), t ), fcost( "cost", Form( "%s_proton_costh_true", S ), t ),
    fvx( "vx", "reco_nu_vtx_sce_x", t ), fvy( "vy", "reco_nu_vtx_sce_y", t ), fvz( "vz", "reco_nu_vtx_sce_z", t );

  TFile fout( out, "RECREATE" );
  TTree pa( "pa", "proton-tagged angle friend of stv_tree" );
  double thr, tht, cpr, cpt, chk; int prm;
  pa.Branch( "thpipr_reco", &thr ); pa.Branch( "thpipr_true", &tht );
  pa.Branch( "costhp_reco", &cpr ); pa.Branch( "costhp_true", &cpt );
  pa.Branch( "costhp_chk", &chk ); pa.Branch( "prmatch", &prm );

  long long N = t->GetEntries(), nsel = 0, nbad = 0, nmatch[4] = {0,0,0,0};
  for ( long long i = 0; i < N; ++i ) {
    t->GetEntry( i );
    fsel.GetNdata(); fsig.GetNdata();
    bool sel = fsel.EvalInstance() != 0., sig = fsig.EvalInstance() != 0.;
    thr = tht = cpr = cpt = chk = BOGUS; prm = -1;
    if ( sel || sig ) {
      // ---- truth: leading charged pion and leading proton
      int ipi = -1, ipr = -1; double lpi = 0., lpr = 0.;
      for ( size_t k = 0; k < mpdg->size(); ++k ) {
        double p = std::sqrt( mpx->at(k)*mpx->at(k) + mpy->at(k)*mpy->at(k) + mpz->at(k)*mpz->at(k) );
        if ( std::abs( mpdg->at(k) ) == 211 && p > lpi ) { lpi = p; ipi = (int)k; }
        if ( mpdg->at(k) == 2212 && p > lpr ) { lpr = p; ipr = (int)k; }
      }
      TVector3 tpr;
      if ( ipr >= 0 ) {
        tpr.SetXYZ( mpx->at(ipr), mpy->at(ipr), mpz->at(ipr) );
        fcost.GetNdata(); cpt = fcost.EvalInstance();
        if ( ipi >= 0 ) tht = TVector3( mpx->at(ipi), mpy->at(ipi), mpz->at(ipi) ).Angle( tpr );
      }
      // ---- reco: repeat the proton pick
      int imu = (int)fmu.EvalInstance(), ipc = (int)fpi.EvalInstance(), best = -1;
      double bl = 1e9;
      for ( size_t k = 0; k < llr->size(); ++k ) {
        if ( (int)k == imu || (int)k == ipc ) continue;
        if ( ts->at(k) < TSCUT ) continue;
        if ( llr->at(k) >= LLRCUT ) continue;
        double KE = kep->at(k), pm = KE*KE + 2.*MP*KE; pm = pm > 0. ? std::sqrt( pm ) : 0.;
        if ( pm < 0.3 ) continue;
        if ( !( FV[0] < ex->at(k) && ex->at(k) < FV[1] && FV[2] < ey->at(k) && ey->at(k) < FV[3]
             && FV[4] < ez->at(k) && ez->at(k) < FV[5] ) ) continue;
        if ( llr->at(k) < bl ) { bl = llr->at(k); best = (int)k; }
      }
      if ( sel && best >= 0 && ipc >= 0 ) {
        ++nsel;
        TVector3 dpr( dx->at(best), dy->at(best), dz->at(best) ), dpi( dx->at(ipc), dy->at(ipc), dz->at(ipc) );
        thr = dpi.Angle( dpr );
        cpr = dpr.Unit().Dot( NuMIBeam::nu_dir_from_vertex( fvx.EvalInstance(), fvy.EvalInstance(), fvz.EvalInstance() ) );
        fcos.GetNdata(); chk = std::abs( cpr - fcos.EvalInstance() );
        if ( chk > 1e-4 ) ++nbad;
        if ( best < (int)bpdg->size() ) {
          if ( bpdg->at(best) != 2212 ) prm = 0;
          else {
            TVector3 bp( bpx->at(best), bpy->at(best), bpz->at(best) );
            prm = ( ipr >= 0 && bp.Angle( tpr ) < 0.1 && std::abs( bp.Mag() - lpr ) < 0.1*lpr ) ? 2 : 1;
          }
        }
        ++nmatch[prm+1];
      }
      else if ( sel ) ++nbad;  // selected but the pick found no proton: replication failure
    }
    pa.Fill();
  }
  fout.cd(); pa.Write(); fout.Close();
  printf( "[PAFRIEND] %s entries %lld selected %lld mismatched %lld | prmatch none %lld notp %lld subleading %lld leading %lld\n",
          in, N, nsel, nbad, nmatch[0], nmatch[1], nmatch[2], nmatch[3] );
}
