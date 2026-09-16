// binning_screen_2d.C -- single-pass pre-screen of 1D AND 2D candidate binnings, evaluated
// against BOTH the 0.68 migration-diagonal criterion and the 0.50 alternative under study.
//
// Same quantities and conventions as macros/binning_screen.C (column-normalised diagonal per
// true bin, efficiency per true bin, selected sig+bkg per reco bin, per-run scales, safe CV
// weight, first/last edge semantics), but every candidate is filled in ONE pass over each
// file instead of four TTree::Draw calls per bin, so a long candidate list is affordable.
//
// A 2D candidate is a Y axis (outer edges) and, per Y slice, its own X edges. Analysis bins
// are flattened slice by slice (Y slice 0 X bins, then Y slice 1, ...), which is the order a
// 2D bin config must use so that the per-slice A_C block in UnfolderNuMI stays contiguous.
// For 2D candidates the off-diagonal migration of selected signal is split into
//   X-leak : reconstructed in the right Y slice but the wrong X bin
//   Y-leak : reconstructed in the wrong Y slice
// so it is visible which axis limits the binning.
//
//   root -l -b -q 'macros/binning_screen_2d.C+("fhc","incl")'   (rhc; incl | 1p)
// Output: printed table + logs/c50/screen_<mode>_<sample>.tsv (one line per true bin).
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace scr {
  struct Src { std::string file; double scale; };
  struct Cand {
    std::string sample, obs, label;
    std::string xt, xr;                          // X true / reco expressions
    std::vector<std::vector<double>> xe;         // X edges per Y slice (size 1 -> shared)
    bool xOpenTop;
    std::string yt = "", yr = "";                // empty -> 1D
    std::vector<double> ye = {};
    bool yOpenTop = false;
    int nslice() const { return yt.empty() ? 1 : (int)ye.size() - 1; }
    const std::vector<double>& xedges( int j ) const { return xe.size() == 1 ? xe[0] : xe.at(j); }
    int nbins() const { int n = 0; for ( int j = 0; j < nslice(); ++j ) n += (int)xedges(j).size() - 1; return n; }
    int offset( int j ) const { int n = 0; for ( int k = 0; k < j; ++k ) n += (int)xedges(k).size() - 1; return n; }
  };
  // bin of v in edges; the last bin has no upper edge when openTop (both truth and reco,
  // as in binning_screen.C); values below the first edge or above the last are out (-1)
  inline int find( double v, const std::vector<double>& e, bool openTop ) {
    int n = (int)e.size() - 1;
    for ( int b = 0; b < n; ++b ) {
      if ( v < e[b] ) return -1;
      if ( v < e[b+1] || ( openTop && b == n-1 ) ) return b;
    }
    return -1;
  }
  inline int flat( const Cand& c, double x, double y ) {
    int j = 0;
    if ( !c.yt.empty() ) { j = find( y, c.ye, c.yOpenTop ); if ( j < 0 ) return -1; }
    int b = find( x, c.xedges(j), c.xOpenTop ); if ( b < 0 ) return -1;
    return c.offset(j) + b;
  }
  inline int slice_of( const Cand& c, int idx ) {
    for ( int j = 0; j < c.nslice(); ++j ) if ( idx < c.offset(j) + (int)c.xedges(j).size() - 1 ) return j;
    return -1;
  }
}

void binning_screen_2d( const char* mode = "fhc", const char* sample = "incl" ) {
  using namespace scr;
  const std::string S = std::string(sample) == "1p" ? "CC1mu1pi1p" : "CC1mu1piXp";
  const std::string P = std::string("/data/uboone/processed/") + ( std::string(sample) == "1p" ? "w/" : "" );
  std::vector<Src> mc;
  if ( std::string(mode) == "fhc" ) {
    const char* rn[4] = {"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
                         "Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
    double sc[4] = {0.14101,0.05085,0.07323,0.11560};
    for ( int i = 0; i < 4; ++i ) mc.push_back( { P+"xsec-ana-"+rn[i]+".root", sc[i] } );
  } else {
    const char* rn[5] = {"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
    double sc[5] = {0.06728,0.04478,0.08847,0.08847,0.08847};
    for ( int i = 0; i < 5; ++i ) mc.push_back( { P+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root", sc[i] } );
    for ( auto s : {"aa","ab","ac","ad","ae"} )
      mc.push_back( { P+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root", 0.09066 } );
  }
  const std::string CVW = "(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  auto X = [&]( const char* v ){ return S + "_" + v; };
  const std::string PMT = X("candidate_muon_mom_true"),  PMR = X("candidate_muon_mom_reco");
  const std::string CMT = X("candidate_muon_costh_true"), CMR = X("candidate_muon_costh_reco");
  const std::string CPT = X("candidate_pion_costh_true"), CPR = X("candidate_pion_costh_reco");
  const std::string TMPT = X("true_mu_pi_opening_angle"), TMPR = X("mu_pi_opening_angle");
  const std::string PPT = X("candidate_pion_mom_true");
  const std::string PPR = "sqrt(pow(sqrt(pow(" + X("candidate_pion_mom_reco") + ",2)+0.011164)-0.10566+0.13957,2)-0.019480)";
  const std::string TH  = "TMath::ACos(" + CMT + ")", THR = "TMath::ACos(" + CMR + ")";
  const std::string THP = "TMath::ACos(" + X("proton_costh_true") + ")", THPR = "TMath::ACos(" + X("proton_costh_reco") + ")";
  const double PI = 3.1416;

  std::vector<Cand> C = {
    // ------------------------------------------------------------------ inclusive, 1D
    {"incl","pmu","released 7",PMT,PMR,{{0.15,0.35,0.55,0.75,0.95,1.25,1.75,3.0}},true},
    {"incl","pmu","cand 8",PMT,PMR,{{0.15,0.35,0.55,0.75,0.95,1.25,1.5,1.75,3.0}},true},
    {"incl","pmu","cand 9",PMT,PMR,{{0.15,0.30,0.45,0.60,0.75,0.95,1.2,1.5,2.0,3.0}},true},
    {"incl","costhmu","released 5",CMT,CMR,{{-1,0.45,0.65,0.8,0.9,1}},false},
    {"incl","costhmu","cand 6",CMT,CMR,{{-1,0.3,0.55,0.7,0.8,0.9,1}},false},
    {"incl","costhmu","cand 7",CMT,CMR,{{-1,0.2,0.45,0.6,0.7,0.8,0.9,1}},false},
    {"incl","costhmu","cand 8",CMT,CMR,{{-1,0.2,0.4,0.55,0.65,0.75,0.85,0.93,1}},false},
    {"incl","thetamu","released 5",TH,THR,{{0,0.43,0.62,0.84,1.17,3.15}},false},
    {"incl","thetamu","cand 7",TH,THR,{{0,0.35,0.5,0.65,0.8,1.0,1.3,3.15}},false},
    {"incl","thetamu","cand 8",TH,THR,{{0,0.3,0.42,0.55,0.68,0.84,1.05,1.4,3.15}},false},
    {"incl","costhpi","old released 4",CPT,CPR,{{-1,-0.1,0.35,0.75,1}},false},
    {"incl","costhpi","5 bwd split",CPT,CPR,{{-1,-0.3,0.2,0.5,0.75,1}},false},
    {"incl","costhpi","released 5 (fwd split)",CPT,CPR,{{-1,-0.1,0.35,0.6,0.8,1}},false},
    {"incl","costhpi","cand 6",CPT,CPR,{{-1,-0.4,0,0.3,0.55,0.75,1}},false},
    {"incl","costhpi","cand 6 fwd",CPT,CPR,{{-1,-0.1,0.25,0.5,0.7,0.85,1}},false},
    {"incl","costhpi","cand 7",CPT,CPR,{{-1,-0.4,0,0.3,0.5,0.7,0.85,1}},false},
    {"incl","thmupi","released 5",TMPT,TMPR,{{0,0.6,0.85,1.3,1.85,2.6}},false},
    {"incl","thmupi","cand 6",TMPT,TMPR,{{0,0.55,0.75,0.95,1.3,1.85,2.6}},false},
    {"incl","thmupi","cand 7",TMPT,TMPR,{{0,0.5,0.7,0.9,1.15,1.45,1.9,2.6}},false},
    {"incl","thmupi","cand 8",TMPT,TMPR,{{0,0.45,0.6,0.75,0.9,1.1,1.35,1.8,2.6}},false},
    {"incl","ppi","released 2",PPT,PPR,{{0.175,0.205,1.0}},true},
    {"incl","ppi","cand 3",PPT,PPR,{{0.175,0.205,0.30,1.0}},true},
    {"incl","ppi","cand 3b",PPT,PPR,{{0.175,0.205,0.26,1.0}},true},
    {"incl","ppi","cand 4",PPT,PPR,{{0.175,0.205,0.26,0.34,1.0}},true},
    // ------------------------------------------------------------------ inclusive, 2D (X | Y)
    {"incl","costhpi|costhmu","3x2",CPT,CPR,{{-1,0.35,1}},false,CMT,CMR,{-1,0.65,0.85,1},false},
    {"incl","costhpi|costhmu","3x3",CPT,CPR,{{-1,-0.1,0.5,1}},false,CMT,CMR,{-1,0.65,0.85,1},false},
    {"incl","costhpi|costhmu","4x2",CPT,CPR,{{-1,0.35,1}},false,CMT,CMR,{-1,0.45,0.7,0.88,1},false},
    {"incl","costhpi|costhmu","2x3",CPT,CPR,{{-1,-0.1,0.5,1}},false,CMT,CMR,{-1,0.75,1},false},
    {"incl","costhmu|pmu","2x2",CMT,CMR,{{-1,0.8,1}},false,PMT,PMR,{0.15,0.55,3.0},true},
    {"incl","costhmu|pmu","2x3",CMT,CMR,{{-1,0.7,0.9,1}},false,PMT,PMR,{0.15,0.55,3.0},true},
    {"incl","costhmu|pmu","3x2",CMT,CMR,{{-1,0.8,1}},false,PMT,PMR,{0.15,0.45,0.85,3.0},true},
    {"incl","costhmu|pmu","3x3",CMT,CMR,{{-1,0.7,0.9,1}},false,PMT,PMR,{0.15,0.45,0.85,3.0},true},
    {"incl","costhpi|thmupi","3x2",CPT,CPR,{{-1,0.35,1}},false,TMPT,TMPR,{0,0.85,1.5,2.6},false},
    {"incl","costhpi|thmupi","3x3",CPT,CPR,{{-1,-0.1,0.5,1}},false,TMPT,TMPR,{0,0.85,1.5,2.6},false},
    {"incl","costhpi|thmupi","2x3",CPT,CPR,{{-1,-0.1,0.5,1}},false,TMPT,TMPR,{0,1.1,2.6},false},
    {"incl","costhpi|ppi","2x2",CPT,CPR,{{-1,0.35,1}},false,PPT,PPR,{0.175,0.205,1.0},true},
    {"incl","costhpi|ppi","2x3",CPT,CPR,{{-1,-0.1,0.5,1}},false,PPT,PPR,{0.175,0.205,1.0},true},
    // ------------------------------------------------------------------ proton-tagged, 1D
    // released two-bin TKI / W, then the maximin-scan edges of macros/tki_binning_scan.C
    // (FHC Run1+2+4, >=400 events/bin) at the largest bin count passing 0.68 and 0.50
    {"1p","deltaPt","released 2",X("deltaPt_true"),X("deltaPt_reco"),{{0,0.3,2.5}},true},
    {"1p","deltaPt","c50 3",X("deltaPt_true"),X("deltaPt_reco"),{{0,0.225,0.45,2.5}},true},
    {"1p","dalphat","released 2",X("deltaAlphaT_true"),X("deltaAlphaT_reco"),{{0,120,180}},false},
    {"1p","dalphat","c50 3",X("deltaAlphaT_true"),X("deltaAlphaT_reco"),{{0,85,150,180}},false},
    {"1p","dphit","released 2",X("deltaPhiT_true"),X("deltaPhiT_reco"),{{0,50,180}},false},
    {"1p","dphit","c50 4",X("deltaPhiT_true"),X("deltaPhiT_reco"),{{0,15,50,90,180}},false},
    {"1p","pn","released 2",X("pn_true"),X("pn_reco"),{{0,0.375,2.0}},true},
    {"1p","pn","c50 4",X("pn_true"),X("pn_reco"),{{0,0.2,0.375,0.625,2.0}},true},
    {"1p","Wpipr","released 2",X("W_pipr_true"),X("W_pipr_reco"),{{1.08,1.15,2.90}},true},
    {"1p","Wpipr","c50 2",X("W_pipr_true"),X("W_pipr_reco"),{{1.08,1.19,2.90}},true},
    {"1p","Whad","c50 2",X("W_had_true"),X("W_had_reco"),{{0,1.18,2.74}},true},
    {"1p","costhp","c68 4",X("proton_costh_true"),X("proton_costh_reco"),{{-1,0.275,0.575,0.825,1}},false},
    {"1p","costhp","c50 6",X("proton_costh_true"),X("proton_costh_reco"),{{-1,0.275,0.575,0.7,0.8,0.9,1}},false},
    {"1p","thetap","c68 5",THP,THPR,{{0,0.377,0.691,0.942,1.225,PI}},false},
    {"1p","thetap","c50 7",THP,THPR,{{0,0.283,0.471,0.659,0.816,1.005,1.287,PI}},false},
    {"1p","thpipr","c68 3","pa.thpipr_true","pa.thpipr_reco",{{0,1.193,2.010,PI}},false},
    {"1p","thpipr","c50 7","pa.thpipr_true","pa.thpipr_reco",{{0,0.628,0.911,1.162,1.413,1.664,1.978,PI}},false},
    // ------------------------------------------------------------------ proton-tagged, 2D (X | Y)
    {"1p","thetap|deltaPt","2x2",THP,THPR,{{0,0.8,PI}},false,X("deltaPt_true"),X("deltaPt_reco"),{0,0.3,2.5},true},
    {"1p","thetap|deltaPt","3x2",THP,THPR,{{0,0.6,1.1,PI}},false,X("deltaPt_true"),X("deltaPt_reco"),{0,0.3,2.5},true},
    {"1p","thpipr|deltaPt","2x2","pa.thpipr_true","pa.thpipr_reco",{{0,1.2,PI}},false,X("deltaPt_true"),X("deltaPt_reco"),{0,0.3,2.5},true},
    {"1p","thpipr|deltaPt","3x2","pa.thpipr_true","pa.thpipr_reco",{{0,1.0,1.7,PI}},false,X("deltaPt_true"),X("deltaPt_reco"),{0,0.3,2.5},true},
    {"1p","thetap|Wpipr","2x2",THP,THPR,{{0,0.8,PI}},false,X("W_pipr_true"),X("W_pipr_reco"),{1.08,1.19,2.90},true},
    {"1p","thpipr|Wpipr","2x2","pa.thpipr_true","pa.thpipr_reco",{{0,1.2,PI}},false,X("W_pipr_true"),X("W_pipr_reco"),{1.08,1.19,2.90},true},
    {"1p","costhp|deltaPt","2x2",X("proton_costh_true"),X("proton_costh_reco"),{{-1,0.575,1}},false,X("deltaPt_true"),X("deltaPt_reco"),{0,0.3,2.5},true},
  };
  // Optional override: candidates appended from a text file (one per line), so scan-derived
  // edges can be screened without editing the macro:
  //   sample|obs|label|xt|xr|x edges (space separated)[;slice2 edges...]|xOpenTop|yt|yr|y edges|yOpenTop
  // SCREEN_EXTRA_ONLY=1 drops the built-in list and screens only the file's candidates
  if ( gSystem->Getenv( "SCREEN_EXTRA_ONLY" ) ) C.clear();
  if ( const char* extra = gSystem->Getenv( "SCREEN_EXTRA" ) ) {
    std::ifstream in( extra ); std::string line;
    while ( std::getline( in, line ) ) {
      if ( line.empty() || line[0] == '#' ) continue;
      std::vector<std::string> f; size_t p0 = 0, p;
      while ( ( p = line.find( '|', p0 ) ) != std::string::npos ) { f.push_back( line.substr( p0, p-p0 ) ); p0 = p+1; }
      f.push_back( line.substr( p0 ) );
      if ( f.size() != 7 && f.size() != 11 ) { printf( "SCREEN_EXTRA: bad line %s\n", line.c_str() ); continue; }
      auto nums = []( const std::string& s ){ std::vector<double> v; std::istringstream is( s ); double x; while ( is >> x ) v.push_back( x ); return v; };
      Cand c; c.sample = f[0]; c.obs = f[1]; c.label = f[2]; c.xt = f[3]; c.xr = f[4];
      size_t q0 = 0, q; while ( ( q = f[5].find( ';', q0 ) ) != std::string::npos ) { c.xe.push_back( nums( f[5].substr( q0, q-q0 ) ) ); q0 = q+1; }
      c.xe.push_back( nums( f[5].substr( q0 ) ) ); c.xOpenTop = f[6] == "1";
      if ( f.size() == 11 ) { c.yt = f[7]; c.yr = f[8]; c.ye = nums( f[9] ); c.yOpenTop = f[10] == "1"; }
      C.push_back( c );
    }
  }
  std::vector<Cand> use; for ( auto& c : C ) if ( c.sample == sample ) use.push_back( c );
  const bool needFriend = std::string(sample) == "1p";

  struct Acc { std::vector<double> gen, sel, diag, reco, xleak, yleak; };
  std::vector<Acc> A( use.size() );
  for ( size_t k = 0; k < use.size(); ++k ) { int n = use[k].nbins();
    A[k] = { std::vector<double>(n,0.), std::vector<double>(n,0.), std::vector<double>(n,0.),
             std::vector<double>(n,0.), std::vector<double>(n,0.), std::vector<double>(n,0.) }; }

  for ( auto& s : mc ) {
    TFile f( s.file.c_str() ); TTree* t = (TTree*)f.Get( "stv_tree" );
    if ( !t ) { printf( "  missing %s\n", s.file.c_str() ); continue; }
    if ( needFriend ) {
      TString b = gSystem->BaseName( s.file.c_str() ); b.ReplaceAll( ".root", ".pa.root" );
      TString ff = TString( gSystem->DirName( s.file.c_str() ) ) + "/friends_pa/" + b;
      if ( gSystem->AccessPathName( ff ) ) { printf( "  FRIEND MISSING %s -- theta_pipr candidates will be empty\n", ff.Data() ); }
      else t->AddFriend( "pa", ff );
    }
    std::map<std::string, std::unique_ptr<TTreeFormula>> F;
    auto form = [&]( const std::string& e ) -> TTreeFormula* {
      auto it = F.find( e ); if ( it != F.end() ) return it->second.get();
      auto* tf = new TTreeFormula( Form( "f%zu", F.size() ), e.c_str(), t ); tf->SetQuickLoad( false );
      F[e].reset( tf ); return tf; };
    TTreeFormula* fw = form( CVW ), *fsel = form( S+"_Selected" ), *fsig = form( S+"_MC_Signal" );
    for ( auto& c : use ) { form( c.xt ); form( c.xr ); if ( !c.yt.empty() ) { form( c.yt ); form( c.yr ); } }
    std::map<TTreeFormula*, double> val;
    long long N = t->GetEntries();
    for ( long long i = 0; i < N; ++i ) {
      t->LoadTree( i );
      fsel->GetNdata(); fsig->GetNdata();
      bool sel = fsel->EvalInstance() != 0., sig = fsig->EvalInstance() != 0.;
      if ( !sel && !sig ) continue;
      fw->GetNdata(); double w = fw->EvalInstance() * s.scale;
      val.clear();
      auto get = [&]( const std::string& e ){ TTreeFormula* tf = F[e].get(); auto it = val.find( tf );
        if ( it != val.end() ) return it->second; tf->GetNdata(); double v = tf->EvalInstance(); val[tf] = v; return v; };
      for ( size_t k = 0; k < use.size(); ++k ) {
        const Cand& c = use[k]; Acc& a = A[k]; bool twoD = !c.yt.empty();
        int tb = sig ? flat( c, get( c.xt ), twoD ? get( c.yt ) : 0. ) : -1;
        int rb = sel ? flat( c, get( c.xr ), twoD ? get( c.yr ) : 0. ) : -1;
        if ( tb >= 0 ) { a.gen[tb] += w;
          if ( sel ) { a.sel[tb] += w;
            if ( rb == tb ) a.diag[tb] += w;
            else if ( rb >= 0 ) { if ( twoD && slice_of( c, rb ) != slice_of( c, tb ) ) a.yleak[tb] += w; else a.xleak[tb] += w; } } }
        if ( rb >= 0 ) a.reco[rb] += w;
      }
    }
    printf( "  done %s\n", s.file.c_str() ); fflush( stdout );
  }

  gSystem->mkdir( "../logs/c50", true );
  FILE* tsv = fopen( Form( "../logs/c50/screen_%s_%s%s.tsv", mode, sample, gSystem->Getenv( "SCREEN_EXTRA_ONLY" ) ? "_extra" : "" ), "w" );
  fprintf( tsv, "obs\tlabel\tbin\tslice\tdiag\teff\tsel_sig\treco_all\txleak\tyleak\n" );
  printf( "\n==== %s %s: col-normalised diagonal [diag | eff | selected sig+bkg in reco bin]; 2D adds X-leak/Y-leak\n", mode, sample );
  for ( size_t k = 0; k < use.size(); ++k ) {
    const Cand& c = use[k]; Acc& a = A[k]; int nb = c.nbins();
    printf( "%-16s %-24s", c.obs.c_str(), c.label.c_str() );
    double mind = 9, minsel = 1e18, sx = 0, sy = 0, ss = 0;
    for ( int b = 0; b < nb; ++b ) {
      double d = a.sel[b] > 0 ? a.diag[b]/a.sel[b] : 0.;
      mind = std::min( mind, d ); minsel = std::min( minsel, a.sel[b] );
      sx += a.xleak[b]; sy += a.yleak[b]; ss += a.sel[b];
      printf( " [%.2f e%.0f%% n%.0f]", d, a.gen[b] > 0 ? 100*a.sel[b]/a.gen[b] : 0., a.reco[b] );
      fprintf( tsv, "%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.1f\t%.1f\t%.4f\t%.4f\n", c.obs.c_str(), c.label.c_str(), b, slice_of( c, b ), d,
               a.gen[b] > 0 ? a.sel[b]/a.gen[b] : 0., a.sel[b], a.reco[b],
               a.sel[b] > 0 ? a.xleak[b]/a.sel[b] : 0., a.sel[b] > 0 ? a.yleak[b]/a.sel[b] : 0. );
    }
    printf( "  min diag %.2f  0.68:%s 0.50:%s  min sel-sig %.0f%s", mind, mind > 0.68 ? "PASS" : "fail",
            mind > 0.50 ? "PASS" : "fail", minsel, minsel < 100 ? " (<100!)" : "" );
    if ( !c.yt.empty() && ss > 0 ) printf( "  X-leak %.1f%% Y-leak %.1f%%", 100*sx/ss, 100*sy/ss );
    printf( "\n" );
  }
  fclose( tsv );
}
