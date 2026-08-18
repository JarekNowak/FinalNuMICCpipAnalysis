#include <functional>
// ppi_binning_scan.C -- exhaustive scan over p_pi binning choices.
//
// For every candidate set of bin edges it evaluates the migration diagonal per true bin,
// i.e. the fraction of SELECTED SIGNAL events in a true bin whose reconstructed p_pi lands
// in the matching reco bin. That is the quantity the BNB CC1pi criterion is stated in: a
// bin is usable when its diagonal exceeds 0.68.
//
// The scan answers two questions the note needs to settle:
//   (a) can MORE than two bins be made to work, at any choice of edges?
//   (b) is 0.205 the best first boundary, or merely the one inherited from the inclusive
//       optimisation and reused for the proton-tagged selection?
//
// Method. Fill ONE fine (true, reco) histogram, then form its 2D cumulative sum so the
// counts inside any candidate bin rectangle are an O(1) lookup. Without that, scanning
// ~10^5 four-bin configurations would mean ~10^10 histogram reads.
//
// The reco axis starts at 0, not at the analysis threshold: reconstructed p_pi frequently
// falls BELOW the true threshold (that is the whole problem), and those events must be
// counted in the lowest reco bin rather than dropped.
//
//   usage: root -l -b -q 'macros/ppi_binning_scan.C("incl","fhc")'
//          selection "incl" (CC1mu1piXp) or "1p" (CC1mu1pi1p);  beam "fhc"/"rhc"/"comb"

#include <vector>
#include <algorithm>

namespace {
  constexpr double LO      = 0.175;   // analysis threshold = first true edge
  constexpr double HI      = 2.000;   // top of the scanned range; above this -> overflow
  constexpr double STEP    = 0.005;   // 5 MeV/c grid for the fine histogram
  constexpr double MINW    = 0.030;   // narrowest bin with precedent (BNB CC1pi)
  constexpr double CRIT    = 0.68;    // a bin is usable above this diagonal
  constexpr int    MINEVT  = 100;     // BNB also requires >=100 predicted events per bin
}

void ppi_binning_scan( const char* sel = "incl", const char* beam = "fhc" ) {

  const bool onep = ( std::string(sel) == "1p" );
  const char* S = onep ? "CC1mu1pi1p" : "CC1mu1piXp";

  std::vector<std::string> mc;
  const std::string bm( beam );
  if ( bm == "fhc" || bm == "comb" )
    mc = { "Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
           "Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
           "Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
           "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc" };
  if ( bm == "rhc" || bm == "comb" ) {
    const std::vector<std::string> r = {
           "Run1_rhc_new_numi_flux_rhc_pandora_ntuple",
           "Run2_rhc_new_numi_flux_rhc_pandora_ntuple",
           "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa",
           "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab",
           "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac",
           "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad",
           "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae",
           "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple",
           "Run4b_rhc_new_numi_flux_rhc_pandora_ntuple",
           "Run4c_rhc_new_numi_flux_rhc_pandora_ntuple" };
    mc.insert( mc.end(), r.begin(), r.end() );
  }

  // The two productions carry DIFFERENT selections: processed/ holds the inclusive
  // CC1mu1piXp branches, processed/w/ the proton-tagged CC1mu1pi1p ones. Reading the wrong
  // directory gives "Bad numerical expression" on the signal branch, not an empty result.
  const char* dir = onep ? "/data/uboone/processed/w/" : "/data/uboone/processed/";
  TChain ch( "stv_tree" );
  int nf = 0;
  for ( const auto& m : mc ) {
    TString p = Form( "%sxsec-ana-%s.root", dir, m.c_str() );
    if ( gSystem->AccessPathName(p) ) { printf("  [missing] %s\n", m.c_str()); continue; }
    ch.Add( p ); ++nf;
  }
  if ( !nf ) { printf("  no inputs\n"); return; }

  // ---- one fine (true, reco) histogram --------------------------------------------
  const int NT = (int)std::lround( (HI-LO)/STEP );   // true axis starts at threshold
  const int NR = (int)std::lround(  HI    /STEP );   // reco axis starts at 0
  TString PC = Form( "sqrt(pow(sqrt(pow(%s_candidate_pion_mom_reco,2)+0.011164)"
                     "-0.10566+0.13957,2)-0.019480)", S );
  TString cut = Form( "%s_MC_Signal && %s_Selected", S, S );
  TString draw = Form( "%s : %s_candidate_pion_mom_true >> h2(%d,%g,%g,%d,%g,%g)",
                       PC.Data(), S, NT, LO, HI, NR, 0., HI );
  ch.Draw( draw, cut, "goff" );
  TH2D* h2 = (TH2D*) gDirectory->Get("h2");
  if ( !h2 || h2->GetEntries() == 0 ) { printf("  no entries\n"); return; }

  // include the reco overflow in the top reco row so nothing is silently lost
  for ( int i = 1; i <= NT; ++i )
    h2->SetBinContent( i, NR, h2->GetBinContent(i,NR) + h2->GetBinContent(i,NR+1) );

  // ---- 2D cumulative sum: C[i][j] = events with true index < i and reco index < j ----
  std::vector<std::vector<double>> C( NT+1, std::vector<double>( NR+1, 0. ) );
  for ( int i = 1; i <= NT; ++i )
    for ( int j = 1; j <= NR; ++j )
      C[i][j] = h2->GetBinContent(i,j) + C[i-1][j] + C[i][j-1] - C[i-1][j-1];
  auto box = [&]( int t0, int t1, int r0, int r1 ) {   // half-open index ranges
    return C[t1][r1] - C[t0][r1] - C[t1][r0] + C[t0][r0];
  };
  auto ti = [&]( double e ) { return (int)std::lround( (e-LO)/STEP ); };  // true index
  auto ri = [&]( double e ) { return (int)std::lround(  e    /STEP ); };  // reco index

  const double total = box( 0, NT, 0, NR );
  printf( "\n  %s, %s: %.0f selected signal events, %d files\n", S, beam, total, nf );
  printf( "  scanning edges on a %.0f MeV/c grid, minimum bin width %.0f MeV/c,"
          " criterion diagonal > %.2f\n", STEP*1000, MINW*1000, CRIT );

  // evaluate one edge set: returns per-bin diagonal and per-bin true count
  auto evaluate = [&]( const std::vector<double>& e,
                       std::vector<double>& diag, std::vector<double>& n ) {
    const int nb = e.size() - 1;
    diag.assign( nb, 0. ); n.assign( nb, 0. );
    for ( int b = 0; b < nb; ++b ) {
      int t0 = ti(e[b]), t1 = ( b==nb-1 ? NT : ti(e[b+1]) );
      // the lowest reco bin is open below the threshold
      int r0 = ( b==0 ? 0 : ri(e[b]) ), r1 = ( b==nb-1 ? NR : ri(e[b+1]) );
      double tot = box( t0, t1, 0, NR );
      n[b] = tot;
      diag[b] = tot > 0. ? box( t0, t1, r0, r1 ) / tot : 0.;
    }
  };

  // ---- reference schemes, to check this method against the numbers already quoted ----
  {
    std::vector<double> d, n;
    std::vector<std::vector<double>> refs = {
      { LO, 0.205, HI },
      { LO, 0.205, 0.235, HI },          // two 30 MeV/c bins, everything else in the last
      { LO, 0.250, 0.320, 0.420, 0.550, HI }
    };
    const char* rn[3] = { "adopted 2-bin", "3-bin 30+30+rest", "original 5-bin" };
    for ( int r = 0; r < 3; ++r ) {
      evaluate( refs[r], d, n );
      printf( "  %-15s column-norm (efficiency-like, the criterion):", rn[r] );
      for ( double x : d ) printf( " %5.1f%%", 100*x );
      printf( "\n  %-15s events per true bin                       :", "" );
      for ( double x : n ) printf( " %6.0f", x );
      printf( "\n  %-15s row-norm    (purity of each reco bin)      :", "" );
      const std::vector<double>& e = refs[r];
      const int nb = e.size()-1;
      for ( int b = 0; b < nb; ++b ) {
        int t0 = ti(e[b]), t1 = ( b==nb-1 ? NT : ti(e[b+1]) );
        int r0 = ( b==0 ? 0 : ri(e[b]) ), r1 = ( b==nb-1 ? NR : ri(e[b+1]) );
        double col = box( 0, NT, r0, r1 );
        printf( " %5.1f%%", col>0 ? 100*box(t0,t1,r0,r1)/col : 0. );
      }
      printf( "\n" );
    }
  }

  // ---- two-bin boundary trade-off, the table the note quotes -----------------------
  {
    printf( "\n  two-bin boundary trade-off (first-bin width vs diagonals)\n" );
    printf( "    boundary  width      bin 1    bin 2   events b1\n" );
    std::vector<double> d, n;
    for ( double b : { 0.185, 0.190, 0.195, 0.200, 0.205, 0.210, 0.220, 0.250, 0.300 } ) {
      evaluate( { LO, b, HI }, d, n );
      printf( "    %6.3f   %3.0f MeV/c   %5.1f%%   %5.1f%%   %6.0f\n",
              b, (b-LO)*1000, 100*d[0], 100*d[1], n[0] );
    }
  }

  // ---- scan 2, 3 and 4 bins -------------------------------------------------------
  // candidate internal edges: from LO+MINW up to 0.8 (above the saturation point extra
  // edges cannot help, and the scan demonstrates that rather than assuming it)
  std::vector<double> cand;
  for ( double e = LO+MINW; e <= 0.800001; e += STEP ) cand.push_back( e );

  for ( int nb = 2; nb <= 4; ++nb ) {
    std::vector<double> best_e; double best_min = -1., best_worst_n = 0.;
    std::vector<double> bd, bn, d, n;
    std::vector<int> idx( nb-1, 0 );
    // odometer over combinations of (nb-1) internal edges, strictly increasing
    std::function<void(int,int)> rec = [&]( int pos, int start ) {
      if ( pos == nb-1 ) {
        std::vector<double> e{ LO };
        for ( int k = 0; k < nb-1; ++k ) e.push_back( cand[idx[k]] );
        e.push_back( HI );
        for ( size_t k = 1; k < e.size(); ++k )
          if ( e[k]-e[k-1] < MINW-1e-9 ) return;
        evaluate( e, d, n );
        double mn = *std::min_element( d.begin(), d.end() );
        double wn = *std::min_element( n.begin(), n.end() );
        if ( wn < MINEVT ) return;                 // reject bins with too few events
        if ( mn > best_min ) { best_min = mn; best_e = e; bd = d; bn = n; best_worst_n = wn; }
        return;
      }
      for ( int i = start; i < (int)cand.size(); ++i ) { idx[pos] = i; rec( pos+1, i+1 ); }
    };
    rec( 0, 0 );

    printf( "\n  --- %d bins ---\n", nb );
    if ( best_min < 0. ) { printf( "    no configuration satisfies the >=%d events/bin floor\n", MINEVT ); continue; }
    printf( "    best edges :" );
    for ( size_t k = 0; k+1 < best_e.size(); ++k ) printf( " %.3f", best_e[k] );
    printf( " -> open\n    diagonals  :" );
    for ( double x : bd ) printf( " %5.1f%%", 100*x );
    printf( "   (worst %5.1f%%, %s)\n", 100*best_min, best_min > CRIT ? "PASSES" : "FAILS" );
    printf( "    events/bin :" );
    for ( double x : bn ) printf( " %6.0f", x );
    printf( "\n" );
  }
}
