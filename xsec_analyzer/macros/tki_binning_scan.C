// tki_binning_scan.C -- re-derive the binning for the four transverse-kinematic-imbalance
// observables after the beam-frame correction (NuMIBeamFrame.hh). The previous edges were
// equal-population quantiles of the uncorrected, detector-z variables and are unusable on
// the corrected ones: delta_pT would put 80.8% of events in the first bin and 0.2% in the
// last, delta_phiT 63.6% in the first.
//
// Same criterion as macros/ppi_binning_scan.C and macros/w_binning_scan.C -- migration
// diagonal above 0.68 with a floor on the events per bin -- but the search is a maximin
// dynamic program rather than an exhaustive edge scan, so it reaches six bins instead of
// four. best[k][i] = max over j < i of min(best[k-1][j], diag(j..i)); the first bin's reco
// lower edge and the last bin's reco upper edge are open, matching the bin configs.
//
// The events floor is deliberately stricter than the 100 used for W_pipr: that scheme
// passed the diagonal criterion with 137 events in its first bin and still failed the
// extraction, A_C collapsing to a rank-one matrix, while the 469-event p_pi case worked.
//
// Extended 2026-09-16 for the 0.50-criterion study: the scan also covers W_pipr, W_had and
// the proton angles (cos theta_p, theta_p, and the pion-proton opening angle theta_pipr from
// the friend trees of macros/proton_angle_friend.C), and ONLY selects a subset by name.
// Results are printed against both 0.68 and CRIT; the scan still stops at the largest bin
// count passing CRIT.
//   usage: root -l -b -q 'macros/tki_binning_scan.C("/data/uboone/processed/beta")'
//          root -l -b -q 'macros/tki_binning_scan.C("/data/uboone/processed/w",400,0.50,"thpipr,costhp")'
#include <vector>
#include <algorithm>
#include <string>

namespace {
  // br: CC1mu1pi1p_<br>_{reco,true}; when rexp/texp are given they are used verbatim instead
  struct Spec { const char* br; const char* name; double lo, hi, step, minw;
                const char* rexp = nullptr; const char* texp = nullptr; };
}

void tki_binning_scan( const char* dir = "/data/uboone/processed/beta",
                       int MINEVT = 400, double CRIT = 0.68, const char* ONLY = "", int KMAX = 6 ) {

  const char* S = "CC1mu1pi1p";
  std::vector<Spec> obs = {
    { "deltaPt",     "delta_pT  [GeV/c]", 0.,   1.50, 0.025, 0.075 },
    { "deltaAlphaT", "delta_alphaT [deg]", 0., 180.,  5.,   15.    },
    { "deltaPhiT",   "delta_phiT   [deg]", 0., 180.,  5.,   15.    },
    { "pn",          "p_n       [GeV/c]", 0.,   1.50, 0.025, 0.075 },
    { "W_pipr",      "W_pipr  [GeV/c2]",  1.08, 2.60, 0.01,  0.05  },
    { "W_had",       "W_had   [GeV/c2]",  0.,   2.74, 0.02,  0.10  },
    { "costhp",      "cos theta_p",       -1.,  1.,   0.025, 0.10,
      "CC1mu1pi1p_proton_costh_reco", "CC1mu1pi1p_proton_costh_true" },
    { "thetap",      "theta_p [rad]",      0.,  3.1416, 0.0314, 0.157,
      "TMath::ACos(CC1mu1pi1p_proton_costh_reco)", "TMath::ACos(CC1mu1pi1p_proton_costh_true)" },
    { "thpipr",      "theta_pipr [rad]",   0.,  3.1416, 0.0314, 0.157,
      "pa.thpipr_reco", "pa.thpipr_true" },
  };

  TChain ch("stv_tree");
  int nf = ch.Add( Form("%s/xsec-ana-Run*_fhc_*.root", dir) );
  printf("\n  chained %d file(s) from %s\n", nf, dir);
  if ( ch.GetEntries() <= 0 ) { printf("  no entries -- has stage 1 finished?\n"); return; }
  // friend trees (theta_pipr), file by file in the same order as the main chain
  TChain fr("pa"); int nfr = 0;
  for ( auto* el : *ch.GetListOfFiles() ) {
    TString fn = el->GetTitle(); TString b = gSystem->BaseName( fn ); b.ReplaceAll( ".root", ".pa.root" );
    TString ff = TString( gSystem->DirName( fn ) ) + "/friends_pa/" + b;
    if ( !gSystem->AccessPathName( ff ) ) { fr.Add( ff ); ++nfr; }
  }
  if ( nfr == nf ) { ch.AddFriend( &fr, "pa" ); printf("  friend trees attached (%d)\n", nfr); }
  std::string only = std::string(",") + ONLY + ",";

  for ( auto& o : obs ) {
    if ( only != ",," && only.find( std::string(",")+o.br+"," ) == std::string::npos ) continue;
    if ( std::string(o.br)=="thpipr" && nfr != nf ) { printf("\n  %s: friend trees missing, skipped\n", o.name); continue; }
    int N = (int)std::lround( (o.hi-o.lo)/o.step );
    std::string rx = o.rexp ? o.rexp : Form("%s_%s_reco",S,o.br);
    std::string tx = o.texp ? o.texp : Form("%s_%s_true",S,o.br);
    gDirectory->Delete("h2;*");
    ch.Draw( Form("%s : %s >> h2(%d,%g,%g,%d,%g,%g)",
                  rx.c_str(),tx.c_str(),N,o.lo,o.hi,N,o.lo,o.hi),
             Form("%s_Selected && %s_MC_Signal",S,S), "goff" );
    TH2D* h2 = (TH2D*)gDirectory->Get("h2");
    if ( !h2 || h2->GetEntries()<=0 ) { printf("\n  %s: no entries\n", o.name); continue; }

    // 2D cumulative sums for O(1) box queries
    std::vector<std::vector<double>> C(N+1, std::vector<double>(N+1,0.));
    for ( int i=1;i<=N;++i ) for ( int j=1;j<=N;++j )
      C[i][j]=h2->GetBinContent(i,j)+C[i-1][j]+C[i][j-1]-C[i-1][j-1];
    auto box=[&](int t0,int t1,int r0,int r1){ return C[t1][r1]-C[t0][r1]-C[t1][r0]+C[t0][r0]; };

    double tot = box(0,N,0,N);
    printf("\n  === %s : %.0f selected signal events ===\n", o.name, tot);

    // diagonal and population of the true-bin [t0,t1)
    auto diag=[&](int t0,int t1,bool first,bool last,double& nev){
      nev = box(t0,t1,0,N);
      if ( nev <= 0 ) return 0.;
      int r0 = first ? 0 : t0, r1 = last ? N : t1;
      return box(t0,t1,r0,r1)/nev;
    };
    int MINW = (int)std::lround(o.minw/o.step);

    for ( int K=KMAX; K>=2; --K ) {
      // best[k][i] : best achievable worst-diagonal using k bins covering [0,i)
      std::vector<std::vector<double>> best(K+1, std::vector<double>(N+1,-1.));
      std::vector<std::vector<int>>    prev(K+1, std::vector<int>(N+1,-1));
      for ( int i=MINW; i<=N; ++i ) {
        double nev; double d = diag(0,i,true,i==N,nev);
        if ( nev >= MINEVT ) { best[1][i]=d; prev[1][i]=0; }
      }
      for ( int k=2;k<=K;++k )
        for ( int i=k*MINW; i<=N; ++i )
          for ( int j=(k-1)*MINW; j<=i-MINW; ++j ) {
            if ( best[k-1][j] < 0 ) continue;
            double nev; double d = diag(j,i,false,i==N,nev);
            if ( nev < MINEVT ) continue;
            double v = std::min(best[k-1][j], d);
            if ( v > best[k][i] ) { best[k][i]=v; prev[k][i]=j; }
          }

      if ( best[K][N] < 0 ) { printf("    %d bins : impossible at >=%d events/bin\n", K, MINEVT); continue; }
      std::vector<int> cut; int i=N;
      for ( int k=K;k>=1;--k ){ cut.push_back(i); i=prev[k][i]; }
      cut.push_back(0); std::reverse(cut.begin(),cut.end());
      printf("    %d bins : worst diagonal %5.1f%%  %s at %.2f, %s at 0.68\n", K, 100*best[K][N],
             best[K][N] > CRIT ? "PASSES" : "fails", CRIT, best[K][N] > 0.68 ? "passes" : "fails");
      printf("        edges  :"); for ( size_t k=1;k+1<cut.size();++k ) printf(" %7.3f", o.lo+cut[k]*o.step); printf("\n");
      printf("        diag   :");
      for ( size_t k=0;k+1<cut.size();++k ){ double nev; printf(" %6.1f%%",100*diag(cut[k],cut[k+1],k==0,k+2==cut.size(),nev)); }
      printf("\n        events :");
      for ( size_t k=0;k+1<cut.size();++k ){ double nev; diag(cut[k],cut[k+1],k==0,k+2==cut.size(),nev); printf(" %7.0f",nev); }
      printf("\n");
      if ( best[K][N] > CRIT ) break;   // the largest K that passes is what we want
    }
  }
}
