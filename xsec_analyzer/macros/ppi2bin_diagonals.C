// ppi2bin_diagonals.C -- migration diagonal for the two-bin p_pi scheme.
//
// For each TRUE bin, the fraction of SELECTED SIGNAL events whose reconstructed p_pi lands
// in the matching reco bin. This is the number the BNB CC1pi criterion is stated in (a bin
// is usable above 0.68), and it is what decided both the two-bin scheme and the move of
// p_pi onto the proton-tagged selection.
//
// The reco momentum is mass-corrected inline exactly as the bin configs do it (range
// momentum under the MUON hypothesis -> kinetic energy -> pion). Using the raw branch
// instead overstates the bias; that has caught us out before.
//
// Measured results (2026-08-17/18):
//   inclusive CC1mu1piXp, FHC : 97 / 62 %   -> bin 2 FAILS
//   proton-tagged CC1mu1pi1p, FHC : 88.9 / 74.2 %  -> both pass
//
//   usage: root -l -b -q 'macros/ppi2bin_diagonals.C("1p","fhc")'
//          first arg  "1p" (CC1mu1pi1p) or "incl" (CC1mu1piXp)
//          second arg "fhc" or "rhc"

void ppi2bin_diagonals( const char* sel = "1p", const char* beam = "fhc" ) {

  const bool onep = ( std::string(sel) == "1p" );
  const char* S = onep ? "CC1mu1pi1p" : "CC1mu1piXp";
  const double EDGE = 0.205, LO = 0.175;

  // NOTE: the Run5 FHC overlay does NOT follow the Run1/2/4 naming convention. A glob built
  // on that pattern silently drops it and shifts the answer (90.3/74.4 from three runs
  // against 88.9/74.2 from four), so the files are listed explicitly.
  std::vector<std::string> mc;
  const std::string bm( beam );
  if ( bm == "fhc" || bm == "comb" ) {
    mc = { "Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
           "Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
           "Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
           "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc" };
  }
  if ( bm == "rhc" || bm == "comb" ) {
    // combined = the union of both beam file sets, matching file_properties_numi_comb_w.txt
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

  TChain ch( "stv_tree" );
  int nf = 0;
  for ( const auto& m : mc ) {
    // processed/ holds the inclusive branches, processed/w/ the proton-tagged ones
    const char* dir = onep ? "/data/uboone/processed/w/" : "/data/uboone/processed/";
    TString p = Form( "%sxsec-ana-%s.root", dir, m.c_str() );
    if ( gSystem->AccessPathName(p) ) { printf( "  [missing] %s\n", p.Data() ); continue; }
    ch.Add( p ); ++nf;
  }
  if ( !nf ) { printf( "  no input files\n" ); return; }

  TString PC = Form( "sqrt(pow(sqrt(pow(%s_candidate_pion_mom_reco,2)+0.011164)"
                     "-0.10566+0.13957,2)-0.019480)", S );
  TString base = Form( "%s_MC_Signal && %s_Selected", S, S );
  TString t1 = Form( "%s && %s_candidate_pion_mom_true >= %g && %s_candidate_pion_mom_true < %g",
                     base.Data(), S, LO, S, EDGE );
  TString t2 = Form( "%s && %s_candidate_pion_mom_true >= %g", base.Data(), S, EDGE );

  long n1  = ch.GetEntries( t1 );
  long n1d = ch.GetEntries( Form("%s && %s <  %g", t1.Data(), PC.Data(), EDGE) );
  long n2  = ch.GetEntries( t2 );
  long n2d = ch.GetEntries( Form("%s && %s >= %g", t2.Data(), PC.Data(), EDGE) );

  printf( "\n  %s, %s, two-bin p_pi split at %.3f GeV/c (%d files)\n",
          S, beam, EDGE, nf );
  printf( "    true bin 1 [%.3f,%.3f): %6ld selected signal, diagonal %5.1f%%%s\n",
          LO, EDGE, n1, n1 ? 100.*n1d/n1 : 0., (n1 && 100.*n1d/n1 > 68.) ? "  PASS" : "  FAIL" );
  printf( "    true bin 2 (> %.3f)   : %6ld selected signal, diagonal %5.1f%%%s\n",
          EDGE, n2, n2 ? 100.*n2d/n2 : 0., (n2 && 100.*n2d/n2 > 68.) ? "  PASS" : "  FAIL" );
  printf( "\n    criterion: a bin is usable above 68%%\n" );
}
