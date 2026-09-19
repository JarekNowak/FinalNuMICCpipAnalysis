// total_xsec_counting.C — cut-and-count TOTAL flux-averaged cross section, per
// configuration, WITHOUT unfolding, for the inclusive (CC1mu1piXp) and the
// proton-tagged (CC1mu1pi1p) selections. For each config:
//   sigma = (N_sel - N_bkg) / (eff * conv_factor)   [1e-38 cm^2 / Ar]
// with conv_factor = num_Ar * (data_POT * flux_per_pot) / 1e38 (same normalisation
// the framework uses for the differential result), eff = N_sig_sel/N_sig_gen, and
// N_sel = selected data (blind: the per-run fake data == N_sig_sel + N_bkg, so the
// count returns the MC-tune total sigma; at unblinding N_sel is the real beam-on).
// PER-RUN POT scaling (Table tab:pot). MC counts carry the framework's central-value
// weight tune*ppfx_cv*normalisation (see CVW below), sanitised exactly as the
// framework's safe_weight() does (non-finite, negative or >30 -> 1).
// The framework's CV weight for NuMI is tuned_cv_weight * ppfx_cv_weight *
// normalisation_weight (UniverseMaker.cxx, NuMI branch: spline is 1 in these
// ntuples). The PPFX CV averages 0.944 on selected signal, so leaving it out (as the
// earlier version of this macro and the selection's cut-flow histograms did) overstates
// every yield and the cut-and-count total by ~5-6% relative to the unfolding.
static const char* CVW = "((TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight) && tuned_cv_weight*ppfx_cv_weight*normalisation_weight>=0 && tuned_cv_weight*ppfx_cv_weight*normalisation_weight<=30) ? tuned_cv_weight*ppfx_cv_weight*normalisation_weight : 1.0)";
//   Run:  root -l -b -q 'macros/total_xsec_counting.C'
#include <vector>
#include <string>
struct Src { std::string file; double scale; };

double num_Ar_FV() {
  // NuMI FV X[10,246] Y[-101,101] Z[10,986] minus dead slab z[675.1,775.1]
  double vol = (246.-10.)*(101.-(-101.))*(986.-10.)
             - (246.-10.)*(101.-(-101.))*(775.1-675.1);         // cm^3
  return vol * 1.3836 * 6.02214076e23 / 39.948;                  // Ar nuclei
}

// weighted counts of (selected, signal) combinations in one scan of the tree:
//   code = 2*SIG + SEL  ->  0 none, 1 sel&!sig, 2 sig&!sel, 3 sig&sel
static void count_mc( const std::string& file, const TString& SEL, const TString& SIG,
                      double scale, double& gen, double& sel_sig, double& bkg ) {
  TChain c("stv_tree"); c.Add(file.c_str());
  TH1D h("h_code","",4,-0.5,3.5);
  c.Draw("2*("+SIG+")+("+SEL+")>>h_code", CVW, "goff");
  gen     += ( h.GetBinContent(3) + h.GetBinContent(4) ) * scale;
  sel_sig += h.GetBinContent(4) * scale;
  bkg     += h.GetBinContent(2) * scale;
}

struct ExtHelperDummy{};
  // EXT (beam-off cosmic), PER RUN PERIOD (2026-09-02, stage 2): each run's beam-off files
  // (cosmic-only, horn label irrelevant) are scaled by bnb_gates[run] / sum of the files'
  // gate counts (analyser's table) x the 0.98 NuMI occupancy factor -- exactly the
  // framework's per-run scaling. Run 2 has no processed beam-off: the Run-1 file stands in.
  //   gates: run1 4582248.27 | run3b 32649128.65 | run4 (4a+4b 9338778.225, 4c 8060024.70,
  //          4d 17432345.70) = 34831148.625 | run5 19256341.475
  //   beam-on triggers: FHC run1 5748692 run2 3535129 run4 4131149 run5 5154196;
  //                     RHC run1 1458253 run2 5422907 run3 10349610 run4 6304167
  const char* E="/data/uboone/processed/ext_perrun/xsec-ana-";
  const char* R1=  "neutrinoselection_filt_run1_beamoff.root";
  const char* R3B= "neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root",
                     "numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5=  "numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCCX=0.98;
static void add_ext(std::vector<Src>& v, const char* m){
    if(std::string(m)=="fhc"||std::string(m)=="comb"){
      v.push_back({std::string(E)+R1, OCCX*5748692./G1});   // run 1
      v.push_back({std::string(E)+R1, OCCX*3535129./G1});   // run 2 (stand-in)
      for(auto f:R4) v.push_back({std::string(E)+f, OCCX*4131149./G4});
      v.push_back({std::string(E)+R5, OCCX*5154196./G5});
    }
    if(std::string(m)=="rhc"||std::string(m)=="comb"){
      v.push_back({std::string(E)+R1, OCCX*1458253./G1});
      v.push_back({std::string(E)+R1, OCCX*5422907./G1});
      v.push_back({std::string(E)+R3B, OCCX*10349610./G3});
      for(auto f:R4) v.push_back({std::string(E)+f, OCCX*6304167./G4});
    }
  }

void one(const char* mode, const char* selname, const char* P, double dataPOT,
         double flux_per_pot, std::vector<Src>& mc, double sc_ext, double sc_dirt,
         double syst_frac, const char* gtag) {
  TString SEL = TString(selname)+"_Selected", SIG = TString(selname)+"_MC_Signal";
  double Nsig_gen=0, Nsig_sel=0, Nbkg_mc=0;
  for (auto& s : mc) count_mc( s.file, SEL, SIG, s.scale, Nsig_gen, Nsig_sel, Nbkg_mc );
  std::vector<Src> extv; add_ext(extv, std::string(mode)=="FHC"?"fhc":(std::string(mode)=="RHC"?"rhc":"comb"));
  double Next = 0.;
  for (auto& s : extv) { TChain ce("stv_tree"); ce.Add(s.file.c_str()); Next += ce.GetEntries(SEL) * s.scale; }
  double Ndirt = 0., dummy1 = 0., dummy2 = 0.;
  // dirt: every selected dirt event is background (no signal in dirt), weighted like MC
  {
    TChain cd("stv_tree"); cd.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P));
    TH1D h("h_dirt","",1,0,2);
    cd.Draw("1>>h_dirt", TString(CVW)+"*("+SEL+")", "goff");
    Ndirt = h.GetBinContent(1) * sc_dirt;
  }

  double Nbkg = Nbkg_mc + Next + Ndirt;
  double Nsel = Nsig_sel + Nbkg;                 // blind: fake data = prediction
  double eff  = Nsig_sel / Nsig_gen;
  double conv = num_Ar_FV() * (dataPOT * flux_per_pot) / 1e38;   // events per 1e-38 cm^2/Ar
  double sigma = (Nsel - Nbkg) / (eff * conv);
  double sigma_mc = Nsig_gen / conv;             // == sigma for fake data (cross-check)

  printf("\n=== %s %s ===  num_Ar=%.4e  conv=%.2f  (dataPOT=%.4e flux/POT=%.5e)\n",
         selname, mode, num_Ar_FV(), conv, dataPOT, flux_per_pot);
  printf("  N_sig_gen=%.1f  N_sig_sel=%.1f  eff=%.4f\n", Nsig_gen, Nsig_sel, eff);
  printf("  N_bkg=%.1f (MC %.1f + EXT %.1f + dirt %.1f)  N_sel=%.1f  purity=%.3f\n",
         Nbkg, Nbkg_mc, Next, Ndirt, Nsel, Nsig_sel/Nsel);
  // uncertainty: systematic (prediction total of the released p_mu extraction, i.e. the
  // unfolded flux/detector/xsec/reint/MC-stat/EXT-stat/POT quadrature; conservative for a
  // cut-and-count, since it carries the unfolding amplification) (+) data stat.
  double dstat = sqrt(Nsel) / (Nsel - Nbkg);
  double dtot  = sqrt(syst_frac*syst_frac + dstat*dstat);
  printf("  sigma_count = %.4f +/- %.4f  (syst %.1f%% (+) data-stat %.1f%% = %.1f%%)"
         "  [1e-38 cm^2/Ar]   (N_sig_gen/conv = %.4f)\n", sigma, sigma*dtot, 100*syst_frac,
         100*dstat, 100*dtot, sigma_mc);

  // compare to the standalone generator flux-averaged totals (integral of pmu_fte)
  if ( gtag && gtag[0] ) {
    const char* GP="../generator_predictions/newg4/";
    const char* gens[4]={"genie","gibuu","neut","nuwro"};
    printf("  generators [1e-38 cm^2/Ar]:");
    for(int g=0;g<4;g++){ TFile*fg=TFile::Open(Form("%s%s_%s_fte.root",GP,gens[g],gtag));
      if(fg&&!fg->IsZombie()){ TH1D*h=(TH1D*)fg->Get("pmu_fte");
        if(h){ double sg=h->Integral(); printf("  %s=%.3f(%.1f sig)",gens[g],sg,(sigma-sg)/(sigma*dtot)); } fg->Close(); } }
    printf("\n");
  }
}

void total_xsec_counting() {
  // Two processed trees: the inclusive selection lives in processed/, the proton-tagged
  // one in processed/w/ (different POT-storage conventions, irrelevant here because the
  // per-run scale factors are applied explicitly).
  const char* P_incl="/data/uboone/processed/";
  const char* P_1p  ="/data/uboone/processed/w/";
  // per-run POT scales D_run/MC_run (Table tab:pot)
  auto FHC=[&](std::vector<Src>&mc, const char* P){
    const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
      "Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
    double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root", sc[i]}); };
  auto RHC=[&](std::vector<Src>&mc, const char* P){
    const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
    double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root", sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"})
      mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root", 0.09066}); };

  // dirt scaled to the full mode exposure (D_total/summedMC*0.65), as in cutflow_yields.C.
  // EXT scale = beam-on/beam-off GATE ratio x the 2% NuMI beam-occupancy factor. The pooled
  // beam-off file is the Run-1 FHC + Run-3b RHC samples, 4582248.27 + 32649128.65 =
  // 37231376.92 gates (analyser's table, 2026-09-02); the earlier 3821593 was an event
  // count and over-scaled EXT ~9.7x. FHC 22667109/G, RHC 23534937/G, comb 46202046/G.
  // Dirt for "Combined" sums both modes' exposures.
  const double OCC = 0.98;
  // syst_frac = prediction-total fractional uncertainty of the released p_mu extraction
  // (report/current_results.tsv, PredTotal_pct; the same numbers as the note's
  // per-configuration result tables): incl 41.1/44.5/40.6 %, 1p 49.3/37.6/42.0 %.
  const char* sels[2] = {"CC1mu1piXp","CC1mu1pi1p"};
  const char* dirs[2] = {P_incl, P_1p};
  const double sf[2][3] = {{0.411,0.445,0.406},{0.493,0.376,0.420}};
  for ( int k = 0; k < 2; ++k ) {
    std::vector<Src> fhc, rhc, comb;
    FHC(fhc,dirs[k]); RHC(rhc,dirs[k]); FHC(comb,dirs[k]); RHC(comb,dirs[k]);
    const bool incl = ( k == 0 );
    one("FHC",      sels[k], dirs[k], 8.857e20,  6.81159e-10, fhc,  OCC*0.60882,  0.092402*0.65, sf[k][0], incl?"newg4":"");
    one("RHC",      sels[k], dirs[k], 1.1082e21, 6.44646e-10, rhc,  OCC*0.63213,  0.071666*0.65, sf[k][1], incl?"rhc":"");
    one("Combined", sels[k], dirs[k], 1.99390e21,6.60865e-10, comb, OCC*1.24094, (0.092402+0.071666)*0.65, sf[k][2], incl?"comb":"");
  }
}
