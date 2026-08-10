// total_xsec_counting.C — cut-and-count TOTAL flux-averaged cross section, per
// configuration, WITHOUT unfolding. For each config:
//   sigma = (N_sel - N_bkg) / (eff * conv_factor)   [1e-38 cm^2 / Ar]
// with conv_factor = num_Ar * (data_POT * flux_per_pot) / 1e38 (same normalisation
// the framework uses for the differential result), eff = N_sig_sel/N_sig_gen, and
// N_sel = selected data (blind: the per-run fake data == N_sig_sel + N_bkg, so the
// count returns the MC-tune total sigma; at unblinding N_sel is the real beam-on).
// PER-RUN POT scaling (Table tab:pot). Run:  root -l -b -q 'macros/total_xsec_counting.C'
#include <vector>
#include <string>
struct Src { std::string file; double scale; };

double num_Ar_FV() {
  // NuMI FV X[10,246] Y[-101,101] Z[10,986] minus dead slab z[675.1,775.1]
  double vol = (246.-10.)*(101.-(-101.))*(986.-10.)
             - (246.-10.)*(101.-(-101.))*(775.1-675.1);         // cm^3
  return vol * 1.3836 * 6.02214076e23 / 39.948;                  // Ar nuclei
}

void one(const char* mode, double dataPOT, double flux_per_pot,
         std::vector<Src>& mc, double sc_ext, double sc_dirt,
         double syst_frac, const char* gtag) {
  const char* P="/data/uboone/processed/";
  TString SEL="CC1mu1piXp_Selected", SIG="CC1mu1piXp_MC_Signal";
  double Nsig_gen=0, Nsig_sel=0, Nbkg_mc=0;
  for (auto& s : mc) {
    TChain c("stv_tree"); c.Add(s.file.c_str());
    Nsig_gen += c.GetEntries(SIG) * s.scale;
    Nsig_sel += c.GetEntries(SEL+" && "+SIG) * s.scale;
    Nbkg_mc  += c.GetEntries(SEL+" && !"+SIG) * s.scale;
  }
  TChain ce("stv_tree"); ce.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
  double Next = ce.GetEntries(SEL) * sc_ext;
  TChain cd("stv_tree"); cd.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P));
  double Ndirt = cd.GetEntries(SEL) * sc_dirt;

  double Nbkg = Nbkg_mc + Next + Ndirt;
  double Nsel = Nsig_sel + Nbkg;                 // blind: fake data = prediction
  double eff  = Nsig_sel / Nsig_gen;
  double conv = num_Ar_FV() * (dataPOT * flux_per_pot) / 1e38;   // events per 1e-38 cm^2/Ar
  double sigma = (Nsel - Nbkg) / (eff * conv);
  double sigma_mc = Nsig_gen / conv;             // == sigma for fake data (cross-check)

  printf("\n=== %s ===  num_Ar=%.4e  conv=%.2f  (dataPOT=%.4e flux/POT=%.5e)\n",
         mode, num_Ar_FV(), conv, dataPOT, flux_per_pot);
  printf("  N_sig_gen=%.1f  N_sig_sel=%.1f  eff=%.4f\n", Nsig_gen, Nsig_sel, eff);
  printf("  N_bkg=%.1f (MC %.1f + EXT %.1f + dirt %.1f)  N_sel=%.1f\n",
         Nbkg, Nbkg_mc, Next, Ndirt, Nsel);
  // uncertainty: systematic (flux/detector/xsec/reint/MC-stat/EXT-stat/POT, from the
  // systematic breakdown, extraction-method independent) (+) cut-and-count data stat.
  double dstat = sqrt(Nsel) / (Nsel - Nbkg);
  double dtot  = sqrt(syst_frac*syst_frac + dstat*dstat);
  printf("  sigma_count = %.4f +/- %.4f  (syst %.1f%% (+) data-stat %.1f%% = %.1f%%)"
         "  [1e-38 cm^2/Ar]\n", sigma, sigma*dtot, 100*syst_frac, 100*dstat, 100*dtot);

  // compare to the standalone generator flux-averaged totals (integral of pmu_fte)
  const char* GP="../generator_predictions/newg4/";
  const char* gens[4]={"genie","gibuu","neut","nuwro"};
  printf("  generators [1e-38 cm^2/Ar]:");
  for(int g=0;g<4;g++){ TFile*fg=TFile::Open(Form("%s%s_%s_fte.root",GP,gens[g],gtag));
    if(fg&&!fg->IsZombie()){ TH1D*h=(TH1D*)fg->Get("pmu_fte");
      if(h){ double sg=h->Integral(); printf("  %s=%.3f(%.1f sig)",gens[g],sg,(sigma-sg)/(sigma*dtot)); } fg->Close(); } }
  printf("\n");
}

void total_xsec_counting() {
  const char* P="/data/uboone/processed/";
  // per-run POT scales D_run/MC_run (Table tab:pot)
  auto FHC=[&](std::vector<Src>&mc){
    const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple",
      "Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
    double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root", sc[i]}); };
  auto RHC=[&](std::vector<Src>&mc){
    const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
    double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root", sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"})
      mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root", 0.09066}); };

  std::vector<Src> fhc, rhc, comb;
  FHC(fhc); RHC(rhc); FHC(comb); RHC(comb);
  // dirt scaled to the full mode exposure (D_total/summedMC*0.65), as in cutflow_yields.C.
  // syst_frac = quadrature of the systematic-breakdown sources (excl. data stat), pmu column
  //   (Tables tab:systbreak_*): FHC 25.5%, RHC 35.7%, comb 28.9%.
  one("FHC",      8.857e20,  6.81159e-10, fhc,  5.9313,  0.092402*0.65, 0.255, "newg4");
  one("RHC",      1.1082e21, 6.44646e-10, rhc,  6.1584,  0.071666*0.65, 0.357, "rhc");
  one("Combined", 1.99390e21,6.60865e-10, comb, 12.0898, 0.092402*0.65, 0.289, "comb");
}
