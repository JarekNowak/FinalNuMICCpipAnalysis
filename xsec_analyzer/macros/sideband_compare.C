// sideband_compare.C — data/MC yields in the background-control sidebands
// (CC0pi, multi-pi, pi0, cosmic), per configuration, PER-RUN POT-scaled. Each sideband
// is signal-depleted, so this validates (and can constrain) the MC background before
// the signal region is unblinded. Reports, per sideband: N_data, the MC breakdown
// (signal contamination / nu-bkg / EXT / dirt), the data/MC normalisation ratio, and a
// scale factor for the dominant background.
//
// REQUIRES reprocessed files carrying the sb_* + EventCategory branches (add after the
// per-run univmake batch finishes: recompile, then re-run ProcessNTuples).
//
// Data source:  "fake" (blind, default) uses the per-run fake data; "beamon" uses the
// real per-run beam-on files (ONLY when explicitly unblinding the control regions).
//   usage:  root -l -b -q 'macros/sideband_compare.C("fhc","fake")'   // or rhc/comb, beamon
#include "sb_guard.h"   // blinding guard on every beam-on file (2026-09-06)
#include <vector>
#include <string>
#include "TRandom3.h"
struct Src { std::string file; double scale; };
// Per-run-period beam-off files (both selections, swtrig fail-open build) and the
// analyser's gate counts (2026-09-02). Run 2 has no beam-off sample and never will: the
// Run-1 sample stands in, scaled by its own gates.
static void add_ext(std::vector<Src>& v, const char* m){
  const char* E="/data/uboone/processed/ext_perrun/xsec-ana-";
  const char* R1="neutrinoselection_filt_run1_beamoff.root", *R3B="neutrinoselection_filt_run3b_beamoff.root";
  const char* R4[4]={"numi_pelee_ntuple_beam_off_run4a_rhc_ana.root","numi_pelee_ntuple_beam_off_run4b_rhc_ana.root",
                     "numi_pelee_ntuple_beam_off_run4c_fhc_ana.root","numi_pelee_ntuple_beam_off_run4d_fhc_ana.root"};
  const char* R5="numi_pelee_ntuple_beam_off_run5_fhc_ana.root";
  const double G1=4582248.27, G3=32649128.65, G4=34831148.625, G5=19256341.475, OCCX=0.98;
  std::string mm=m;
  if(mm=="fhc"||mm=="comb"){ v.push_back({std::string(E)+R1,OCCX*9846635./G1}); v.push_back({std::string(E)+R1,OCCX*3535129./G1});
    for(auto f:R4) v.push_back({std::string(E)+f,OCCX*4131149./G4}); v.push_back({std::string(E)+R5,OCCX*5154196./G5}); }
  if(mm=="rhc"||mm=="comb"){ v.push_back({std::string(E)+R1,OCCX*1458253./G1}); v.push_back({std::string(E)+R1,OCCX*5422907./G1});
    v.push_back({std::string(E)+R3B,OCCX*10349610./G3}); for(auto f:R4) v.push_back({std::string(E)+f,OCCX*6304167./G4}); }
}

void sideband_compare(const char* mode="fhc", const char* datasrc="fake",
                      const char* dir="/data/uboone/processed/sb/",
                      const char* sel="CC1mu1piXp", bool proton_tag=false) {
  // `sel`        : selection prefix. "CC1mu1piXp" = inclusive; "CC1mu1pi1p" = the
  //                proton-tagged subsample, whose files inherit the same sb_* branches
  //                from the parent selection, so no new definitions are needed.
  // `proton_tag` : additionally require a reconstructed tagged proton
  //                (n_proton_reco >= 1). This turns each inherited region into its
  //                proton-tagged counterpart, which is the appropriate control region
  //                for the W/TKI measurement because it lives in the same phase space.
  //                NOTE: these regions APPLY the proton tag, so they cannot constrain
  //                the proton MIS-TAG background -- that needs an inverted-PID region.
  // `dir` holds the sb_-instrumented reprocessed files (MC + fake data + EXT + dirt).
  // Defaults to the sideband reprocess dir so the running xsec batch's standard files
  // stay untouched; point it back at the standard dir once those carry sb_.
  const char* P=dir;
  std::vector<Src> mc; std::vector<std::string> data;
  double sc_ext, sc_dirt;
  auto FHCmc=[&](){ const char* rn[4]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple",
      "Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple",
      "reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"}; double sc[4]={0.14101,0.05085,0.07323,0.11560};
    for(int i=0;i<4;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+".root",sc[i]}); };
  auto RHCmc=[&](){ const char* rn[5]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
    double sc[5]={0.06728,0.04478,0.08847,0.08847,0.08847};
    for(int i=0;i<5;i++) mc.push_back({std::string(P)+"xsec-ana-"+rn[i]+"_new_numi_flux_rhc_pandora_ntuple.root",sc[i]});
    for(auto s:{"aa","ab","ac","ad","ae"}) mc.push_back({std::string(P)+"xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_"+std::string(s)+".root",0.09066}); };
  auto FHCdata=[&](){ for(auto r:{"run1","run2","run4","run5"}) data.push_back(std::string(P)+"xsec-ana-fakedata_fhc_"+r+".root"); };
  auto RHCdata=[&](){ for(auto r:{"run1","run2","run3","run4"}) data.push_back(std::string(P)+"xsec-ana-fakedata_rhc_"+r+".root"); };
  std::string m=mode;
  // sc_ext carries the 2% NuMI beam-occupancy factor the framework applies on top of the
  // beam-on/beam-off trigger ratio (SystematicsCalculator.cxx: "* 0.98"); it was missing
  // here, over-counting EXT by 2%. sc_dirt for "comb" sums both modes' exposures, matching
  // sc_ext -- the old two-branch ternary silently gave comb the FHC-only dirt scale.
  const double NUMI_EXT_OCC = 0.98;
  if(m=="fhc"){FHCmc();FHCdata();sc_ext=NUMI_EXT_OCC*5.9313;sc_dirt=0.092402*0.65;}
  else if(m=="rhc"){RHCmc();RHCdata();sc_ext=NUMI_EXT_OCC*6.1584;sc_dirt=0.071666*0.65;}
  else {FHCmc();RHCmc();FHCdata();RHCdata();sc_ext=NUMI_EXT_OCC*12.0898;sc_dirt=(0.092402+0.071666)*0.65;}
  // real beam-on files would replace `data` here when unblinding the control regions.
  if(std::string(datasrc)=="beamon"){ printf("  [beamon requested — real-data control-region unblinding; not wired until authorised]\n"); return; }
  // "fakestack": the full-stack technical test requested at review. The pseudo-data are the
  // neutrino-MC throw PLUS an independent Poisson throw (fixed seed) of the EXT and dirt
  // expectations, and are compared with the COMPLETE prediction (nu-MC + EXT + dirt) through
  // the same code path that real beam-on data will take. It tests bookkeeping, not physics.
  const bool fakestack = ( std::string(datasrc)=="fakestack" );
  TRandom3 rng(20260903);

  // Weighted event count passing `cut`. MC (numuMC + dirt) carry the CV weight
  // (tuned_cv x ppfx_cv x normalisation) -- the SAME weight the analysis and the
  // fake-data throw apply -- so a bare GetEntries would under/over-count by the mean
  // CV weight. EXT and (fake) data are unweighted (data-driven / weight-1 throws).
  // Per-event finite/non-negative guard on the CV weight, matching the throw macro's
  // `if(!isfinite(cv)||cv<0) continue` -- otherwise a single inf/nan weight makes the
  // whole Draw sum NaN.
  // Framework safe_weight rule: non-finite, negative or >30 -> 1 (one authoritative rule, 2026-09-03)
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)<=30?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:1)";
  auto wsum=[&](TChain& c,const TString& cut,bool weighted)->double{
    TH1D h("h_ws","",1,-0.5,1.5);
    c.Draw("0.5>>h_ws", (weighted?TString(CVW):TString("1"))+"*("+cut+")","goff");
    double v=h.Integral(0,2); return std::isfinite(v)?v:0.; };

  const char* sb[4]={"sb_cc0pi","sb_multipi","sb_pi0","sb_cosmic"};
  TString PRE = TString(sel)+"_";
  TString PTAG = proton_tag ? TString(" && ")+PRE+"n_proton_reco>=1" : TString("");
  printf("\n==== %s sidebands [%s%s] (data source: %s, per-run scaled, CV-weighted MC) ====\n",
         mode, sel, proton_tag?" + proton tag":"", datasrc);
  printf("%-12s %9s %9s %9s %9s %9s %8s %9s\n","sideband","N_data","sig","nu-bkg","EXT","dirt","sig%%","data/vMC");
  for(auto s:sb){
    TString F=TString("(")+PRE+s+PTAG+")", SIG=PRE+"MC_Signal";
    double Nsig=0,Nbkg=0;
    for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str());
      Nsig+=wsum(c,F+" && "+SIG,true)*x.scale; Nbkg+=wsum(c,F+" && !"+SIG,true)*x.scale; }
    // EXT per RUN PERIOD (2026-09-02): cosmic-only beam-off matched by period, horn label
    // irrelevant; each run's files scaled by bnb_gates[run]/sum(file gates) x 0.98 (see
    // add_ext), exactly as the framework does. The per-run files carry both selections'
    // sb_ flags. sc_ext (pooled sample) is no longer used for EXT.
    std::vector<Src> extv; add_ext(extv, m.c_str());
    double Next=0; for(auto&x:extv){ TChain ce("stv_tree"); ce.Add(x.file.c_str()); Next+=wsum(ce,F,false)*x.scale; }
    TChain cd("stv_tree"); cd.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P));
    double Ndirt=wsum(cd,F,true)*sc_dirt;
    double Ndata=0; for(auto&d:data){ sb_guard_data(d); TChain c("stv_tree"); c.Add(d.c_str()); Ndata+=wsum(c,F,false); }
    // The fake data is a Poisson throw of the neutrino MC only, so the machinery/
    // normalisation check compares it to the CV-weighted neutrino MC (sig+nubkg).
    // EXT + dirt are the additional components real beam-on data will contain (shown
    // for scale); their model is validated against the FULL stack only after unblinding.
    double Nnu=Nsig+Nbkg;
    if(fakestack){ double Nfull=Nnu+Next+Ndirt; double Ndat2=Ndata+rng.Poisson(Next)+rng.Poisson(Ndirt);
      printf("%-12s %9.1f %9.1f %9.1f %9.1f %9.1f %7.1f%% %9.3f  [full stack %9.1f, data/full %.3f +- %.3f]\n",
             s,Ndat2,Nsig,Nbkg,Next,Ndirt,Nnu>0?100*Nsig/Nnu:0,Nfull>0?Ndat2/Nfull:0,Nfull,Nfull>0?Ndat2/Nfull:0,Nfull>0?sqrt(Ndat2)/Nfull:0); continue; }
    printf("%-12s %9.1f %9.1f %9.1f %9.1f %9.1f %7.1f%% %9.3f\n",
           s,Ndata,Nsig,Nbkg,Next,Ndirt,Nnu>0?100*Nsig/Nnu:0,Nnu>0?Ndata/Nnu:0);
  }
  printf("(FAKE data = Poisson throw of the neutrino MC CV; data/vMC ~1 validates the\n"
         " per-run weighting/scaling machinery, and sig%% (of vMC) confirms signal depletion.\n"
         " EXT+dirt are the extra real-data components -- data/full-MC (background model)\n"
         " is only meaningful once the control regions are unblinded with beam-on data.)\n");
}
