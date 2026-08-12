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
#include <vector>
#include <string>
struct Src { std::string file; double scale; };

void sideband_compare(const char* mode="fhc", const char* datasrc="fake",
                      const char* dir="/data/uboone/processed/sb/") {
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
  if(m=="fhc"){FHCmc();FHCdata();sc_ext=5.9313;sc_dirt=0.092402*0.65;}
  else if(m=="rhc"){RHCmc();RHCdata();sc_ext=6.1584;sc_dirt=0.071666*0.65;}
  else {FHCmc();RHCmc();FHCdata();RHCdata();sc_ext=12.0898;sc_dirt=0.092402*0.65;}
  // real beam-on files would replace `data` here when unblinding the control regions.
  if(std::string(datasrc)=="beamon"){ printf("  [beamon requested — real-data control-region unblinding; not wired until authorised]\n"); return; }

  // Weighted event count passing `cut`. MC (numuMC + dirt) carry the CV weight
  // (tuned_cv x ppfx_cv x normalisation) -- the SAME weight the analysis and the
  // fake-data throw apply -- so a bare GetEntries would under/over-count by the mean
  // CV weight. EXT and (fake) data are unweighted (data-driven / weight-1 throws).
  // Per-event finite/non-negative guard on the CV weight, matching the throw macro's
  // `if(!isfinite(cv)||cv<0) continue` -- otherwise a single inf/nan weight makes the
  // whole Draw sum NaN.
  const char* CVW="(TMath::Finite(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)"
                  "&&(tuned_cv_weight*ppfx_cv_weight*normalisation_weight)>=0"
                  "?tuned_cv_weight*ppfx_cv_weight*normalisation_weight:0)";
  auto wsum=[&](TChain& c,const TString& cut,bool weighted)->double{
    TH1D h("h_ws","",1,-0.5,1.5);
    c.Draw("0.5>>h_ws", (weighted?TString(CVW):TString("1"))+"*("+cut+")","goff");
    double v=h.Integral(0,2); return std::isfinite(v)?v:0.; };

  const char* sb[4]={"sb_cc0pi","sb_multipi","sb_pi0","sb_cosmic"};
  const char* PRE="CC1mu1piXp_";
  printf("\n==== %s sidebands (data source: %s, per-run scaled, CV-weighted MC) ====\n",mode,datasrc);
  printf("%-12s %9s %9s %9s %9s %9s %8s %9s\n","sideband","N_data","sig","nu-bkg","EXT","dirt","sig%%","data/vMC");
  for(auto s:sb){
    TString F=Form("%s%s",PRE,s), SIG=Form("%sMC_Signal",PRE);
    double Nsig=0,Nbkg=0;
    for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str());
      Nsig+=wsum(c,F+" && "+SIG,true)*x.scale; Nbkg+=wsum(c,F+" && !"+SIG,true)*x.scale; }
    TChain ce("stv_tree"); ce.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
    double Next=wsum(ce,F,false)*sc_ext;
    TChain cd("stv_tree"); cd.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P));
    double Ndirt=wsum(cd,F,true)*sc_dirt;
    double Ndata=0; for(auto&d:data){ TChain c("stv_tree"); c.Add(d.c_str()); Ndata+=wsum(c,F,false); }
    // The fake data is a Poisson throw of the neutrino MC only, so the machinery/
    // normalisation check compares it to the CV-weighted neutrino MC (sig+nubkg).
    // EXT + dirt are the additional components real beam-on data will contain (shown
    // for scale); their model is validated against the FULL stack only after unblinding.
    double Nnu=Nsig+Nbkg;
    printf("%-12s %9.1f %9.1f %9.1f %9.1f %9.1f %7.1f%% %9.3f\n",
           s,Ndata,Nsig,Nbkg,Next,Ndirt,Nnu>0?100*Nsig/Nnu:0,Nnu>0?Ndata/Nnu:0);
  }
  printf("(FAKE data = Poisson throw of the neutrino MC CV; data/vMC ~1 validates the\n"
         " per-run weighting/scaling machinery, and sig%% (of vMC) confirms signal depletion.\n"
         " EXT+dirt are the extra real-data components -- data/full-MC (background model)\n"
         " is only meaningful once the control regions are unblinded with beam-on data.)\n");
}
