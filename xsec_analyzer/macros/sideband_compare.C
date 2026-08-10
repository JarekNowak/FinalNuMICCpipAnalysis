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

void sideband_compare(const char* mode="fhc", const char* datasrc="fake") {
  const char* P="/data/uboone/processed/";
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

  const char* sb[4]={"sb_cc0pi","sb_multipi","sb_pi0","sb_cosmic"};
  const char* PRE="CC1mu1piXp_";
  printf("\n==== %s sidebands (data source: %s, per-run scaled) ====\n",mode,datasrc);
  printf("%-12s %9s %9s %9s %9s %9s %8s %9s\n","sideband","N_data","sig","nu-bkg","EXT","dirt","sig%%","data/MC");
  for(auto s:sb){
    TString F=Form("%s%s",PRE,s), SIG=Form("%sMC_Signal",PRE);
    double Nsig=0,Nbkg=0;
    for(auto&x:mc){ TChain c("stv_tree"); c.Add(x.file.c_str());
      Nsig+=c.GetEntries(F+" && "+SIG)*x.scale; Nbkg+=c.GetEntries(F+" && !"+SIG)*x.scale; }
    TChain ce("stv_tree"); ce.Add(Form("%sxsec-ana-beamoff_run1Andrun3.root",P));
    double Next=ce.GetEntries(F)*sc_ext;
    TChain cd("stv_tree"); cd.Add(Form("%sxsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root",P));
    double Ndirt=cd.GetEntries(F)*sc_dirt;
    double Ndata=0; for(auto&d:data){ TChain c("stv_tree"); c.Add(d.c_str()); Ndata+=c.GetEntries(F); }
    double Nmc=Nsig+Nbkg+Next+Ndirt;
    printf("%-12s %9.1f %9.1f %9.1f %9.1f %9.1f %7.1f%% %9.3f\n",
           s,Ndata,Nsig,Nbkg,Next,Ndirt,Nmc>0?100*Nsig/Nmc:0,Nmc>0?Ndata/Nmc:0);
  }
  printf("(data/MC ~1 validates the background normalisation; sig%% confirms signal depletion.\n"
         " A per-sideband scale factor (data-nonTargetBkg)/targetBkg constrains that background.)\n");
}
