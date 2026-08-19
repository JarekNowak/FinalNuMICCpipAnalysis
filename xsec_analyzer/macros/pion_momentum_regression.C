// pion_momentum_regression.C -- can a multivariate regression beat range for p_pi?
//
// Every single-quantity estimator tested saturates, because each measures what the pion did
// BEFORE it interacted: range 0.83->0.32 across the momentum range, calorimetry 1.14->0.38
// (and 3-5x noisier), MCS unbiased at 1.01 in the top bin but with 210% RMS. They are biased
// DIFFERENTLY, so a regression may be able to combine them.
//
// Target is p_true; every input is reco-level, so the result is a legitimate estimator and
// not a truth-level cheat. The check that matters is NOT bias or RMS but the MIGRATION
// DIAGONAL: a regression can flatter its own residuals while sculpting the response, which
// is the trap documented for the tight golden-pion cut.
//
//   usage: root -l -b -q 'macros/pion_momentum_regression.C()'

#include <vector>
#include <cmath>

void pion_momentum_regression() {

  const double mmu=0.10566, mpi=0.13957;
  const char* files[4]={
    "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
    "/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root"};

  TString out="/data/uboone/processed/validate/pion_regr_sample.root";
  TFile fo(out,"recreate");
  TTree* tr=new TTree("reg","pion momentum regression sample");
  float p_true,p_range,p_mcs,calo,len,wstd,wmean,wsep;
  float bragg_pion,bragg_mip,bragg_mu,bragg_p,tscore,ndau;
  int run_id, golden;
  tr->Branch("p_true",&p_true); tr->Branch("p_range",&p_range); tr->Branch("p_mcs",&p_mcs);
  tr->Branch("calo",&calo);     tr->Branch("len",&len);        tr->Branch("wstd",&wstd);
  tr->Branch("wmean",&wmean);   tr->Branch("wsep",&wsep);
  tr->Branch("bragg_pion",&bragg_pion); tr->Branch("bragg_mip",&bragg_mip);
  tr->Branch("bragg_mu",&bragg_mu);     tr->Branch("bragg_p",&bragg_p);
  tr->Branch("tscore",&tscore); tr->Branch("ndau",&ndau);
  tr->Branch("run_id",&run_id); tr->Branch("golden",&golden);

  long nsel=0;
  for(int fi=0;fi<4;++fi){
    TFile* f=TFile::Open(files[fi]); if(!f||f->IsZombie()){printf("  [skip] %s\n",files[fi]);continue;}
    TTree* t=(TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    t->SetBranchStatus("*",0);
    for(auto b:{"backtracked_pdg","backtracked_px","backtracked_py","backtracked_pz",
                "backtracked_purity","mc_pdg","mc_px","mc_py","mc_pz","mc_end_p",
                "trk_score_v","trk_len_v","trk_range_muon_mom_v","trk_mcs_muon_mom_v",
                "trk_calo_energy_y_v","trk_avg_deflection_stdev_v","trk_avg_deflection_mean_v",
                "trk_avg_deflection_separation_mean_v","trk_bragg_pion_v","trk_bragg_mip_v",
                "trk_bragg_mu_v","trk_bragg_p_v","pfp_trk_daughters_v","pfp_shr_daughters_v"})
      t->SetBranchStatus(b,1);
    std::vector<int>*bpdg=0,*mpdg=0;
    std::vector<float> *bpx=0,*bpy=0,*bpz=0,*bpur=0,*mpx=0,*mpy=0,*mpz=0,*mep=0;
    std::vector<float> *ts=0,*tl=0,*rm=0,*mcs=0,*ce=0,*ws=0,*wm=0,*wsp=0,*bpi=0,*bmip=0,*bmu=0,*bp=0;
    std::vector<unsigned int> *dtr=0,*dsh=0;
    t->SetBranchAddress("backtracked_pdg",&bpdg);t->SetBranchAddress("backtracked_px",&bpx);
    t->SetBranchAddress("backtracked_py",&bpy);  t->SetBranchAddress("backtracked_pz",&bpz);
    t->SetBranchAddress("backtracked_purity",&bpur);
    t->SetBranchAddress("mc_pdg",&mpdg);t->SetBranchAddress("mc_px",&mpx);
    t->SetBranchAddress("mc_py",&mpy);  t->SetBranchAddress("mc_pz",&mpz);
    t->SetBranchAddress("mc_end_p",&mep);
    t->SetBranchAddress("trk_score_v",&ts); t->SetBranchAddress("trk_len_v",&tl);
    t->SetBranchAddress("trk_range_muon_mom_v",&rm);
    t->SetBranchAddress("trk_mcs_muon_mom_v",&mcs);
    t->SetBranchAddress("trk_calo_energy_y_v",&ce);
    t->SetBranchAddress("trk_avg_deflection_stdev_v",&ws);
    t->SetBranchAddress("trk_avg_deflection_mean_v",&wm);
    t->SetBranchAddress("trk_avg_deflection_separation_mean_v",&wsp);
    t->SetBranchAddress("trk_bragg_pion_v",&bpi);t->SetBranchAddress("trk_bragg_mip_v",&bmip);
    t->SetBranchAddress("trk_bragg_mu_v",&bmu);  t->SetBranchAddress("trk_bragg_p_v",&bp);
    t->SetBranchAddress("pfp_trk_daughters_v",&dtr);
    t->SetBranchAddress("pfp_shr_daughters_v",&dsh);

    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
      if(!bpdg||!mpdg)continue;
      for(size_t j=0;j<bpdg->size();++j){
        if(abs(bpdg->at(j))!=211)continue;
        if(bpur->at(j)<0.5||ts->at(j)<0.5||tl->at(j)<5.)continue;
        double pr=rm->at(j); if(pr<=0.)continue;
        int best=-1;
        for(size_t k=0;k<mpdg->size();++k){ if(mpdg->at(k)!=bpdg->at(j))continue;
          if(fabs(mpx->at(k)-bpx->at(j))+fabs(mpy->at(k)-bpy->at(j))+fabs(mpz->at(k)-bpz->at(j))<1e-4){best=k;break;} }
        if(best<0)continue;
        double pt=sqrt(pow(mpx->at(best),2)+pow(mpy->at(best),2)+pow(mpz->at(best),2));
        if(pt<0.175)continue;
        double E=sqrt(pr*pr+mmu*mmu)-mmu+mpi;
        p_range=sqrt(std::max(0.,E*E-mpi*mpi));
        p_mcs = (mcs&&mcs->at(j)>0.&&mcs->at(j)<5.) ? mcs->at(j) : -1.f;
        double raw = ce? ce->at(j) : -1.;
        calo = (raw>0.&&raw<2000.) ? raw/1000. : -1.f;   // GeV; the branch has outliers
        len=tl->at(j); wstd=ws->at(j); wmean=wm->at(j); wsep=wsp->at(j);
        auto san=[](float v){return std::isfinite(v)?v:-1.f;};
        bragg_pion=san(bpi->at(j)); bragg_mip=san(bmip->at(j));
        bragg_mu=san(bmu->at(j));   bragg_p=san(bp->at(j));
        tscore=ts->at(j); ndau=dtr->at(j)+dsh->at(j);
        p_true=pt; golden=(mep->at(best)<0.01); run_id=fi;
        tr->Fill(); ++nsel;
      } }
    f->Close();
  }
  fo.cd(); tr->Write();
  printf("\n  regression sample: %ld pion tracks\n", nsel);

  TMVA::Tools::Instance();
  TFile* fout=TFile::Open("/data/uboone/processed/validate/pion_regr_tmva.root","RECREATE");
  TMVA::Factory factory("PionMom", fout,
    "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
  TMVA::DataLoader* d=new TMVA::DataLoader("pion_regr");
  for(auto v:{"p_range","p_mcs","calo","len","wstd","wmean","wsep",
              "bragg_pion","bragg_mip","bragg_mu","bragg_p","tscore","ndau"})
    d->AddVariable(v,'F');
  d->AddTarget("p_true");
  d->AddRegressionTree(tr,1.0);
  // Train on runs 0 and 2 only. NOTE the regression form of this call takes (cut, options);
  // the classification form (trainCut, testCut, options) silently misparses and then dies in
  // TMVA with an out_of_range on the target index. The held-out run is evaluated separately
  // with a TMVA::Reader, which also keeps the split by RUN rather than by event so the
  // quoted performance is not inflated by run-correlated detector conditions.
  d->PrepareTrainingAndTestTree(TCut("run_id==0||run_id==2"),
    "SplitMode=Random:NormMode=NumEvents:!V");
  factory.BookMethod(d,TMVA::Types::kBDT,"BDTG",
    "!H:!V:NTrees=600:MinNodeSize=1%:MaxDepth=4:BoostType=Grad:Shrinkage=0.1:"
    "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40");
  factory.TrainAllMethods(); factory.TestAllMethods(); factory.EvaluateAllMethods();
  fout->Close(); fo.Close();
  printf("\n  trained; evaluate with macros/pion_regression_eval.C\n");
}
