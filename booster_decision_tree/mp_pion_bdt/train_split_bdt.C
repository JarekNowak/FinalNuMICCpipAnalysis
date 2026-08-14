// Train two momentum-specialized pion-ID BDTs: mp_pion_bdt_soft (length<20cm, the
// hard-to-ID short-track regime) and mp_pion_bdt_hard (length>=20cm). Same features
// as the single BDT. Applied per-candidate by length in the selection.
void train_one(const char* name, const char* lengthcut){
  TFile*in=TFile::Open("../booster_decision_tree/mp_pion_bdt/bdt_train.root");
  TTree*tr=(TTree*)in->Get("cand");
  TFile*out=TFile::Open(Form("../booster_decision_tree/%s/tmva_out.root",name),"RECREATE");
  TMVA::Factory f("TMVAClassification",out,"!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification");
  TMVA::DataLoader* d=new TMVA::DataLoader(Form("../booster_decision_tree/%s",name));
  for(auto v:{"llr","bragg_p","bragg_mu","bragg_mip","bragg_pion","trk_score","length","dist"}) d->AddVariable(v,'F');
  TCut sig=Form("is_pion==1 && %s",lengthcut), bkg=Form("is_pion==0 && %s",lengthcut);
  d->SetInputTrees(tr,sig,bkg);
  d->PrepareTrainingAndTestTree(sig,bkg,"SplitMode=Random:NormMode=EqualNumEvents:!V");
  f.BookMethod(d,TMVA::Types::kBDT,"BDT",
    "!H:!V:NTrees=600:MinNodeSize=2.5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40");
  f.TrainAllMethods(); f.TestAllMethods(); f.EvaluateAllMethods();
  out->Close(); in->Close();
}
void train_split_bdt(){
  gSystem->mkdir("../booster_decision_tree/mp_pion_bdt_soft",true);
  gSystem->mkdir("../booster_decision_tree/mp_pion_bdt_hard",true);
  train_one("mp_pion_bdt_soft","length<20");
  train_one("mp_pion_bdt_hard","length>=20");
  printf("SPLIT_TRAIN_DONE\n");
}
