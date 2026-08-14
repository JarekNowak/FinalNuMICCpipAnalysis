// Train the multi-pion pion-ID BDT (pion vs muon/proton/other) with TMVA.
// Produces mp_pion_bdt/weights/TMVAClassification_BDT.weights.xml for the selection.
void train_bdt(){
  TMVA::Tools::Instance();
  TFile*in=TFile::Open("../booster_decision_tree/mp_pion_bdt/bdt_train.root");
  TTree*tr=(TTree*)in->Get("cand");
  TFile*out=TFile::Open("../booster_decision_tree/mp_pion_bdt/tmva_train_out.root","RECREATE");
  TMVA::Factory f("TMVAClassification",out,
    "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
  TMVA::DataLoader* d = new TMVA::DataLoader("../booster_decision_tree/mp_pion_bdt");
  for(auto v:{"llr","bragg_p","bragg_mu","bragg_mip","bragg_pion","trk_score","length","dist"})
    d->AddVariable(v,'F');
  TCut sig="is_pion==1", bkg="is_pion==0";
  d->SetInputTrees(tr,sig,bkg);
  d->PrepareTrainingAndTestTree(sig,bkg,
    "SplitMode=Random:NormMode=EqualNumEvents:!V");
  f.BookMethod(d,TMVA::Types::kBDT,"BDT",
    "!H:!V:NTrees=600:MinNodeSize=2.5%:MaxDepth=4:BoostType=Grad:Shrinkage=0.10"
    ":UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40");
  f.TrainAllMethods(); f.TestAllMethods(); f.EvaluateAllMethods();
  out->Close();
  printf("TRAIN_DONE (weights in ../booster_decision_tree/mp_pion_bdt/weights/)\n");
}
