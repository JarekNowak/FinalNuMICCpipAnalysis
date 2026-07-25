#include "WienerSVD.h"
void fw_inputs_indep(){
  // Framework costhmu numuMC inputs from the univmake
  TFile f("/data/uboone/processed/ccpi_Run1_costhmu_univmake.root");
  const char* B="ccpi_CC1mu1piXp_costhmu_1D";
  TDirectoryFile* mc=(TDirectoryFile*)f.Get(Form("%s/+data+uboone+processed+xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",B));
  TH2D* mig=(TH2D*)mc->Get("weight_TunedCentralValue_UBGenie_0_2d");   // X=true(13) Y=reco(12)
  TH1D* tru=(TH1D*)mc->Get("weight_TunedCentralValue_UBGenie_0_true"); // true(13)
  // Build h_smear_sel (12 true x 12 reco) and h_true_gen (12), self-closure data = signal reco
  TH2D* smear=new TH2D("smear","",12,0,12,12,0,12);
  TH1D* tgen=new TH1D("tgen","",12,0,12);
  TH1D* dreco=new TH1D("dreco","",12,0,12);
  TH1D* zbkg=new TH1D("zbkg","",12,0,12);
  for(int t=1;t<=12;t++){ tgen->SetBinContent(t,tru->GetBinContent(t));
    for(int r=1;r<=12;r++){ double m=mig->GetBinContent(t,r); smear->SetBinContent(t,r,m); dreco->AddBinContent(r,m); } }
  printf("framework inputs: smear int=%.1f  true_gen int=%.1f  data_reco int=%.1f\n", smear->Integral(), tgen->Integral(), dreco->Integral());
  double tsum=tgen->Integral();
  for(std::string ct : {std::string("identity"),std::string("smooth2")}){
    WienerSVDResult r=RunWienerSVD_FW(smear,tgen,dreco,zbkg,ct,false);
    int n=r.x_unfolded.GetNrows(); double x=0,ac=0; for(int i=0;i<n;i++){x+=r.x_unfolded[i];ac+=r.Ac[i][i];}
    printf("INDEP-on-FW-inputs  %-9s  xhat/truth=%.4f  Ac_diag=%.3f\n",ct.c_str(),x/tsum,ac/n);
  }
}
