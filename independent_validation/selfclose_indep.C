#include "WienerSVD.h"
void selfclose_indep(){
  TFile f("unfolding_inputs.root");
  const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  printf("%-9s %-16s %11s %9s\n","obs","method","xhat/truth","Ac_diag");
  for(int o=0;o<5;o++){
    TH2D* smear=(TH2D*)f.Get(Form("h_smear_sel_%s",obs[o]));
    TH1D* truth=(TH1D*)f.Get(Form("h_true_gen_%s",obs[o]));
    TH1D* recosig=(TH1D*)f.Get(Form("h_reco_sig_%s",obs[o]));
    if(!smear||!truth||!recosig){printf("%s missing\n",obs[o]);continue;}
    TH1D* zbkg=(TH1D*)recosig->Clone(Form("zb_%s",obs[o])); zbkg->Reset();
    double tru=0; for(int b=1;b<=truth->GetNbinsX();b++) tru+=truth->GetBinContent(b);
    for(std::string ct : {std::string("identity"),std::string("smooth2")}){
      WienerSVDResult rf=RunWienerSVD_FW(smear,truth,recosig,zbkg,ct,false);
      WienerSVDResult rd=RunWienerSVD(smear,truth,recosig,zbkg,ct,false);
      int n=rf.x_unfolded.GetNrows(); double xf=0,xd=0,af=0,ad=0;
      for(int i=0;i<n;i++){xf+=rf.x_unfolded[i];xd+=rd.x_unfolded[i];af+=rf.Ac[i][i];ad+=rd.Ac[i][i];}
      printf("%-9s FW(BNL)_%-8s %11.4f %9.3f\n",obs[o],ct.c_str(),xf/tru,af/n);
      printf("%-9s direct_%-9s %11.4f %9.3f\n",obs[o],ct.c_str(),xd/tru,ad/n);
    }
  }
}
