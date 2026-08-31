// export_curves.C -- release the plotted curves themselves, in measurement space.
//
// A_C is released too, but a user who wants to compare against the models shown in the
// note should not have to re-derive the chain from a raw generator file: the prediction
// passes through a units conversion on load, the A_C multiplication, and a
// 1/(conv_factor x width) transform, and reproducing that exactly from the release alone
// is not currently demonstrated. Releasing the already-smeared curves removes the need.
//
//   root -l -b -q 'macros/export_curves.C("../report/data_release")'
#include <sys/stat.h>

static void dump(const char* path,const char* tag,const char* outdir,FILE* idx){
  TFile* f=TFile::Open(path);
  if(!f||f->IsZombie()) return;
  const char* keys[7]={"h_unfolded_nuwro","h_fakedata_truth","h_genie_tune",
                       "h_gen_GENIE","h_gen_GiBUU","h_gen_NEUT","h_gen_NuWro"};
  const char* lbl[7] ={"unfolded_data","truth_smeared","tune_smeared",
                       "GENIE_smeared","GiBUU_smeared","NEUT_smeared","NuWro_smeared"};
  TH1* ref=nullptr;
  for(int k=0;k<7 && !ref;++k) ref=(TH1*)f->Get(keys[k]);
  if(!ref){ f->Close(); return; }
  int n=ref->GetNbinsX();
  TString out=TString::Format("%s/curves_%s.tsv",outdir,tag);
  FILE* fp=fopen(out.Data(),"w");
  fprintf(fp,"# Measurement-space curves for %s.\n",tag);
  fprintf(fp,"# All model columns are ALREADY smeared by A_C and are directly\n");
  fprintf(fp,"# comparable with unfolded_data bin by bin. Values are dsigma/dx in\n");
  fprintf(fp,"# 1e-38 cm^2/(unit)/Ar; multiply by width to integrate.\n");
  fprintf(fp,"bin\tlow\thigh\twidth\tunfolded_data\tstat_err");
  for(int k=1;k<7;++k) fprintf(fp,"\t%s",lbl[k]);
  fprintf(fp,"\n");
  TH1* hu=(TH1*)f->Get("h_unfolded_nuwro");
  for(int b=1;b<=n;++b){
    fprintf(fp,"%d\t%.6g\t%.6g\t%.6g",b-1,
      ref->GetBinLowEdge(b),ref->GetBinLowEdge(b)+ref->GetBinWidth(b),ref->GetBinWidth(b));
    fprintf(fp,"\t%.8g\t%.8g", hu?hu->GetBinContent(b):0., hu?hu->GetBinError(b):0.);
    for(int k=1;k<7;++k){
      TH1* h=(TH1*)f->Get(keys[k]);
      fprintf(fp,"\t%.8g", h?h->GetBinContent(b):std::nan(""));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  fprintf(idx,"%s\t%d\n",tag,n);
  printf("  wrote curves_%-28s %d bins\n",tag,n);
  f->Close();
}

void export_curves(const char* outdir="../report/data_release"){
  mkdir(outdir,0755);
  TString ip=TString::Format("%s/index_curves.tsv",outdir);
  FILE* idx=fopen(ip.Data(),"w");
  fprintf(idx,"# Measurement-space curves; model columns are already A_C-smeared.\ntag\tn_bins\n");
  const char* cfgs[3]={"FHC5","RHCFULL","COMB"};
  const char* incl[6]={"pmu","ppi2bin","costhmu","costhpi","thmupi","thetamu"};
  const char* p1p[11]={"pmu","ppi2bin","costhmu","costhpi","thmupi","Wpipr","Whad",
                       "dpt2bin","dphit2bin","dalphat2bin","pn2bin"};
  for(auto c:cfgs) for(auto o:incl){
    TString p=TString::Format("/data/uboone/processed/closure_hists_xsec_%s_%s.root",c,o);
    if(!gSystem->AccessPathName(p)) dump(p,TString::Format("incl_%s_%s",c,o).Data(),outdir,idx);
  }
  for(auto c:cfgs) for(auto o:p1p){
    TString p=TString::Format("/data/uboone/processed/closure_hists_xsec_ccpi1p_%s_%s.root",c,o);
    if(!gSystem->AccessPathName(p)) dump(p,TString::Format("1p_%s_%s",c,o).Data(),outdir,idx);
  }
  fclose(idx);
}
