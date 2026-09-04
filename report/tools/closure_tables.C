// closure_tables.C -- unf./truth integral ratio and the truth chi2 (from the unfolder log)
// for all 51 release extractions -> TSV
#include <fstream>
#include <regex>
void closure_tables(const char* out="closure_summary.tsv"){
  const char* PROC="/data/uboone/processed/"; const char* LOG="/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/";
  const char* cfgs[3]={"FHC5","RHCFULL","COMB"}; const char* lcfg[3]={"fhc5","rhcfull","comb"};
  const char* incl[6]={"pmu","ppi2bin","costhmu","costhpi","thmupi","thetamu"}; const char* lincl[6]={"pmu","ppi","costhmu","costhpi","thmupi","thetamu"};
  const char* p1p[11]={"Wpipr","Whad","dpt2bin","dalphat2bin","dphit2bin","pn2bin","pmu","ppi2bin","costhmu","costhpi","thmupi"};
  FILE* fo=fopen(out,"w"); fprintf(fo,"family\tconfig\tobs\tsigma_int\ttruth_int\tunf_over_truth\tchi2\tndf\tpval\n");
  auto one=[&](const char* fam,int c,const char* o,const char* lo){
    TString fn=Form("%sclosure_hists_xsec_%s%s_%s.root",PROC,TString(fam)=="1p"?"ccpi1p_":"",cfgs[c],o);
    TFile f(fn); if(f.IsZombie()){printf("missing %s\n",fn.Data());return;}
    TH1D* u=(TH1D*)f.Get("h_unfolded_nuwro"); TH1D* t=(TH1D*)f.Get("h_fakedata_truth");
    double su=0,st=0; for(int b=1;b<=u->GetNbinsX();b++){double w=u->GetBinWidth(b); su+=u->GetBinContent(b)*w; st+=t->GetBinContent(b)*w;}
    TString ln=Form("%s%s%s_%s.raw",LOG,TString(fam)=="1p"?"1p_":"",lcfg[c],lo);
    std::ifstream in(ln.Data()); std::string line; double chi2=-1,pv=-1; int ndf=-1; std::regex re("^truth: .* = ([0-9.e+-]+)/([0-9]+) bins?, p-value = ([0-9.e+-]+)");
    while(std::getline(in,line)){ std::smatch m; if(std::regex_search(line,m,re)){chi2=atof(m[1].str().c_str()); ndf=atoi(m[2].str().c_str()); pv=atof(m[3].str().c_str()); break;} }
    fprintf(fo,"%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.3f\t%d\t%.3f\n",fam,cfgs[c],o,su,st,su/st,chi2,ndf,pv);
  };
  for(int c=0;c<3;c++){ for(int i=0;i<6;i++) one("incl",c,incl[i],lincl[i]); for(int i=0;i<11;i++) one("1p",c,p1p[i],p1p[i]); }
  fclose(fo); printf("wrote %s\n",out);
}
