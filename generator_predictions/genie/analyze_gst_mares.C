// analyze_gst_mares.C — like analyze_gst.C, but also reads a grwght1p MaCCRES
// weight file and produces the M_A^RES nominal / -1sigma / +1sigma CC1pi
// predictions. The reweighted differential xsec = nominal * <w>_signal per bin
// (the flux-averaged total sigma cancels between numerator and normalisation),
// so we fill with the per-event MaCCRES weight and normalise identically to the
// nominal (analyze_gst) with the nominal sigma_tot. Writes <prefix>_{nom,minus,
// plus}.root, each with the standard histogram names for combine_newg4.C.
//
//   root -l -b -q 'analyze_gst_mares.C("g.gst.root","rw.root",SIGTOT,"pred")'
#include "TArrayF.h"
void analyze_gst_mares(const char* gstfile, const char* wfile,
    double sigma_tot_cc_perAr, const char* prefix){
  TFile f(gstfile); TTree* t=(TTree*)f.Get("gst");
  int nf; int pdgf[400]; double pxf[400],pyf[400],pzf[400];
  bool cc; int nfpip,nfpim,nfpi0,nfkp,nfkm,nfk0; double pxl,pyl,pzl;
  t->SetBranchAddress("nf",&nf); t->SetBranchAddress("cc",&cc);
  t->SetBranchAddress("pdgf",pdgf); t->SetBranchAddress("pxf",pxf);
  t->SetBranchAddress("pyf",pyf); t->SetBranchAddress("pzf",pzf);
  t->SetBranchAddress("pxl",&pxl); t->SetBranchAddress("pyl",&pyl); t->SetBranchAddress("pzl",&pzl);
  t->SetBranchAddress("nfpip",&nfpip); t->SetBranchAddress("nfpim",&nfpim);
  t->SetBranchAddress("nfpi0",&nfpi0); t->SetBranchAddress("nfkp",&nfkp);
  t->SetBranchAddress("nfkm",&nfkm); t->SetBranchAddress("nfk0",&nfk0);
  TFile fw(wfile); TTree* tw=(TTree*)fw.Get("MaCCRES");
  TArrayF* weights=0; tw->SetBranchAddress("weights",&weights);
  TArrayF* twk=0;     tw->SetBranchAddress("twkdials",&twk);
  auto is_meson=[](int pdg){int a=abs(pdg); if(a>=9900000)return false; if((a/1000)%10!=0)return false;
    if((a/100)%10==0)return false; if(a>=901&&a<=930)return false;
    if(a==110||a==990||a==998||a==999||a==100)return false; return true;};
  auto is_kaon=[](int pdg){int a=abs(pdg); return a==321||a==311||pdg==310||pdg==130;};
  std::vector<double> EPMU={0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.30,2.60,3.00};
  std::vector<double> EPPI={0.175,0.20,0.30,0.40,0.50,1.00};
  std::vector<double> ECMU={-1.0,-0.5,0.0,0.2,0.4,0.55,0.65,0.75,0.82,0.88,0.93,0.97,1.0};
  std::vector<double> ECPI={-1.0,-0.7,-0.4,-0.2,0.0,0.2,0.35,0.5,0.65,0.78,0.88,0.95,1.0};
  std::vector<double> ETH={0.0,0.314,0.628,0.942,1.257,1.571,1.885,2.199,2.513,2.6};
  const char* vn[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  std::vector< std::vector<double> > edges={EPMU,EPPI,ECMU,ECPI,ETH};
  // [variation 0=nom,1=minus,2=plus][observable]
  TH1D* h[3][5];
  const char* suf[3]={"nom","minus","plus"};
  for(int v=0;v<3;v++) for(int o=0;o<5;o++)
    h[v][o]=new TH1D(Form("%s_%s",vn[o],suf[v]),"",edges[o].size()-1,edges[o].data());
  long N=t->GetEntries(), Nw=tw->GetEntries(); long sig=0;
  printf("gst N=%ld  weight N=%ld\n", N, Nw);
  int iw_m=-1, iw_p=-1; // indices of -1 and +1 sigma in the tweak array
  for(long i=0;i<N;i++){ t->GetEntry(i); tw->GetEntry(i);
    if(iw_m<0 && twk){ for(int k=0;k<twk->GetSize();k++){ if(fabs(twk->At(k)+1)<1e-3)iw_m=k; if(fabs(twk->At(k)-1)<1e-3)iw_p=k; } }
    if(!cc) continue;
    if(nfpip+nfpim!=1||nfpi0!=0||nfkp+nfkm+nfk0!=0) continue;
    TVector3 pmu(pxl,pyl,pzl); int npi=0,nheavy=0; TVector3 ppi;
    for(int j=0;j<nf;j++){ int pdg=pdgf[j],a=abs(pdg); TVector3 p(pxf[j],pyf[j],pzf[j]);
      if(a==211){npi++;ppi=p;} else if(a!=111&&a!=211&&!is_kaon(pdg)&&is_meson(pdg))nheavy++; }
    if(npi!=1||nheavy!=0) continue;
    double Pmu=pmu.Mag(),Ppi=ppi.Mag(),th=pmu.Angle(ppi);
    if(Pmu<=0.15||Ppi<=0.175||th>=2.6) continue;
    double wm = (iw_m>=0)? weights->At(iw_m):1.0;
    double wp = (iw_p>=0)? weights->At(iw_p):1.0;
    double val[5]={Pmu,Ppi,pmu.CosTheta(),ppi.CosTheta(),th};
    for(int o=0;o<5;o++){ h[0][o]->Fill(val[o],1.0); h[1][o]->Fill(val[o],wm); h[2][o]->Fill(val[o],wp); }
    sig++;
  }
  double sperNuc=sigma_tot_cc_perAr/40.0;
  auto norm=[&](TH1D* hh){ for(int b=1;b<=hh->GetNbinsX();++b){ double c=hh->GetBinContent(b),w=hh->GetBinWidth(b);
    hh->SetBinContent(b,sperNuc*(c/(double)N)/w); hh->SetBinError(b,sperNuc*(sqrt(fabs(c))/(double)N)/w);} };
  for(int v=0;v<3;v++){
    TFile fo(Form("%s_%s.root",prefix,suf[v]),"recreate");
    for(int o=0;o<5;o++){ norm(h[v][o]); h[v][o]->SetName(vn[o]); h[v][o]->Write(); }
    fo.Close();
  }
  printf("MaCCRES: CC1pi-signal=%ld  (twk idx -1sig=%d +1sig=%d)  wrote %s_{nom,minus,plus}.root\n", sig, iw_m, iw_p, prefix);
}
