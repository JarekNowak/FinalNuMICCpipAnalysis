// analyze_gst_1p.C — proton-tagged (CC1mu1pi1p) W/TKI generator prediction from a
// GENIE gst tree. Same inclusive signal as analyze_gst.C plus a leading true proton
// above 0.3 GeV/c; computes the six W/TKI observables (wtki_gen.h) at truth level and
// writes per-nucleon dsigma/dx, later combined+normalised like the inclusive prediction.
//   root -l -b -q 'genie/analyze_gst_1p.C("newg4/g_numu.gst.root",13.90129,"newg4/genie_1p_numu.root")'
#include "../wtki_gen.h"
void analyze_gst_1p(const char* gstfile, double sigma_tot_cc_perAr, const char* outfile){
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
  auto is_meson=[](int pdg){int a=abs(pdg); if(a>=9900000)return false;
    if((a/1000)%10!=0)return false; if((a/100)%10==0)return false;
    if(a>=901&&a<=930)return false; if(a==110||a==990||a==998||a==999||a==100)return false; return true;};
  auto is_kaon=[](int pdg){int a=abs(pdg); return a==321||a==311||pdg==310||pdg==130;};
  const char* OBS[6]={"Wpipr","Whad","dpt","dalphat","dphit","pn"};
  TH1D* h[6]; for(int o=0;o<6;o++){ auto e=wtki::edges(OBS[o]); h[o]=new TH1D(OBS[o],"",e.size()-1,e.data()); }
  long N=t->GetEntries(); long sig=0;
  for(long i=0;i<N;i++){ t->GetEntry(i); if(!cc) continue;
    if(nfpip+nfpim!=1||nfpi0!=0||nfkp+nfkm+nfk0!=0) continue;
    TVector3 pmu(pxl,pyl,pzl);
    int npi=0,nheavy=0; TVector3 ppi,ppr; double leadp=0;
    for(int j=0;j<nf;j++){ int pdg=pdgf[j],a=abs(pdg); TVector3 p(pxf[j],pyf[j],pzf[j]);
      if(a==211){npi++;ppi=p;}
      else if(pdg==2212){ if(p.Mag()>leadp){leadp=p.Mag(); ppr=p;} }
      else if(a!=111&&a!=211&&!is_kaon(pdg)&&is_meson(pdg))nheavy++; }
    if(npi!=1||nheavy!=0) continue;
    if(leadp<=0.3) continue;                                  // require leading proton > 0.3 GeV/c
    double Pmu=pmu.Mag(),Ppi=ppi.Mag(),th=pmu.Angle(ppi);
    if(Pmu<=0.15||Ppi<=0.175||th>=2.6) continue;
    sig++;
    wtki::Obs w=wtki::compute(pmu,ppi,ppr);
    double v[6]={w.Wpipr,w.Whad,w.dpt,w.dalphat,w.dphit,w.pn};
    for(int o=0;o<6;o++){ auto e=wtki::edges(OBS[o]); h[o]->Fill(wtki::clamp(v[o],e)); }
  }
  double sperNuc=sigma_tot_cc_perAr/40.0;
  auto norm=[&](TH1D* hh){ for(int b=1;b<=hh->GetNbinsX();++b){ double c=hh->GetBinContent(b),w=hh->GetBinWidth(b);
    hh->SetBinContent(b,sperNuc*(c/(double)N)/w); hh->SetBinError(b,sperNuc*(sqrt(c)/(double)N)/w);} };
  for(int o=0;o<6;o++) norm(h[o]);
  TFile fo(outfile,"recreate"); for(int o=0;o<6;o++) h[o]->Write(); fo.Close();
  printf("GENIE gst 1p: N=%ld  proton-tagged W/TKI signal=%ld (%.3f%%)  wrote %s\n",N,sig,100.*sig/N,outfile);
}
