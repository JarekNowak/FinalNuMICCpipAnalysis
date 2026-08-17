// Apply the CC1mu1piXp signal to a GENIE gst tree, compute per-nucleon
// dsigma/dx. sigma_tot_cc = flux-averaged total CC per Ar (passed in); the
// per-nucleon dsdx = (sigma_tot_cc/40) * (N_bin/N) / width.  (GENIE tot_cc is
// per Ar; /40 to match the NuWro per-nucleon convention for later combine.)
void analyze_gst_2bin(const char* gstfile, double sigma_tot_cc_perAr, const char* outfile){
  TFile f(gstfile); TTree* t=(TTree*)f.Get("gst");
  int nf; int pdgf[400]; double pxf[400],pyf[400],pzf[400];
  bool cc;   // cc is Bool_t in the gst
  int nfpip,nfpim,nfpi0,nfkp,nfkm,nfk0;
  double pxl,pyl,pzl;   // primary lepton (muon for CC numu)
  t->SetBranchAddress("nf",&nf); t->SetBranchAddress("cc",&cc);
  t->SetBranchAddress("pdgf",pdgf); t->SetBranchAddress("pxf",pxf);
  t->SetBranchAddress("pyf",pyf); t->SetBranchAddress("pzf",pzf);
  t->SetBranchAddress("pxl",&pxl); t->SetBranchAddress("pyl",&pyl);
  t->SetBranchAddress("pzl",&pzl);
  t->SetBranchAddress("nfpip",&nfpip); t->SetBranchAddress("nfpim",&nfpim);
  t->SetBranchAddress("nfpi0",&nfpi0); t->SetBranchAddress("nfkp",&nfkp);
  t->SetBranchAddress("nfkm",&nfkm); t->SetBranchAddress("nfk0",&nfk0);
  auto is_meson=[](int pdg){int a=abs(pdg); if(a>=9900000)return false;
    if((a/1000)%10!=0)return false; if((a/100)%10==0)return false;
    if(a>=901&&a<=930)return false; if(a==110||a==990||a==998||a==999||a==100)return false; return true;};
  auto is_kaon=[](int pdg){int a=abs(pdg); return a==321||a==311||pdg==310||pdg==130;};
  std::vector<double> EPMU={0.150,0.350,0.550,0.750,0.950,1.250,1.750,3.000};
  std::vector<double> EPPI={0.175,0.205,99.0};
  std::vector<double> ECMU={-1.0,0.45,0.65,0.80,0.90,1.0};
  std::vector<double> ECPI={-1.0,-0.10,0.35,0.55,0.75,1.0};
  std::vector<double> ETH={0.0,0.60,0.85,1.10,1.30,1.52,1.85,2.60};
  TH1D* hpmu=new TH1D("pmu","",EPMU.size()-1,EPMU.data());
  TH1D* hppi=new TH1D("ppi","",EPPI.size()-1,EPPI.data());
  TH1D* hcmu=new TH1D("costhmu","",ECMU.size()-1,ECMU.data());
  TH1D* hcpi=new TH1D("costhpi","",ECPI.size()-1,ECPI.data());
  TH1D* hth=new TH1D("thmupi","",ETH.size()-1,ETH.data());
  long N=t->GetEntries(); long sig=0;
  for(long i=0;i<N;i++){ t->GetEntry(i); if(!cc) continue;
    // Topology from the count branches (robust). Muon = primary lepton
    // (exactly 1 for CC numu). Charged pion momentum from the hadron arrays.
    if(nfpip+nfpim!=1||nfpi0!=0||nfkp+nfkm+nfk0!=0) continue;
    TVector3 pmu(pxl,pyl,pzl);
    int npi=0,nheavy=0; TVector3 ppi;
    for(int j=0;j<nf;j++){ int pdg=pdgf[j],a=abs(pdg); TVector3 p(pxf[j],pyf[j],pzf[j]);
      if(a==211){npi++;ppi=p;}
      else if(a!=111&&a!=211&&!is_kaon(pdg)&&is_meson(pdg))nheavy++; }
    if(npi!=1||nheavy!=0) continue;
    double Pmu=pmu.Mag(),Ppi=ppi.Mag(),th=pmu.Angle(ppi);
    if(Pmu<=0.15||Ppi<=0.175||th>=2.6) continue;      // ppi relaxed to 0.113 (low bin)
    sig++; hppi->Fill(Ppi);                            // ppi extends to 0.113
    if(Ppi>0.175){ hpmu->Fill(Pmu);hcmu->Fill(pmu.CosTheta());hcpi->Fill(ppi.CosTheta());hth->Fill(th); } // others: standard 0.175 phase space
  }
  double sperNuc=sigma_tot_cc_perAr/40.0;
  auto norm=[&](TH1D* h){ for(int b=1;b<=h->GetNbinsX();++b){ double c=h->GetBinContent(b),w=h->GetBinWidth(b);
    h->SetBinContent(b,sperNuc*(c/(double)N)/w); h->SetBinError(b,sperNuc*(sqrt(c)/(double)N)/w);} };
  norm(hpmu);norm(hppi);norm(hcmu);norm(hcpi);norm(hth);
  TFile fo(outfile,"recreate"); hpmu->Write();hppi->Write();hcmu->Write();hcpi->Write();hth->Write(); fo.Close();
  printf("GENIE gst: N=%ld CC1pi-signal=%ld (%.2f%% of all, %.2f%% of CC)\n",N,sig,100.*sig/N,100.*sig/N);
  printf("  sigma_tot_cc(perAr)=%.3e  -> per-nucleon CC1pi dsdx written\n",sigma_tot_cc_perAr);
}
