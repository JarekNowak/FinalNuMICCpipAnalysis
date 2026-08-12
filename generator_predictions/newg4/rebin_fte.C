// rebin_fte.C — bring the generator FTE predictions onto the CURRENT (coarsened)
// analysis binning so every generator overlays on every observable. The analysis
// binning is a clean merge/subset of the generator binning, and the FTE bin content
// is the per-bin cross section (dsigma/dx * width), so merging = summing contents
// (errors in quadrature), dropping = omitting. Applied conditionally on bin count so
// an already-coarse generator (e.g. gibuu ppi=5) is left untouched.
//   ppi     : 6 -> 5  (drop the generator's first bin [0.113,0.175], below threshold)
//   costhpi : 5 -> 4  (merge old bins 3+4 -> [0.35,0.75])
//   thmupi  : 7 -> 5  (merge old 3+4 and 5+6)
//   pmu,costhmu : unchanged.  Reads all histos into memory, then recreates the file.
//   usage: root -l -b -q 'rebin_fte.C("nuwro","newg4")'
#include <vector>
static TH1D* remap(TH1D* h, const char* name, std::vector<std::vector<int>> groups){
  int nn=groups.size();
  TH1D* o=new TH1D(name,"",nn,0,nn); o->SetDirectory(0);
  for(int j=0;j<nn;j++){ double c=0,e2=0; for(int ob:groups[j]){ c+=h->GetBinContent(ob); e2+=h->GetBinError(ob)*h->GetBinError(ob);}
    o->SetBinContent(j+1,c); o->SetBinError(j+1,std::sqrt(e2)); }
  return o;
}
void rebin_fte(const char* gen, const char* gtag){
  const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  TString fn=Form("%s_%s_fte.root",gen,gtag);
  // 1) read everything into memory
  std::vector<TH1D*> keep;
  { TFile fi(fn,"read");
    for(auto o:obs){
      TH1D* h=(TH1D*)fi.Get(Form("%s_fte",o));
      if(!h){ printf("  %-6s %-8s : missing\n",gen,o); continue; }
      int nb=h->GetNbinsX();
      std::vector<std::vector<int>> g;
      if(!strcmp(o,"ppi") && nb==6)          g={{2},{3},{4},{5},{6}};
      else if(!strcmp(o,"costhpi") && nb==5) g={{1},{2},{3,4},{5}};
      else if(!strcmp(o,"thmupi") && nb==7)  g={{1},{2},{3,4},{5,6},{7}};
      if(g.empty()){ TH1D* c=(TH1D*)h->Clone(Form("%s_fte",o)); c->SetDirectory(0); keep.push_back(c); printf("  %-6s %-8s : nb=%d unchanged\n",gen,o,nb); }
      else { TH1D* nh=remap(h,Form("%s_fte",o),g); keep.push_back(nh); printf("  %-6s %-8s : %d -> %d\n",gen,o,nb,nh->GetNbinsX()); }
    }
  }
  // 2) recreate the file with the corrected histos
  TFile fo(fn,"recreate");
  for(auto h:keep){ h->SetDirectory(&fo); h->Write(); }
  fo.Close();
}
