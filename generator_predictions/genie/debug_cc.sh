set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
root -l -b -q -e '
TFile fs("splines_numu_ar.root");
TDirectory* d=(TDirectory*)fs.Get("nu_mu_Ar40");
TIter it(d->GetListOfKeys()); TKey* k;
double qel=0,res=0,dis=0,mec=0,coh=0,other=0; int nq=0,nr=0,nd=0,nm=0,nc=0;
double E=1.0;
while((k=(TKey*)it())){ TString n=k->GetName();
  if(!n.Contains("cc")||n.Contains("nc")) continue;
  TGraph* g=(TGraph*)d->Get(n); double v=g->Eval(E); if(v<0)v=0;
  if(n.BeginsWith("qel")){qel+=v;nq++;}
  else if(n.BeginsWith("res")){res+=v;nr++;}
  else if(n.BeginsWith("dis")){dis+=v;nd++;}
  else if(n.BeginsWith("mec")){mec+=v;nm++;}
  else if(n.BeginsWith("coh")){coh+=v;nc++;}
  else{other+=v; printf("  OTHER cc graph: %s = %.3e\n",n.Data(),v);}
}
printf("At E=1GeV [1e-38 cm2/Ar]:\n");
printf("  QE (%d graphs) = %.3e\n",nq,qel);
printf("  RES(%d graphs) = %.3e\n",nr,res);
printf("  DIS(%d graphs) = %.3e\n",nd,dis);
printf("  MEC(%d graphs) = %.3e\n",nm,mec);
printf("  COH(%d graphs) = %.3e\n",nc,coh);
printf("  sum = %.3e\n", qel+res+dis+mec+coh+other);
' 2>&1 | grep -iE "GeV|QE|RES|DIS|MEC|COH|sum|OTHER"
