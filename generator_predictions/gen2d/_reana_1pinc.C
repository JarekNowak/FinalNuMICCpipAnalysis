// _reana_1pinc.C -- NEUT part of regen_1pinc.sh (run by reana_neut_1pinc.sh inside the SL7 container, ROOT 5.34):
// the proton-tagged reader with the muon/pion observables at the proton-tagged binnings (obs_1pinc.h).
{ gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_1p_inc.C+");
  gROOT->ProcessLine("neut_1p_inc(\"../newg4/neut/neutvect_numu.root\",\"neut_1p_inc_numu.root\")");
  gROOT->ProcessLine("neut_1p_inc(\"../newg4/neut/neutvect_numubar.root\",\"neut_1p_inc_numubar.root\")"); }
