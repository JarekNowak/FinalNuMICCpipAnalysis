{ gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi.C+");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numu.root\",\"pred_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi(\"neutvect_numubar.root\",\"pred_numubar.root\")"); }
