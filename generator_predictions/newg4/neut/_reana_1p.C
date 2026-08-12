{ gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi1p.C+");
  gROOT->ProcessLine("neut_cc1pi1p(\"neutvect_numu.root\",\"neut_1p_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi1p(\"neutvect_numubar.root\",\"neut_1p_numubar.root\")"); }
