{ gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_incl.C+");
  gROOT->ProcessLine("neut_incl(\"../newg4/neut/neutvect_numu.root\",\"neut_ext_numu.root\")");
  gROOT->ProcessLine("neut_incl(\"../newg4/neut/neutvect_numubar.root\",\"neut_ext_numubar.root\")");
  gROOT->ProcessLine(".L neut_1p.C+");
  gROOT->ProcessLine("neut_1p(\"../newg4/neut/neutvect_numu.root\",\"neut_1p_ext_numu.root\")");
  gROOT->ProcessLine("neut_1p(\"../newg4/neut/neutvect_numubar.root\",\"neut_1p_ext_numubar.root\")"); }
