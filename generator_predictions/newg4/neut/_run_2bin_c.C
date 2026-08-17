{
  gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi_2bin.C+");
  gROOT->ProcessLine("neut_cc1pi_2bin(\"neutvect_numu.root\",\"../../ppi2bin/neut_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi_2bin(\"neutvect_numubar.root\",\"../../ppi2bin/neut_numubar.root\")");
}
