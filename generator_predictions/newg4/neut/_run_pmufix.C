{
  gSystem->AddIncludePath("-I/home/t2k/nowak/generators/neut/src/neutclass");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvtx.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutpart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsipart.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutfsivert.so");
  gSystem->Load("/home/t2k/nowak/generators/neut/src/neutclass/neutvect.so");
  gROOT->ProcessLine(".L neut_cc1pi_pmufix.C+");
  gROOT->ProcessLine("neut_cc1pi_pmufix(\"neutvect_numu.root\",\"../../pmufix/neut_numu.root\")");
  gROOT->ProcessLine("neut_cc1pi_pmufix(\"neutvect_numubar.root\",\"../../pmufix/neut_numubar.root\")");
}
