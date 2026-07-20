void flux_script() {
std::string filename = "flux_newg4.root";

TFile* file = TFile::Open(filename.c_str(), "READ");

TH1D* histNumu = file->Get<TH1D>("numu_cv_fhc");
TH1D* histNumubar = file->Get<TH1D>("numubar_cv_fhc");
int E_min = 13; // bin for 60 MeV cut == 13

double totalPOT = 3.28e20; // Run 1 FHC NuMI
//double totalPOT = 1.27e20; // Run 2 FHC NuMI
//double totalPOT = 2.08e20; // Run 4 FHC NuMI
//double totalPOT = 2.23e20; // Run 4 FHC NuMI




float flux_numu = histNumu->Integral(E_min, histNumu->GetNbinsX()) / totalPOT;
float flux_numubar = histNumubar->Integral(E_min, histNumubar->GetNbinsX()) / totalPOT;
float flux_both = flux_numu + flux_numubar;

std::cout << "Neutrino energy cut: " << histNumu->GetBinLowEdge(E_min) << " GeV" << std::endl;
std::cout << "Integrated Flux:" << std::endl
<< " Muon neutrino = " << flux_numu << " nu / POT / cm2" << std::endl
<< " Muon antineutrino = " << flux_numubar << " nu / POT / cm2" << std::endl
<< " Both = " << flux_both << " nu / POT / cm2" << std::endl;
}
