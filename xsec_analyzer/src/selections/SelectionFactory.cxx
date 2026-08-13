
// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0pi.hh"
#include "XSecAnalyzer/Selections/CC1muNp0pi.hh"
#include "XSecAnalyzer/Selections/NuMICC1e.hh"
#include "XSecAnalyzer/Selections/DummySelection.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"
#include "XSecAnalyzer/Selections/CC1mu1pi1p.hh"
#include "XSecAnalyzer/Selections/CC1mu2pi.hh"
#include "XSecAnalyzer/Selections/CC1mu3pi.hh"


SelectionFactory::SelectionFactory() {
}

SelectionBase* SelectionFactory::CreateSelection(
  const std::string& selection_name )
{
  SelectionBase* sel;
  if ( selection_name == "CC1mu1p0pi" ) {
    sel = new CC1mu1p0pi;
  }
  else if ( selection_name == "CC1mu2p0pi" ) {
    sel = new CC1mu2p0pi;
  }
  else if ( selection_name == "CC1muNp0pi" ) {
    sel = new CC1muNp0pi;
  }
  else if ( selection_name == "NuMICC1e" ) {
    sel = new NuMICC1e;
  }
  else if ( selection_name == "Dummy" ) {
    sel = new DummySelection;
  }
else if ( selection_name == "CC1mu1piXp" ) {
    sel = new CC1mu1piXp;
  }
  else if ( selection_name == "CC1mu1pi1p" ) {
    sel = new CC1mu1pi1p;
  }
  else if ( selection_name == "CC1mu2pi" ) {
    sel = new CC1mu2pi;
  }
  else if ( selection_name == "CC1mu3pi" ) {
    sel = new CC1mu3pi;
  }

  else {
    std::cerr << "Selection name requested: " << selection_name
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  // Ensure that the owned map of category definitions is set up
  sel->define_category_map();

  return sel;
}
