
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
#include "XSecAnalyzer/Selections/NewPIDSelections.hh"


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
  // study selections with the particle classifier as pion identification
  else if ( selection_name == "CC1mu1piXpNewPID" ) {
    sel = new CC1mu1piXpNewPID;
  }
  else if ( selection_name == "CC1mu2piNewPID" ) {
    sel = new CC1mu2piNewPID;
  }
  else if ( selection_name == "CC1mu3piNewPID" ) {
    sel = new CC1mu3piNewPID;
  }
  // released single-pion selection with the per-candidate PID diagnostic
  else if ( selection_name == "CC1mu1piXpPIDDiag" ) {
    sel = new CC1mu1piXpPIDDiag;
  }
  // released single-pion selections without the Bragg-pion cut
  else if ( selection_name == "CC1mu1piXpNoBragg" ) {
    sel = new CC1mu1piXpNoBragg;
  }
  else if ( selection_name == "CC1mu1pi1pNoBragg" ) {
    sel = new CC1mu1pi1pNoBragg;
  }
  // current two-pion selection with the event-level track dump (Phase 2 of the multi-pion plan)
  else if ( selection_name == "CC1mu2piEvt" ) {
    sel = new CC1mu2piEvt;
  }
  // threshold variants "<NewPID selection>_tNN": P(pi) > NN/100
  else if ( selection_name.rfind( "NewPID_t" ) != std::string::npos ) {
    const size_t pos = selection_name.rfind( "_t" );
    const std::string base = selection_name.substr( 0, pos );
    const double cut = std::stod( selection_name.substr( pos + 2 ) ) / 100.;
    if ( base == "CC1mu1piXpNewPID" ) sel = new CC1mu1piXpNewPID( cut, selection_name );
    else if ( base == "CC1mu2piNewPID" ) sel = new CC1mu2piNewPID( cut, selection_name );
    else if ( base == "CC1mu3piNewPID" ) sel = new CC1mu3piNewPID( cut, selection_name );
    else {
      std::cerr << "Selection name requested: " << selection_name << " has no NewPID base\n";
      throw;
    }
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
