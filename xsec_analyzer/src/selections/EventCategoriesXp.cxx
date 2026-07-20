// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

std::map< int, std::pair< std::string, int > > CC1muXp_MAP = {
  

  /*{ kUnknown, { "Unknown", kGray } },
  { kNuMuCC0p0pi_CCQE, { "CCmu0p0pi (CCQE)", kBlue - 2 } },
  { kNuMuCC0p0pi_CCMEC, { "CCmu0p0pi (CCMEC)", kBlue - 6 } },
  { kNuMuCC0p0pi_CCRES, { "CCmu0p0pi (CCRES)", kBlue - 9 } },
  { kNuMuCC0p0pi_Other, { "CCmu0p0pi (Other)", kBlue - 10 } },
  { kNuMuCC1p0pi_CCQE, { "CCmu1p0pi (CCQE)", kOrange + 4 } },
  { kNuMuCC1p0pi_CCMEC, { "CCmu1p0pi (CCMEC)", kOrange + 5  } },
  { kNuMuCC1p0pi_CCRES, { "CCmu1p0pi (CCRES)", kOrange + 6 } },
  { kNuMuCC1p0pi_Other, { "CCmu1p0pi (Other)", kOrange + 7 } },
  { kNuMuCC2p0pi_CCQE, { "CCmu2p0pi (CCQE)", kCyan - 3 } },
  { kNuMuCC2p0pi_CCMEC, { "CCmu2p0pi (CCMEC)", kCyan - 6 } },
  { kNuMuCC2p0pi_CCRES, { "CCmu2p0pi (CCRES)", kCyan - 4 } },
  { kNuMuCC2p0pi_Other, { "CCmu2p0pi (Other)", kCyan - 9 } },
  { kNuMuCCMp0pi_CCQE, { "CCmuMp0pi (CCQE)", kGreen } },
  { kNuMuCCMp0pi_CCMEC, { "CCmuMp0pi (CCMEC)", kGreen + 1 } },
  { kNuMuCCMp0pi_CCRES, { "CCmuMp0pi (CCRES)", kGreen + 2 } },
  { kNuMuCCMp0pi_Other, { "CCmuMp0pi (Other)", kGreen + 3 } },
  { kNuMuCCNpi, { "#nu_{#mu} CCN#pi", kAzure - 2 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }*/

{ kNumuCC_sig, {"#nu_{#mu} CC (Signal)", kMagenta-9 } },

{ kEXT, { "EXT", kGray + 2 } },

{ kOOFV, {"Out FV", kSpring -2 } },

// { kNCPi0, { "NC #pi^{0}", kCyan+2 } },

 //{ kNCOther, { "NC", kYellow-7 } },

// { kNueCCPi0, { "#nu_{e} CC #pi^{0}", kRed+1 } },

 //{ kNueCCOther, {"#nu_{e} CC Other", kOrange-3 } },

//{ kNumuCCPi0, { "#nu_{#mu} CC #pi^{0}", kBlue-2 } },

{ kUnknown, { "Unknown", kBlack } },

// { kOOFV, {"Out FV", kSpring -2 } },

// { kNumuCCOther, { "#nu_{#mu} CC Other", kBlue-7 } },

 { kNuECC, { "#nu_{e} CC ", kGreen+3 } },

 { kNC, { "NC", kCyan-4 } },

 { kOther, { "Other", kAzure+7 } }

/*  { kNueCCOther, {"#nu_{e} CC Other", kOrange-3 } },
  { kNumuCCPi0, { "#nu_{#mu} CC #pi^{0}", kBlue-2 } },
  { kNumuCCOther, { "#nu_{#mu} CC Other", kBlue-7 } },
  { kNueCCPi0, { "#nu_{e} CC #pi^{0}", kRed+1 } },
  { kNuECC, { "#nu_{e} CC ", kGreen+3 } },
  { kNCpi0, { "NC #pi^{0}", kCyan+2 } },

 { kOther, { "Other", kAzure+7 } },

  { kNC, { "NC", kCyan-4 } },
  { kNCOther, { "NC", kYellow-7 } },
  { kOOFV, {"Out FV", kSpring -2 } },
  { kEXT, { "EXT", kGray + 2 } },
  { kUnknown, { "Unknown", kBlack } }*/
};
