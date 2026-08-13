// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"
#include "TStyle.h"
#include "XSecAnalyzer/Selections/CC1mu1piXp.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

#include <TMVA/Reader.h>
#include "TFile.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "TLine.h"
#include <TH1D.h>
// NOTE: RooUnfold includes removed — unfolding is handled by the dedicated
// Unfolder stage, not inside the selection.
#include "TCanvas.h"
#include "TLegend.h"

namespace {
  // z dead region excluded from the neutrino-vertex fiducial volume, matching
  // the custom selection (FV_new.h deadzmin/deadzmax). Applies to the neutrino
  // VERTEX FV only — not to track containment.
  constexpr double DEAD_Z_MIN = 675.1; // cm
  constexpr double DEAD_Z_MAX = 775.1; // cm

  // point_inside_FV plus the z dead-region veto, for neutrino-vertex checks.
  inline bool in_vertex_FV( const FiducialVolume& fv,
                            double x, double y, double z ) {
    if ( !point_inside_FV( fv, x, y, z ) ) return false;
    if ( z > DEAD_Z_MIN && z < DEAD_Z_MAX ) return false;
    return true;
  }
}

CC1mu1piXp::CC1mu1piXp() : SelectionBase( "CC1mu1piXp" ) {
  CalcType = kOpt1;
  fPP = new TF1( "fPP", "29.354172*x-14.674918", 0.3, 0.5 );

/*double bins[27] = {
        0.100,0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,
        0.400,0.425,0.450,0.475,0.500,0.550,0.600,0.650,0.700,0.750,
        0.800,0.850,0.900,0.950,1.000,1.100,1.200
    };*/

/*double bins[31] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.2, 2.3, 2.4,
    2.5, 2.6, 2.7, 2.8, 2.9,
    3.0
};*/
double bins[21] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9,
    2.0
};

double bins_mu_shifted[21] = {
    0.025, 0.125, 0.225, 0.325, 0.425,
    0.525, 0.625, 0.725, 0.825, 0.925,
    1.025, 1.125, 1.225, 1.325, 1.425,
    1.525, 1.625, 1.725, 1.825, 1.925,
    2.025
};

double bins_pi_shifted[21] = {
    0.075, 0.175, 0.275, 0.375, 0.475,
    0.575, 0.675, 0.775, 0.875, 0.975,
    1.075, 1.175, 1.275, 1.375, 1.475,
    1.575, 1.675, 1.775, 1.875, 1.975,
    2.075
};

    //h_selected   = new TH1D("h_selected",   "Selected",   30, bins);
    h_selected = new TH1D("h_selected","Selected Events; p_{#mu}^{reco} [GeV/c]; Events",20, bins);
    //h_background = new TH1D("h_background", "Background", 26, bins);
    h_background = new TH1D("h_background","Selected Background; p_{#mu}^{reco} [GeV/c]; Events",20, bins);

    h_selected->Sumw2();
    h_background->Sumw2();

    //h_eff_num = new TH1D("h_eff_num", "Efficiency Numerator", 26, bins);
    h_eff_num = new TH1D("h_eff_num","Efficiency Numerator; p_{#mu}^{true} [GeV/c]; Events",20, bins);
    h_eff_den = new TH1D("h_eff_den","Efficiency Denominator; p_{#mu}^{true} [GeV/c]; Events",20, bins);
    //h_eff_den = new TH1D("h_eff_den", "Efficiency Denominator", 26, bins);

    h_eff_num->Sumw2();
    h_eff_den->Sumw2();

    h_all_signal   = new TH1D("h_all_signal",   "All True Signal; p_{#mu}^{true} [GeV/c]; Events", 20, bins);
    h_all_selected = new TH1D("h_all_selected", "All Selected Events; p_{#mu}^{reco} [GeV/c]; Events", 20, bins);
    h_sel_signal   = new TH1D("h_sel_signal",   "Selected Signal Events; p_{#mu}^{reco} [GeV/c]; Events", 20, bins);
    h_sel_bkg      = new TH1D("h_sel_bkg",      "Selected Background Events; p_{#mu}^{reco} [GeV/c]; Events", 20, bins);
    h_purity       = new TH1D("h_purity",       "Purity; p_{#mu}^{reco} [GeV/c]; Purity", 20, bins);
    ///h_response_full = new TH2D("h_response_full","Response Matrix; p_{#mu}^{true} [GeV/c]; p_{#mu}^{reco} [GeV/c]",20, bins,20, bins);
 h_response_full = new TH2D(
    "h_response_full",
    "Response Matrix; p_{#mu}^{reco} [GeV/c]; p_{#mu}^{true} [GeV/c]",
    20, bins,
    20, bins
);
    h_all_signal->Sumw2();
    h_all_selected->Sumw2();
    h_sel_signal->Sumw2();
    h_sel_bkg->Sumw2();
    h_purity->Sumw2();

    // =====================================
// Muon efficiency histograms
// =====================================

h_mu_eff_num = new TH1D(
    "h_mu_eff_num",
    "Muon Efficiency Numerator; p_{#mu}^{true} [GeV/c]; Events",
    20, bins_mu_shifted
);

h_mu_eff_den = new TH1D(
    "h_mu_eff_den",
    "Muon Efficiency Denominator; p_{#mu}^{true} [GeV/c]; Events",
    20, bins_mu_shifted
);

h_mu_eff_num->Sumw2();
h_mu_eff_den->Sumw2();


// =====================================
// Pion efficiency histograms
// =====================================

h_pi_eff_num = new TH1D(
    "h_pi_eff_num",
    "Pion Efficiency Numerator; p_{#pi}^{true} [GeV/c]; Events",
    20, bins_pi_shifted
);

h_pi_eff_den = new TH1D(
    "h_pi_eff_den",
    "Pion Efficiency Denominator; p_{#pi}^{true} [GeV/c]; Events",
    20, bins_pi_shifted
);

h_pi_eff_num->Sumw2();
h_pi_eff_den->Sumw2();

h_mu_eff = new TH1D(
    "h_mu_eff",
    "Muon Efficiency; p_{#mu}^{true} [GeV/c]; Efficiency",
    20, bins_mu_shifted
);

h_mu_eff->Sumw2();

h_pi_eff = new TH1D(
    "h_pi_eff",
    "Pion Efficiency; p_{#pi}^{true} [GeV/c]; Efficiency",
    20, bins_pi_shifted
);

h_pi_eff->Sumw2();

h_mu_cut0 = new TH1D(
"h_mu_cut0",
"No Cuts;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut0->Sumw2();

h_mu_cut1_vertex =
new TH1D(
"h_mu_cut1_vertex",
"Vertex Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut1_vertex->Sumw2();

h_mu_cut2_topology =
new TH1D(
"h_mu_cut2_topology",
"Topology Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut2_topology->Sumw2();

h_mu_cut3_tracklike =
new TH1D(
"h_mu_cut3_tracklike",
"Tracklike Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut3_tracklike->Sumw2();

h_mu_cut4_pioncontained =
new TH1D(
"h_mu_cut4_pioncontained",
"Pion Contained Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut4_pioncontained->Sumw2();

h_mu_cut5_muongap =
new TH1D(
"h_mu_cut5_muongap",
"Muongap Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut5_muongap->Sumw2();

h_mu_cut6_piongap =
new TH1D(
"h_mu_cut6_piongap",
"Piongap Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut6_piongap->Sumw2();

h_mu_cut7_shower =
new TH1D(
"h_mu_cut7_shower",
"Shower Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut7_shower->Sumw2();

h_mu_cut8_opening =
new TH1D(
"h_mu_cut8_opening",
"Openingangle Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut8_opening->Sumw2();

h_mu_cut9_final =
new TH1D(
"h_mu_cut9_final",
"Final Cut;p_{#mu}^{true} [GeV/c];Events",
20,bins_mu_shifted);

h_mu_cut9_final->Sumw2();

h_pi_cut0 = new TH1D(
"h_pi_cut0",
"No Cuts;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut0->Sumw2();

h_pi_cut1_vertex =
new TH1D(
"h_pi_cut1_vertex",
"Vertex Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut1_vertex->Sumw2();

h_pi_cut2_topology =
new TH1D(
"h_pi_cut2_topology",
"Topology Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut2_topology->Sumw2();

h_pi_cut3_tracklike =
new TH1D(
"h_pi_cut3_tracklike",
"Tracklike Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut3_tracklike->Sumw2();

h_pi_cut4_pioncontained =
new TH1D(
"h_pi_cut4_pioncontained",
"Pion Contained Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut4_pioncontained->Sumw2();

h_pi_cut5_muongap =
new TH1D(
"h_pi_cut5_muongap",
"Muongap Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut5_muongap->Sumw2();

h_pi_cut6_piongap =
new TH1D(
"h_pi_cut6_piongap",
"Piongap Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut6_piongap->Sumw2();

h_pi_cut7_shower =
new TH1D(
"h_pi_cut7_shower",
"Shower Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut7_shower->Sumw2();

h_pi_cut8_opening =
new TH1D(
"h_pi_cut8_opening",
"Openingangle Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut8_opening->Sumw2();

h_pi_cut9_final =
new TH1D(
"h_pi_cut9_final",
"Final Cut;p_{#pi}^{true} [GeV/c];Events",
20,bins_pi_shifted);

h_pi_cut9_final->Sumw2();

// Per-cut RECO cut-flow histograms: 5 observables x 10 stages, all events.
{
  const char* cf_obs[5]   = {"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* cf_xttl[5]  = {"p_{#mu}^{reco} [GeV/c]","p_{#pi}^{reco} [GeV/c]",
                             "cos#theta_{#mu}^{reco}","cos#theta_{#pi}^{reco}",
                             "#theta_{#mu#pi}^{reco} [rad]"};
  const char* cf_stage[10]= {"cut0_none","cut1_vertex","cut2_topology","cut3_tracklike",
                             "cut4_pioncontained","cut5_muongap","cut6_piongap",
                             "cut7_shower","cut8_opening","cut9_final"};
  for (int o = 0; o < 5; ++o) {
    for (int c = 0; c < 10; ++c) {
      TString hn  = Form("h_cf_%s_%s", cf_obs[o], cf_stage[c]);
      TString ttl = Form("%s;%s;Events", cf_stage[c], cf_xttl[o]);
      if      (o == 0) h_cf[o][c] = new TH1D(hn, ttl, 20, bins_mu_shifted);
      else if (o == 1) h_cf[o][c] = new TH1D(hn, ttl, 20, bins_pi_shifted);
      else if (o == 2 || o == 3) h_cf[o][c] = new TH1D(hn, ttl, 20, -1.0, 1.0);
      else             h_cf[o][c] = new TH1D(hn, ttl, 20, 0.0, 3.15);
      h_cf[o][c]->Sumw2();
    }
  }
}

// Cut-flow yield counters (10 stages: 0=none,1=vertex,2=topology,3=tracklike,
// 4=pioncontained,5=muongap,6=piongap,7=shower,8=opening,9=final)
h_cutflow_tot = new TH1D("h_cutflow_tot", "cut-flow total;cut stage;events", 10, 0, 10);
h_cutflow_sig = new TH1D("h_cutflow_sig", "cut-flow signal;cut stage;events", 10, 0, 10);
h_cutflow_tot->Sumw2();
h_cutflow_sig->Sumw2();

// Selection-diagnostic histograms (signal [0] / background [1])
{
  const char* sb[2] = {"sig","bkg"};
  for (int k = 0; k < 2; ++k) {
    h_nm1_topo[k]  = new TH1D(Form("h_nm1_topo_%s",sb[k]),  ";topological score;events", 50, 0, 1);
    h_nm1_oa[k]    = new TH1D(Form("h_nm1_oa_%s",sb[k]),    ";#theta_{#mu#pi} [rad];events", 40, 0, 3.15);
    h_fin_mupid[k] = new TH1D(Form("h_fin_mupid_%s",sb[k]), ";muon LLR PID;events", 40, -1, 1);
    h_fin_pipid[k] = new TH1D(Form("h_fin_pipid_%s",sb[k]), ";pion LLR PID;events", 40, -1, 1);
    h_fin_mulen[k] = new TH1D(Form("h_fin_mulen_%s",sb[k]), ";muon length [cm];events", 50, 0, 300);
    h_fin_pilen[k] = new TH1D(Form("h_fin_pilen_%s",sb[k]), ";pion length [cm];events", 50, 0, 200);
    h_nm1_topo[k]->Sumw2(); h_nm1_oa[k]->Sumw2();
    h_fin_mupid[k]->Sumw2(); h_fin_pipid[k]->Sumw2();
    h_fin_mulen[k]->Sumw2(); h_fin_pilen[k]->Sumw2();
  }
  h_bkgcat = new TH1D("h_bkgcat", ";event category;events", 30, 0, 30);
  h_bkgcat->Sumw2();
}


}

void CC1mu1piXp::define_constants() {
  this->define_true_FV( 10., 246., -101., 101., 10., 986. );
  this->define_reco_FV( 10., 246., -101., 101., 10., 986. );
  // Containment volume matches the custom selection's isContained() box
  // (FV_new.h: FVx/y/z_1 boundaries), which is larger than the fiducial box
  // and is the volume used to decide pion containment and muon range-vs-MCS.
  this->define_containment_FV( 2.0, 254.35, -113.53, 107.47, 2.1, 1034.9 );
       
	tmvaReader = new TMVA::Reader();
        tmvaReader->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva);
        tmvaReader->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva);
        tmvaReader->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva);
        tmvaReader->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva);
        tmvaReader->AddVariable("trk_score_v", &trk_score_v_tmva);
        //      tmvaReader->AddVariable("trk_len_v", &trk_len_v_tmva);
        tmvaReader->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva);
        tmvaReader->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva);
        tmvaReader->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva);
        //tmvaReader->AddVariable("isContained_out", &isContained_tmva);
        //
        tmvaReader_mu = new TMVA::Reader();
        tmvaReader_mu->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_score_v", &trk_score_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_len_v", &trk_len_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_mu);
        tmvaReader_mu->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_mu);

        tmvaReader_pi = new TMVA::Reader();
        tmvaReader_pi->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_score_v", &trk_score_v_tmva_pi);
        //      tmvaReader_pi->AddVariable("trk_len_v", &trk_len_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_pi);
        tmvaReader_pi->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_pi);

        // BDT weight files. These were hardcoded under /home/lar/ipophale/,
        // i.e. another user's home directory, so the selection only ran for
        // whoever could read that path. Resolve them instead from, in order:
        //
        //   1. $CC1MU1PIXP_BDT_WEIGHTS_DIR, if set
        //   2. <XSEC_ANALYZER_DIR>/../booster_decision_tree  (the in-tree copy;
        //      verified byte-identical to the original files)
        //   3. the original absolute path, so existing setups keep working
        //
        // and fail with a clear message rather than letting TMVA book an empty
        // reader if none of them resolve.
        auto resolve_bdt_dir = []() -> std::string {
          if ( const char* env_dir = std::getenv("CC1MU1PIXP_BDT_WEIGHTS_DIR") ) {
            return std::string( env_dir );
          }
          if ( const char* xa_dir = std::getenv("XSEC_ANALYZER_DIR") ) {
            std::string candidate = std::string( xa_dir )
              + "/../booster_decision_tree";
            struct stat sb;
            if ( stat( candidate.c_str(), &sb ) == 0 && S_ISDIR( sb.st_mode ) ) {
              return candidate;
            }
          }
          return "/home/lar/ipophale/booster_decision_tree";
        };

        const std::string bdt_dir = resolve_bdt_dir();

        auto book_bdt = [ &bdt_dir ]( TMVA::Reader* reader,
          const std::string& dataset )
        {
          std::string path = bdt_dir + "/" + dataset
            + "/weights/TMVAClassification_BDT.weights.xml";
          std::ifstream test( path );
          if ( !test.good() ) {
            throw std::runtime_error( "CC1mu1piXp: could not read BDT weight"
              " file \"" + path + "\". Set CC1MU1PIXP_BDT_WEIGHTS_DIR to the"
              " directory containing the dataset_* folders." );
          }
          reader->BookMVA( "BDT", path.c_str() );
        };

        book_bdt( tmvaReader,    "dataset_MIP_BDT_no_len" );
        book_bdt( tmvaReader_mu, "dataset_muon_BDT" );
        book_bdt( tmvaReader_pi, "dataset_pion_BDT_no_len" );

}

int CC1mu1piXp::categorize_event( AnalysisEvent* Event ) {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  //



 // is this EXT only? can change name
   /* int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
    Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO ||
      abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
    if ( !Event->is_mc_ ) {
      return kEXT;
    }*/

  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

 // We shouldn't ever get here, but return "unknown" just in case
   // std::cout << "Warning: Unknown event! Check the categorization logic." << std::endl;
   // return kUnknown;

//std::cout<<"This is fine 1" <<std::endl;

 bool MCVertexInFV = in_vertex_FV( this->true_FV(),
    Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
  if ( !MCVertexInFV ) {
    return kOOFV;
  }
   // return kNuMuCCOther;
   //
  if ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
       // if(Event->mc_npi0_ > 0) return kNCPi0;
        return kNC;
    }

if ( std::abs(Event->mc_nu_pdg_) == ELECTRON_NEUTRINO ) {

	return kNuECC;

//        if(Event->mc_npi0_ > 0) return kNueCCPi0;
        //else return kNueCCOther;
    }

// Must match the isSignal boolean in define_signal(): include the
// true-vertex-in-FV and heavy-meson veto terms so this kNumuCC_sig category
// agrees with is_event_mc_signal().
bool Is_CC1mu1piXp_event = sig_truevertex_in_fv_ && sig_ccnc_ && sig_is_numu_ && sig_one_muon_above_thresh_ && sig_one_charged_pion_
&& sig_no_pions_  && sig_no_kaons_ && sig_no_heavy_mesons_ && TrueCandidateMuonP.Mag() > 0.15 && TrueCandidatePionP.Mag() > 0.175
&& TrueCandidateMuonP.Angle( TrueCandidatePionP ) < 2.6;  // measured phase space: pmu>0.15, ppi>0.175 GeV/c, theta_mupi<2.6 rad

if (std::abs(Event->mc_nu_pdg_) == MUON_NEUTRINO) {
        // signal events
        if ( Is_CC1mu1piXp_event) return kNumuCC_sig;
        // non-signal nues
        //else return kNumuCCOther;
   
}
//bool isNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );
  // DB Currently only one NC category is supported so test first.
  // Will likely want to change this in the future
 // if (isNC) return kNC;

  //if ( Event->mc_nu_pdg_ == ELECTRON_NEUTRINO ) {
   // return kNuECC;
  //}
 /* if ( !(Event->mc_nu_pdg_ == MUON_NEUTRINO) ) {
    return kOther;
  }*/
 //std::cout << "Warning: Unknown event! Check the categorization logic." << std::endl;
   return kUnknown;


// Boolean which basically MC Signal selection without requesting a
  // particular number of protons (N >= 1)
  /*bool Is_CC1muNp0pi_Event = (sig_mc_n_threshold_proton >= 1)
    && sig_no_pions_ && sig_one_muon_above_thresh_ && sig_no_heavy_mesons_;

  if ( Is_CC1muNp0pi_Event ) {
    if (sig_mc_n_threshold_proton == 1) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC1p0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCC1p0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCC1p0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
} else if (sig_mc_n_threshold_proton == 2) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC2p0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCC2p0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCC2p0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
 } else { // i.e. >=3
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCCMp0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCCMp0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCCMp0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
    }
  }
else if (!sig_no_pions_) {
    return kNuMuCCNpi;
  } else if (sig_mc_n_threshold_proton == 0) {
    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
    else if ( Event->mc_nu_interaction_type_ == 10 ) {
      return kNuMuCC0p0pi_CCMEC; // MEC
    }
    else if ( Event->mc_nu_interaction_type_ == 1 ) {
      return kNuMuCC0p0pi_CCRES; // RES
    }
    else return kNuMuCC0p0pi_Other;
  }
  return kNuMuCCOther;
*/

// Boolean which basically MC Signal selection without requesting a
  // particular number of protons (N >= 1)
  //bool Is_CC1mu1piXp_Event = (sig_mc_n_threshold_proton >= 1)
    //&& sig_has_a_pion_ && sig_one_muon_above_thresh_ && sig_no_heavy_mesons_;
//what interaction type to use here? There is nO res with 1pi
//std::cout<<"This is fine 2"<< std::endl;
}//categorise event loop

void CC1mu1piXp::compute_reco_observables( AnalysisEvent* Event ) {

if(CandidateMuonIndex != -1){


candidate_muon_mom_mcs = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);

// Reco muon cos(theta) w.r.t. the detector z-axis from the track direction.
TVector3 mu_dir( Event->track_dirx_->at(CandidateMuonIndex),
                 Event->track_diry_->at(CandidateMuonIndex),
                 Event->track_dirz_->at(CandidateMuonIndex) );
candidate_muon_costh_reco = mu_dir.CosTheta();

}

if(CandidatePionIndex != -1){

// Reco pion momentum: muon-hypothesis range momentum (trk_range_muon_mom_v).
// The muon mass (105.7 MeV) is close to the pion mass (139.6 MeV), so the
// range->momentum conversion is a good proxy, and -- crucially -- it is a
// *momentum*, directly comparable to candidate_pion_mom_true below.
//
// This previously used trk_energy_proton_v, a proton-hypothesis range kinetic
// energy. That is a different physical quantity from the true momentum it is
// binned against in ccpi_ppi_bin_config.txt (identical bin edges), so the
// migration matrix was systematically shifted rather than near-diagonal: a
// true 0.2 GeV/c pion gives KE ~0.104 GeV, two bins low, and a true
// 0.175 GeV/c pion falls below the first reco bin entirely.
//
// NOTE: range momentum is only meaningful for contained tracks. Uncontained
// pions are handled by the sel_pion_contained requirement in the selection.
candidate_pion_mom_reco = Event->track_range_mom_mu_->at(CandidatePionIndex);

// Reco pion cos(theta) w.r.t. the detector z-axis from the track direction.
TVector3 pi_dir( Event->track_dirx_->at(CandidatePionIndex),
                 Event->track_diry_->at(CandidatePionIndex),
                 Event->track_dirz_->at(CandidatePionIndex) );
candidate_pion_costh_reco = pi_dir.CosTheta();

}

}

void CC1mu1piXp::compute_true_observables( AnalysisEvent* Event ) {

if(truemuonindex != -1){

double CandidateMuonPx = Event->mc_nu_daughter_px_->at(truemuonindex);
double CandidateMuonPy = Event->mc_nu_daughter_py_->at(truemuonindex);
double CandidateMuonPz = Event->mc_nu_daughter_pz_->at(truemuonindex);


//double CandidateMuonPx = Event->pfp_true_px_->at(truemuonindex);

//double CandidateMuonPy = Event->pfp_true_py_->at(truemuonindex);
//double CandidateMuonPz = Event->pfp_true_pz_->at(truemuonindex);
TVector3 BackTrackCandidateMuonP(CandidateMuonPx, CandidateMuonPy, CandidateMuonPz );
candidate_muon_mom_true = BackTrackCandidateMuonP.Mag();
candidate_muon_costh_true = BackTrackCandidateMuonP.CosTheta();

}

// True pion momentum / cos(theta) and the true mu-pi opening angle, computed
// from the truth candidate vectors filled in define_signal() (which runs before
// this method).  Guarded so they are only filled when a true muon and a true
// charged pion exist (i.e. signal-like events); non-signal true bins do not use
// these observable values.  NOTE: true_mu_pi_opening_angle was previously left
// at 0 here, which made every signal event fall in the first theta_mu-pi true
// bin — fixed now.
if ( truemuonindex != -1 && sig_mc_n_threshold_pionpm >= 1 ) {
  candidate_pion_mom_true   = TrueCandidatePionP.Mag();
  candidate_pion_costh_true = TrueCandidatePionP.CosTheta();
  true_mu_pi_opening_angle  = TrueCandidateMuonP.Angle( TrueCandidatePionP );
}

}

bool CC1mu1piXp::define_signal( AnalysisEvent* Event ) {
  // Event weight for the diagnostic histograms below. These were all filled
  // with raw counts, so the efficiency and purity printed by finalize(), and
  // the response matrix, were untuned-CV numbers not comparable to anything
  // the Unfolder produces (which applies this same product via UniverseMaker).
  // Both members default to DEFAULT_WEIGHT = 1, so data and EXT are unaffected.
  const double evt_w = Event->spline_weight_ * Event->tuned_cv_weight_;


//TMVA stuff

//std::cout << "Signal called\n";



//std::cout<<"Signal is being called" <<std::endl;

  // =====================================
  // DB Calculate the values which we need
  sig_mc_n_threshold_muon = 0;
  sig_mc_n_threshold_proton = 0;
  sig_mc_n_threshold_pion0 = 0;
  sig_mc_n_threshold_pionpm = 0;
  sig_mc_n_heaviermeson = 0;
  sig_mc_n_kaons =0;
  //TVector3 TrueCandidateMuonP;
  //TVector3 TrueCandidatePionP;

//std::cout<<"This is fine 3"<<std::endl;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    TVector3 MCParticle( Event->mc_nu_daughter_px_->at(p),
      Event->mc_nu_daughter_py_->at(p), Event->mc_nu_daughter_pz_->at(p) );
    double CandidateTrueMuonPx = Event->mc_nu_daughter_px_->at(p);
    double CandidateTrueMuonPy = Event->mc_nu_daughter_py_->at(p);
    double CandidateTrueMuonPz = Event->mc_nu_daughter_pz_->at(p);
   // TVector3 TrueCandidateMuonP;//(CandidateTrueMuonPx, CandidateTrueMuonPy, CandidateTrueMuonPz );

//   std::cout<<"	Signal A "<< std::endl;

    double CandidateTruePionPx = Event->mc_nu_daughter_px_->at(p);
    double CandidateTruePionPy = Event->mc_nu_daughter_py_->at(p);
    double CandidateTruePionPz = Event->mc_nu_daughter_pz_->at(p);
   // TVector3 TrueCandidatePionP;//(CandidateTruePionPx, CandidateTruePionPy, CandidateTruePionPz );

     double mc_muon_momentum = -999.9;
     double mc_pion_momentum = -999.9;
     //TVector3 mc_muon_mom;
     //TVector3 mc_pion_mom;
// std::cout<<" Signal B "<< std::endl;


    double ParticleMomentum = MCParticle.Mag();
//std::cout<<"This is fine 4" <<std::endl;

    if ( std::abs(pdg) == MUON && ParticleMomentum > 0.0 ) {

// std::cout<<" Signal C "<< std::endl;

      sig_mc_n_threshold_muon++;
      truemuonindex = p;
      //TrueCandidateMuonP.SetXYZ(CandidateTrueMuonPx, CandidateTrueMuonPy, CandidateTrueMuonPz );
      TrueCandidateMuonP.SetXYZ(Event->mc_nu_daughter_px_->at(p),Event->mc_nu_daughter_py_->at(p),Event->mc_nu_daughter_pz_->at(p));
//std::cout<<"True Candidate Muon " << TrueCandidateMuonP.Mag() << std::endl;

    }
    else if ( pdg == PROTON && ParticleMomentum >= 0.0 ) {
      sig_mc_n_threshold_proton++;
      trueprotonindex = p;
    }
    else if ( pdg == PI_ZERO ) {
      sig_mc_n_threshold_pion0++;
    }
    else if ( std::abs(pdg) == PI_PLUS && ParticleMomentum > 0.0 ) { //changed from 0.07 to 0.1
      sig_mc_n_threshold_pionpm++;
      TrueCandidatePionP.SetXYZ(Event->mc_nu_daughter_px_->at(p),Event->mc_nu_daughter_py_->at(p),Event->mc_nu_daughter_pz_->at(p));


    }

    else if ( (std::abs(pdg) == 321) || (pdg == 310) || (pdg == 130) || (pdg == 311)){
	sig_mc_n_kaons++; //look at PDG code list
    }
	
    else if ( pdg != PI_ZERO && std::abs(pdg) != PI_PLUS
      && is_meson_or_antimeson(pdg) )
    {
      sig_mc_n_heaviermeson++;
    }
//    if((sig_mc_n_threshold_muon ==1) && ( sig_mc_n_threshold_pionpm == 1)){
//	true_mu_pi_opening_angle = TrueCandidateMuonP.Angle(TrueCandidatePionP);
  //  }
// std::cout<<" Signal D "<< std::endl;


  }

  //  if((sig_mc_n_threshold_muon ==1) && ( sig_mc_n_threshold_pionpm == 1)){

//        true_mu_pi_opening_angle = TrueCandidateMuonP.Angle(TrueCandidatePionP);
 //   }
	
	 
// ===========================================================
  // Calculate the booleans related to the different signal cuts

  sig_truevertex_in_fv_ = in_vertex_FV( this->true_FV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );

  sig_recovertex_in_fv_ = in_vertex_FV( this->reco_FV(), Event->nu_vx_,
    Event->nu_vy_, Event->nu_vz_ );
//std::cout<<"This is fine 5" << std::endl;
//
// std::cout<<" Signal E "<< std::endl;


  sig_ccnc_= ( Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  sig_is_numu_ = (std::abs( Event->mc_nu_pdg_) == MUON_NEUTRINO );
  //std::cout<<"This is the mc_nu_pdg "<< Event->mc_nu_pdg_ <<std::endl;
  sig_one_muon_above_thresh_ = ( sig_mc_n_threshold_muon == 1 );
  sig_one_proton_above_thresh_ = ( sig_mc_n_threshold_proton > 0 ); //the name is misleading here we are saying any number of protons are allowed
  sig_no_pions_ =  (sig_mc_n_threshold_pion0 == 0);   //Misleading, has 1 pion
   // && (sig_mc_n_threshold_pionpm == 1) );
  sig_one_charged_pion_ =(sig_mc_n_threshold_pionpm == required_charged_pions());
  sig_no_heavy_mesons_ = ( sig_mc_n_heaviermeson == 0 );
  sig_no_kaons_ = (sig_mc_n_kaons ==0);

// ====================
  // Is the event signal?
  //
//  std::cout<<" Signal F "<< std::endl;



  bool isSignal =  sig_truevertex_in_fv_  && sig_ccnc_ && sig_is_numu_ && sig_one_muon_above_thresh_ && sig_one_charged_pion_
     && sig_no_pions_  && sig_no_kaons_ && TrueCandidateMuonP.Mag() > 0.15 && TrueCandidatePionP.Mag() > 0.175
     && TrueCandidateMuonP.Angle( TrueCandidatePionP ) < 2.6 && sig_no_heavy_mesons_; // measured phase space (pmu>0.15, ppi>0.175, theta_mupi<2.6)
     // && sig_mc_n_threshold_pionpm

     //
    // std::cout<<"Exiting signal here" <<std::endl;
    //
//   std::cout << "Signal breakdown: "
//          << sig_truevertex_in_fv_ << " "
//          << sig_ccnc_ << " "
//          << sig_is_numu_ << " "
//          << sig_one_muon_above_thresh_
//          << std::endl;

//   std::cout << "isSignal = " << isSignal << std::endl;
if(isSignal){

    if(sig_truevertex_in_fv_) sig_truevertex_fv++;
    if(sig_ccnc_)             sig_ccnc++;
    if(sig_is_numu_)          sig_numu++;
    if(sig_one_muon_above_thresh_) sig_one_muon++;

}

  //return isSignal; 
   //return true;
// std::cout<<" Signal G "<< std::endl;

//std::cout<<"!!!!!!!!!!!"<<std::endl;

if (isSignal) {

//std::cout<<"bfeore muon index "<<std::endl;

    if (truemuonindex != -1) {

        double muon_true = TrueCandidateMuonP.Mag();
       // double muon_true = TrueCandidateMuonP.Mag();
        double pion_true = TrueCandidatePionP.Mag();

	h_all_signal->Fill(muon_true, evt_w);
        h_eff_den->Fill(muon_true, evt_w);
         // denominator = all true signal muons
        h_mu_eff_den->Fill(muon_true, evt_w);
        h_pi_eff_den->Fill(pion_true, evt_w);
	h_mu_cut0->Fill(muon_true, evt_w);
        h_pi_cut0->Fill(pion_true, evt_w);
    }
//std::cout<<"after muon index" <<std::endl;
}

 return isSignal;

}


//std::cout <<"This is fine 6"<< std::endl;

bool CC1mu1piXp::selection( AnalysisEvent* Event) {
  // Event weight for the diagnostic histograms below. These were all filled
  // with raw counts, so the efficiency and purity printed by finalize(), and
  // the response matrix, were untuned-CV numbers not comparable to anything
  // the Unfolder produces (which applies this same product via UniverseMaker).
  // Both members default to DEFAULT_WEIGHT = 1, so data and EXT are unaffected.
  const double evt_w = Event->spline_weight_ * Event->tuned_cv_weight_;


//std::cout << "Selection called\n";

// std::cout<<" Selection A "<< std::endl;


//std::cout<< "This is fine 6" << std::endl;

  // Requirement for exactly one neutrino slice
  if (Event->nslice_ == 1) sel_nslice_eq_1_ = true;
  sel_topo_cut_passed_ = Event->topological_score_ > 0.67; //changed here

 // ======================================================================
  // Requirement for exactly 2 tracks and 0 showers for a CC1p0pi selection

  int reco_shower_count = 0;
  int reco_track_count = 0;
  std::vector<int> CandidateIndex;

// std::cout<<" Selection B "<< std::endl;


  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) {
      ++reco_shower_count;
    } else {
      ++reco_track_count;
      CandidateIndex.push_back(p);
    }

  } 

// std::cout<<" Selection C "<< std::endl;


  if (reco_shower_count == 0) sel_nshower_eq_0_ = true;
  // Require at least the muon+pion (2 tracks). Was ">2" (>=3 tracks), which
  // dropped the 0-proton 1mu+1pi topology; the custom selection
  // ccpi_selection.C has no minimum-track gate (only nPrimaryTracks<5 upper).
  if (reco_track_count >= 2) sel_ntrack_gt_2_ = true;

 // ==============================================
  // Identify candidate muon & proton
  // Muon = the one with the highest LLR PID Score
  // Proton = the one with the lowest LLR PID Score

  //int CandidateMuonIndex = -1;
  //int CandidateProtonIndex = -1;
  //int CandidatePionIndex = -1;

// std::cout<<"CandidateMuonIndex1 is  "<< CandidateMuonIndex << std::endl;

 
 //std::cout<<"This is fine 7" << std::endl;

  TVector3 reco_primary_vtx(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

  if (sel_ntrack_gt_2_) {

// std::cout<<" Selection D "<< std::endl;

    // Highest muon-BDT score seen so far, for the BNB-style muon-candidate choice.
    double best_muon_bdt_score = -1e30;

    for (size_t i_pfp = 0; i_pfp < Event->track_length_->size(); i_pfp++) {

	TVector3 track_start_ia (Event->track_startx_->at(i_pfp),Event->track_starty_->at(i_pfp),Event->track_startz_->at(i_pfp));
	      TVector3 nu_to_track_dist_ia (track_start_ia - reco_primary_vtx);

	if((Event->pfp_generation_->at(i_pfp) ==2)&&(Event->pfp_track_score_->at(i_pfp) >= 0.5)
	      && (Event->track_llr_pid_score_->at(i_pfp) > -1) && (Event->track_llr_pid_score_->at(i_pfp) < 2)
	      && (Event->trk_bragg_p_v->at(i_pfp)  > 0) && (Event->trk_bragg_p_v->at(i_pfp)  < 500) 
	       && (Event->trk_bragg_mu_v->at(i_pfp) > 0) && (Event->trk_bragg_mu_v->at(i_pfp) < 500)
	       && (Event->trk_bragg_mip_v->at(i_pfp)> 0) && (Event->trk_bragg_mip_v->at(i_pfp)< 500) 
	       && (Event->track_length_->at(i_pfp) < 1e6) && (Event->track_length_->at(i_pfp) > 0)
	       ){

	      trk_bragg_p_v_tmva = Event->trk_bragg_p_v->at(i_pfp);
	      trk_bragg_mu_v_tmva = Event->trk_bragg_mu_v->at(i_pfp);
	      trk_bragg_mip_v_tmva = Event->trk_bragg_mip_v->at(i_pfp);
	      trk_llr_pid_score_v_tmva = Event->track_llr_pid_score_->at(i_pfp);
	      trk_score_v_tmva = Event->pfp_track_score_->at(i_pfp);
	      trk_len_v_tmva = Event->track_length_->at(i_pfp);
	      trk_sce_end_x_v_tmva = Event->track_endx_->at(i_pfp);
	      trk_sce_end_y_v_tmva = Event->track_endy_->at(i_pfp);
	      trk_sce_end_z_v_tmva = Event->track_endz_->at(i_pfp);
	      tmvaOutput_mip = tmvaReader->EvaluateMVA("BDT");

	      // Muon BDT: previously booked (dataset_muon_BDT) but never evaluated.
	      // Activated here (BNB-style) to choose the muon candidate by the
	      // highest muon-BDT score instead of by longest track length. The muon
	      // BDT already includes track length among its input features.
	      trk_bragg_p_v_tmva_mu = Event->trk_bragg_p_v->at(i_pfp);
	      trk_bragg_mu_v_tmva_mu = Event->trk_bragg_mu_v->at(i_pfp);
	      trk_bragg_mip_v_tmva_mu = Event->trk_bragg_mip_v->at(i_pfp);
	      trk_llr_pid_score_v_tmva_mu = Event->track_llr_pid_score_->at(i_pfp);
	      trk_score_v_tmva_mu = Event->pfp_track_score_->at(i_pfp);
	      trk_len_v_tmva_mu = Event->track_length_->at(i_pfp);
	      trk_sce_end_x_v_tmva_mu = Event->track_endx_->at(i_pfp);
	      trk_sce_end_y_v_tmva_mu = Event->track_endy_->at(i_pfp);
	      trk_sce_end_z_v_tmva_mu = Event->track_endz_->at(i_pfp);
	      double muon_bdt_score = tmvaReader_mu->EvaluateMVA("BDT");

	// std::cout<<"CandidateMuonIndex2 is "<< CandidateMuonIndex << std::endl;


	if((Event->track_llr_pid_score_->at(i_pfp) > 0.2) && ( Event->track_length_->at(i_pfp) > 10)
	    && ( Event->pfp_track_score_->at(i_pfp) > 0.8) && (nu_to_track_dist_ia.Mag() <= 4)&&  (tmvaOutput_mip >= -0.1) ){
	//	if ( Event->pfp_reco_pdg_->at(i_pfp) == MUON ) {
      			sel_muoncandidate_tracklike_ = true;
    			//}
		// std::cout<<"CandidateMuonIndex3 is "<< CandidateMuonIndex << std::endl;

		// Pick the muon candidate with the HIGHEST muon-BDT score (BNB-style,
		// dataset_muon_BDT). Was previously the longest qualifying track.
		if(CandidateMuonIndex == -1
		   || (muon_bdt_score > best_muon_bdt_score)){
			CandidateMuonIndex = i_pfp;
			best_muon_bdt_score = muon_bdt_score;
		}
	}
	//std::cout<<"CandidateMuonIndex5 is "<< CandidateMuonIndex <<std::endl;
 }
 }

if(sel_muoncandidate_tracklike_) {
    muon_candidate_counter++;
}

//  std::cout<<" Selection E "<< std::endl;


 //std::cout<<"CandidateMuonIndex is "<< CandidateMuonIndex << std::endl;

 //std::cout<<"This is fine 8" << std::endl;
 //


if(CandidateMuonIndex != -1){

//std::cout<<"CandidateMuonIndex is "<< CandidateMuonIndex << std::endl;
//std::cout<<"This is 8.1 "<< std::endl;

muon_in_gap = ((Event->pfp_hitsU_->at(CandidateMuonIndex)>=1) && (Event->pfp_hitsV_->at(CandidateMuonIndex)>=1)&& (Event->pfp_hitsY_->at(CandidateMuonIndex)>=1));

//std::cout<<"This is muon in gap "<< muon_in_gap <<std::endl;
//muon momentum calc
//
// std::cout<<" Selection F "<< std::endl;


//candidate_muon_mom_mcs = Event->track_mcs_mom_mu_->at(CandidateMuonIndex); //GeV

//std::cout<<"This is 8.2 "<< std::endl;

//std::cout<<"This is muon index "<< CandidateMuonIndex <<std::endl;

//std::cout<<"This is the size "<< Event->track_mcs_mom_mu_->size() <<std::endl;

/*if (CandidateMuonIndex >= 0 &&
    CandidateMuonIndex < Event->track_mcs_mom_mu_->size()) {

    candidate_muon_mom_mcs =
        Event->track_mcs_mom_mu_->at(CandidateMuonIndex);
}
else {
    // handle missing MCS momentum
    candidate_muon_mom_mcs = -999;
}*/

/*double CandidateMuonPx = Event->pfp_true_px_->at(CandidateMuonIndex);

double CandidateMuonPy = Event->pfp_true_py_->at(CandidateMuonIndex);
double CandidateMuonPz = Event->pfp_true_pz_->at(CandidateMuonIndex);
TVector3 BackTrackCandidateMuonP(CandidateMuonPx, CandidateMuonPy, CandidateMuonPz );
candidate_muon_mom_true = BackTrackCandidateMuonP.Mag(); // GeV
*/
//std::cout<<"This is 8.3 "<< std::endl;
// For contained muons
//
 double MuonTrackEndX = Event->track_endx_->at(CandidateMuonIndex);
 double MuonTrackEndY = Event->track_endy_->at(CandidateMuonIndex);
 double MuonTrackEndZ = Event->track_endz_->at(CandidateMuonIndex);

//std::cout<<"This is 8.4 "<< std::endl;

 CandidateMuonTrackEndContainment = point_inside_FV( this->containment_FV(),
        MuonTrackEndX, MuonTrackEndY, MuonTrackEndZ );

 //if(CandidateMuonTrackEndContainment){

//std::cout<<"This is 8.5 "<<std::endl;

//candidate_muon_mom_range = Event->track_range_mom_mu_->at(CandidateMuonIndex); 

//}
if (CandidateMuonTrackEndContainment) {

    candidate_muon_mom_range =
        Event->track_range_mom_mu_->at(CandidateMuonIndex);

    candidate_muon_mom_reco =
        candidate_muon_mom_range;

}
else {

    candidate_muon_mom_reco =
        Event->track_mcs_mom_mu_->at(CandidateMuonIndex);

}

//std::cout<<"This is 8.6 "<< std::endl;I
//
// std::cout<<" Selection G "<< std::endl;


}
// No early "return false" when there is no muon candidate: fall through so the
// per-cut cut-flow histograms below are filled for these events too. They still
// fail pass_tracklike (sel_muoncandidate_tracklike_ is false here), so Passed is
// false and the final return is unchanged.  Previously this early return removed
// such events from EVERY cut-flow stage, including vertex/topology, biasing the
// per-stage efficiencies.

//std::cout<<"Muon in gap "<< muon_in_gap << std::endl;
 
 //std::cout<<"This is fine 9" << std::endl;

for (size_t i_pfp_2 = 0; i_pfp_2 < Event->track_length_->size(); i_pfp_2++) {

 //std::cout<<"This is fine 9" << std::endl;


        TVector3 track_start_ib (Event->track_startx_->at(i_pfp_2),Event->track_starty_->at(i_pfp_2),Event->track_startz_->at(i_pfp_2));
              TVector3 nu_to_track_dist_ib (track_start_ib - reco_primary_vtx);

	 double PionTrackEndX = Event->track_endx_->at(i_pfp_2);
         double PionTrackEndY = Event->track_endy_->at(i_pfp_2);
         double PionTrackEndZ = Event->track_endz_->at(i_pfp_2);

         bool CandidatePionTrackEndContainment_1 = point_inside_FV( this->containment_FV(),
         PionTrackEndX, PionTrackEndY, PionTrackEndZ );


	 if((Event->pfp_generation_->at(i_pfp_2) ==2)&&(Event->pfp_track_score_->at(i_pfp_2) >= 0.5)
              && (Event->track_llr_pid_score_->at(i_pfp_2) > -1) && (Event->track_llr_pid_score_->at(i_pfp_2) < 2)
              && (Event->trk_bragg_p_v->at(i_pfp_2)  > 0) && (Event->trk_bragg_p_v->at(i_pfp_2)  < 500)
               && (Event->trk_bragg_mu_v->at(i_pfp_2) > 0) && (Event->trk_bragg_mu_v->at(i_pfp_2) < 500)
               && (Event->trk_bragg_mip_v->at(i_pfp_2)> 0) && (Event->trk_bragg_mip_v->at(i_pfp_2)< 500)
               && (Event->track_length_->at(i_pfp_2) < 1e6) && (Event->track_length_->at(i_pfp_2) > 0)&& (i_pfp_2 != CandidateMuonIndex)
               ){

              trk_bragg_p_v_tmva = Event->trk_bragg_p_v->at(i_pfp_2);
              trk_bragg_mu_v_tmva = Event->trk_bragg_mu_v->at(i_pfp_2);
              trk_bragg_mip_v_tmva = Event->trk_bragg_mip_v->at(i_pfp_2);
              trk_llr_pid_score_v_tmva = Event->track_llr_pid_score_->at(i_pfp_2);
              trk_score_v_tmva = Event->pfp_track_score_->at(i_pfp_2);
              trk_len_v_tmva = Event->track_length_->at(i_pfp_2);
              trk_sce_end_x_v_tmva = Event->track_endx_->at(i_pfp_2);
              trk_sce_end_y_v_tmva = Event->track_endy_->at(i_pfp_2);
              trk_sce_end_z_v_tmva = Event->track_endz_->at(i_pfp_2);
              tmvaOutput = tmvaReader->EvaluateMVA("BDT");


	      trk_bragg_p_v_tmva_pi = Event->trk_bragg_p_v->at(i_pfp_2);
              trk_bragg_mu_v_tmva_pi = Event->trk_bragg_mu_v->at(i_pfp_2);
              trk_bragg_mip_v_tmva_pi = Event->trk_bragg_mip_v->at(i_pfp_2);
              trk_llr_pid_score_v_tmva_pi = Event->track_llr_pid_score_->at(i_pfp_2);
              trk_score_v_tmva_pi = Event->pfp_track_score_->at(i_pfp_2);
              trk_len_v_tmva_pi = Event->track_length_->at(i_pfp_2);
              trk_sce_end_x_v_tmva_pi = Event->track_endx_->at(i_pfp_2);
              trk_sce_end_y_v_tmva_pi = Event->track_endy_->at(i_pfp_2);
              trk_sce_end_z_v_tmva_pi = Event->track_endz_->at(i_pfp_2);
              tmvaOutput_pi= tmvaReader_pi->EvaluateMVA("BDT");

	      if((!CandidatePionTrackEndContainment_1)) continue;

	// Proton rejection on the pion candidate (matches custom selection):
	// require trk_bragg_pion >= 0.08. Fail-open (1.0) if the branch is absent.
	double bragg_pi = 1.0;
	if ( i_pfp_2 < Event->trk_bragg_pion_v->size() )
	  bragg_pi = Event->trk_bragg_pion_v->at(i_pfp_2);

	if((Event->track_llr_pid_score_->at(i_pfp_2) > 0.1) && (nu_to_track_dist_ib.Mag() <= 4) && ( Event->track_length_->at(i_pfp_2) > 20)
          && (tmvaOutput    > -0.1)
          && (tmvaOutput_pi > -0.1)
          && (bragg_pi >= 0.08)){
                     // && (tmvaOutput    > -0.1)
  

	//	if( Event->pfp_reco_pdg_->at(i_pfp_2) == PI_PLUS){

	//		sel_pioncandidate_tracklike_ = true;

	//	}
		pion_number++;
		// Pick the HIGHEST-LLR qualifying pion candidate (matches the custom
		// selection ccpi_selection.C). The previous code assigned
		// CandidatePionIndex=i_pfp_2 unconditionally and then compared the track
		// to itself, a no-op that effectively kept the LAST qualifying track.
		if(CandidatePionIndex == -1
		   || (Event->track_llr_pid_score_->at(i_pfp_2) > Event->track_llr_pid_score_->at(CandidatePionIndex))){
			CandidatePionIndex = i_pfp_2;
		}



	}
	}

}


// std::cout<<"This is fine 9" << std::endl;

if(CandidatePionIndex != -1){

	pion_in_gap = ((Event->pfp_hitsU_->at(CandidatePionIndex)>0) && (Event->pfp_hitsV_->at(CandidatePionIndex)>0)&& 
	(Event->pfp_hitsY_->at(CandidatePionIndex)>0));
	    //TVector3 pion_endpoint(Event->track_endx_->at(CandidatePionIndex),track_endy_->at(CandidatePionIndex),track_endz_->at(CandidatePionIndex));
	double PionTrackEndX = Event->track_endx_->at(CandidatePionIndex);
        double PionTrackEndY = Event->track_endy_->at(CandidatePionIndex);
        double PionTrackEndZ = Event->track_endz_->at(CandidatePionIndex);    
       
        bool CandidatePionTrackEndContainment = point_inside_FV( this->containment_FV(),
        PionTrackEndX, PionTrackEndY, PionTrackEndZ );

	if((CandidatePionTrackEndContainment) && (pion_number ==1)){
		sel_pion_contained = true;
	}
}

if((CandidateMuonIndex!= -1) && (CandidatePionIndex!= -1)){
	    TVector3 MU(Event->track_dirx_->at(CandidateMuonIndex),Event->track_diry_->at(CandidateMuonIndex),Event->track_dirz_->at(CandidateMuonIndex));
	    TVector3 PI(Event->track_dirx_->at(CandidatePionIndex),Event->track_diry_->at(CandidatePionIndex),Event->track_dirz_->at(CandidatePionIndex));
	    mu_pi_opening_angle = MU.Angle(PI);
	    if(mu_pi_opening_angle<2.6) opening_angle_cut = true;
	  }

//Shower cut
//
for (size_t i = 0; i < Event->track_length_->size(); i++) {


if(Event->pfp_generation_->at(i) !=2) continue;

if(Event->pfp_track_score_->at(i) >= 0.5){
              nPrimaryTracks++;

if(Event->track_llr_pid_score_->at(i) > 0.1  )  nonproton++;

}

else{
	      nPrimaryShowers++;
	      shower_index = i;
}
}//for i

//std::cout<<"This is the number of nonprotons " <<nonproton <<std::endl;
//std::cout<<"This is the number of primary tracks "<<nPrimaryTracks <<std::endl;

if(nPrimaryShowers == 1){


TVector3 shower_start (Event->track_startx_->at(shower_index),Event->track_starty_->at(shower_index),Event->track_startz_->at(shower_index));
TVector3 nu_to_shower_dist (shower_start - reco_primary_vtx);

if((Event->pfp_hitsY_->at(shower_index)<50) && (Event->pfp_hitsY_->at(shower_index)>0) && (nu_to_shower_dist.Mag()>10)){

	      shower_cut = true ;
	    }

}

else if (nPrimaryShowers == 0)  shower_cut = true ;
} //gt2 if 

    /*float first_pid_score = Event->track_llr_pid_score_
      ->at( CandidateIndex.at(0) );
    float second_pid_score = Event->track_llr_pid_score_
      ->at( CandidateIndex.at(1) );

    if (first_pid_score > second_pid_score && first_pid_score > 0.8) {
      CandidateMuonIndex = CandidateIndex.at(0);
      CandidateProtonIndex = CandidateIndex.at(1);
    } else {
      CandidateMuonIndex = CandidateIndex.at(1);
      CandidateProtonIndex = CandidateIndex.at(0);
    }

    if ( Event->pfp_reco_pdg_->at(CandidateMuonIndex) == MUON ) {
      sel_muoncandidate_tracklike_ = true;
    }
    if ( Event->pfp_reco_pdg_->at(CandidateProtonIndex) == MUON ) {
      sel_protoncandidate_tracklike_ = true;
    }*/
  



// ======================
  // Neutrino vertex in FV?
//  std::cout<<"This is fine 10" << std::endl; 

 sel_nuvertex_contained_ = in_vertex_FV( this->reco_FV(),
    Event->nu_vx_,Event->nu_vy_,Event->nu_vz_ );

 sel_MCVertexInFV = in_vertex_FV( this->true_FV(),
  Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );

// std::cout<<" Selection H "<< std::endl;


  // =====/===================================
  // Containment check on the muon and proton
  //
  //
  bool pass_vertex        = sel_nuvertex_contained_;
bool pass_topology      = sel_topo_cut_passed_;
bool pass_tracklike     = sel_muoncandidate_tracklike_;
bool pass_pioncontained = sel_pion_contained;
bool pass_muongap       = muon_in_gap;
bool pass_piongap       = pion_in_gap;
bool pass_shower        = shower_cut;
bool pass_opening       = opening_angle_cut;
bool pass_final         = (nonproton < 4 && nPrimaryTracks < 5);
// Software-trigger requirement (matches custom selection). MC overlay is not
// pre-filtered on swtrig, so without this the MC is over-normalised relative
// to data/EXT. Defaults to pass when the branch is absent. This affects only
// the reco selection (Passed) and NOT the efficiency denominator filled in
// define_signal(), so swtrig acceptance is correctly folded into efficiency.
bool pass_swtrig        = (Event->swtrig_ == 1);

//bool Passed = sel_nuvertex_contained_ && sel_topo_cut_passed_  && sel_muoncandidate_tracklike_  && sel_pion_contained && muon_in_gap && pion_in_gap
 // && shower_cut && opening_angle_cut&& (nonproton < 4 && nPrimaryTracks <5); //&& sel_nuvertex_contained_ ;//&& sel_MCVertexInFV && sel_nslice_eq_1_ && sel_nuvertex_contained_ && sel_topo_cut_passed_ && sel_muoncandidate_tracklike_&& muon_in_gap && sel_pion_contained  && pion_in_gap && shower_cut && opening_angle_cut && (nonproton < 4 && nPrimaryTracks <5);

bool Passed =
    pass_swtrig &&
    pass_vertex &&
    pass_topology &&
    pass_tracklike &&
    pass_pioncontained &&
    pass_muongap &&
    pass_piongap &&
    pass_shower &&
    pass_opening &&
    pass_final;

  // --- background-control sidebands (signal-depleted; see CC1mu1piXp.hh) ---
  // Common preselection: software trigger, contained vertex, topology, muon track.
  bool sb_pre = pass_swtrig && pass_vertex && pass_topology && pass_tracklike;
  //  CC0pi: muon selection, zero charged-pion candidates, shower veto passed.
  sb_cc0pi_   = sb_pre && pass_shower && pass_final && ( pion_number == 0 );
  //  Multi-pi: two or more charged-pion candidates.
  sb_multipi_ = sb_pre && pass_final && ( pion_number >= 2 );
  //  Pi0/shower: a reco shower is present and the shower veto FAILED (pi0-like).
  sb_pi0_     = sb_pre && !pass_shower && ( nPrimaryShowers >= 1 );
  //  Cosmic: the full signal selection but with the opening angle inverted
  //  (theta_mupi > 2.6 rad), the EXT-cosmic--dominated region.
  sb_cosmic_  = sb_pre && pass_pioncontained && pass_muongap && pass_piongap
                && pass_shower && pass_final && !pass_opening;

  // std::cout<<"This is Passed " << Passed <<std::endl;
  //
  //std::cout<<"Exiting selection" <<std::endl;
  //
  //
  //
  //
  if (this->is_event_mc_signal() && truemuonindex != -1) {

    double muon_true = TrueCandidateMuonP.Mag();
    double pion_true = TrueCandidatePionP.Mag();

    // Cut 0 (all signal events) is the efficiency denominator and is filled
    // once per signal event in define_signal() (h_mu_cut0 / h_pi_cut0 /
    // h_eff_den).  It is intentionally NOT re-filled here: doing so double-
    // counted the denominator and roughly halved every per-cut efficiency.

    // ------------------------------
    // Cut 1 = vertex containment
    // ------------------------------
    if (pass_vertex) {

        h_mu_cut1_vertex->Fill(muon_true, evt_w);
        h_pi_cut1_vertex->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 2 = topology
    // ------------------------------
    if (pass_vertex &&
        pass_topology) {

        h_mu_cut2_topology->Fill(muon_true, evt_w);
        h_pi_cut2_topology->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 3 = muon track-like
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike) {

        h_mu_cut3_tracklike->Fill(muon_true, evt_w);
        h_pi_cut3_tracklike->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 4 = pion contained
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained) {

        h_mu_cut4_pioncontained->Fill(muon_true, evt_w);
        h_pi_cut4_pioncontained->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 5 = muon gap
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained &&
        pass_muongap) {

        h_mu_cut5_muongap->Fill(muon_true, evt_w);
        h_pi_cut5_muongap->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 6 = pion gap
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained &&
        pass_muongap &&
        pass_piongap) {

        h_mu_cut6_piongap->Fill(muon_true, evt_w);
        h_pi_cut6_piongap->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 7 = shower cut
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained &&
        pass_muongap &&
        pass_piongap &&
        pass_shower) {

        h_mu_cut7_shower->Fill(muon_true, evt_w);
        h_pi_cut7_shower->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 8 = opening angle cut
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained &&
        pass_muongap &&
        pass_piongap &&
        pass_shower &&
        pass_opening) {

        h_mu_cut8_opening->Fill(muon_true, evt_w);
        h_pi_cut8_opening->Fill(pion_true, evt_w);
    }

    // ------------------------------
    // Cut 9 = final multiplicity cut
    // ------------------------------
    if (pass_vertex &&
        pass_topology &&
        pass_tracklike &&
        pass_pioncontained &&
        pass_muongap &&
        pass_piongap &&
        pass_shower &&
        pass_opening &&
        pass_final) {

        h_mu_cut9_final->Fill(muon_true, evt_w);
        h_pi_cut9_final->Fill(pion_true, evt_w);
    }
}

  // ---- per-observable RECO cut-flow (ALL events; mirrors the cumulative cuts) ----
  // Guard against the non-finite (inf/nan) event weights present in a small number
  // of overlay events (they otherwise poison the histogram integral).
  if ( std::isfinite(evt_w) ) {
    bool mu_ok = (CandidateMuonIndex != -1);
    bool pi_ok = (CandidatePionIndex != -1);
    double cf_rv[5] = { candidate_muon_mom_reco, candidate_pion_mom_reco,
                        candidate_muon_costh_reco, candidate_pion_costh_reco,
                        mu_pi_opening_angle };
    bool cf_ok[5] = { mu_ok, pi_ok, mu_ok, pi_ok, (mu_ok && pi_ok) };
    bool cf_pass[10];
    cf_pass[0] = true;
    cf_pass[1] = cf_pass[0] && pass_vertex;
    cf_pass[2] = cf_pass[1] && pass_topology;
    cf_pass[3] = cf_pass[2] && pass_tracklike;
    cf_pass[4] = cf_pass[3] && pass_pioncontained;
    cf_pass[5] = cf_pass[4] && pass_muongap;
    cf_pass[6] = cf_pass[5] && pass_piongap;
    cf_pass[7] = cf_pass[6] && pass_shower;
    cf_pass[8] = cf_pass[7] && pass_opening;
    cf_pass[9] = cf_pass[8] && pass_final;
    bool cf_is_sig = this->is_event_mc_signal();
    for (int c = 0; c < 10; ++c) {
      if (!cf_pass[c]) continue;
      // cut-flow yields (all events; signal split for the bkg = tot - sig)
      h_cutflow_tot->Fill(c, evt_w);
      if (cf_is_sig) h_cutflow_sig->Fill(c, evt_w);
      for (int o = 0; o < 5; ++o)
        if (cf_ok[o] && std::isfinite(cf_rv[o])) h_cf[o][c]->Fill(cf_rv[o], evt_w);
    }

    // ---- selection diagnostics (N-1 event-level, final-cut candidate dists, bkg cat) ----
    int cf_k = cf_is_sig ? 0 : 1;
    auto mu_llr = [&](){ return (mu_ok && CandidateMuonIndex < (int)Event->track_llr_pid_score_->size())
      ? Event->track_llr_pid_score_->at(CandidateMuonIndex) : -2.; };
    auto mu_len = [&](){ return (mu_ok && CandidateMuonIndex < (int)Event->track_length_->size())
      ? Event->track_length_->at(CandidateMuonIndex) : -1.; };
    auto pi_llr = [&](){ return (pi_ok && CandidatePionIndex < (int)Event->track_llr_pid_score_->size())
      ? Event->track_llr_pid_score_->at(CandidatePionIndex) : -2.; };
    auto pi_len = [&](){ return (pi_ok && CandidatePionIndex < (int)Event->track_length_->size())
      ? Event->track_length_->at(CandidatePionIndex) : -1.; };
    // N-1 topology: all cuts except topology (candidates present)
    if (mu_ok && pi_ok && pass_vertex && pass_tracklike && pass_pioncontained
        && pass_muongap && pass_piongap && pass_shower && pass_opening && pass_final)
      h_nm1_topo[cf_k]->Fill(Event->topological_score_, evt_w);
    // N-1 opening angle: all cuts except opening
    if (mu_ok && pi_ok && pass_vertex && pass_topology && pass_tracklike && pass_pioncontained
        && pass_muongap && pass_piongap && pass_shower && pass_final)
      h_nm1_oa[cf_k]->Fill(mu_pi_opening_angle, evt_w);
    // final-cut candidate distributions + background decomposition
    if (cf_pass[9]) {
      if (mu_ok) { h_fin_mupid[cf_k]->Fill(mu_llr(), evt_w); h_fin_mulen[cf_k]->Fill(mu_len(), evt_w); }
      if (pi_ok) { h_fin_pipid[cf_k]->Fill(pi_llr(), evt_w); h_fin_pilen[cf_k]->Fill(pi_len(), evt_w); }
      h_bkgcat->Fill(this->categorize_event(Event), evt_w);
    }
  }


  if(CandidateMuonIndex != -1){
 

    double muon_mcs = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);

    if (Passed) {
        //h_selected->Fill(muon_mcs); CHANGED HERE
        h_selected->Fill(candidate_muon_mom_reco, evt_w);
        h_all_selected->Fill(candidate_muon_mom_reco, evt_w);
       // h_all_selected->Fill(muon_mcs);
    }

    if (Passed && !this->is_event_mc_signal()) {
      //  h_background->Fill(muon_mcs);
	//h_sel_bkg->Fill(muon_mcs);
	h_background->Fill(candidate_muon_mom_reco, evt_w);
	h_sel_bkg->Fill(candidate_muon_mom_reco, evt_w);
    }


if (Passed && this->is_event_mc_signal()) {

    double muon_true = TrueCandidateMuonP.Mag();

    h_eff_num->Fill(muon_true, evt_w);
    h_sel_signal->Fill(candidate_muon_mom_reco, evt_w);
    //h_sel_signal->Fill(muon_mcs); CHANGED HERE


    //double muon_true = TrueCandidateMuonP.Mag();
    double pion_true = TrueCandidatePionP.Mag();

    h_mu_eff_num->Fill(muon_true, evt_w);
    h_pi_eff_num->Fill(pion_true, evt_w);
}
} 

 if (Passed) {
//    std::cout << "PASSED\n";
} else {
//    std::cout << "FAILED: "
//              << sel_nuvertex_contained_ << " "
//              << sel_topo_cut_passed_ << " "
//              << sel_muoncandidate_tracklike_ << " "
//              << sel_pion_contained << " "
//              << muon_in_gap << " "
//              << pion_in_gap << " "
//              << shower_cut << " "
//              << opening_angle_cut << " "
//              << (nonproton < 4 && nPrimaryTracks < 5)
//
//              << std::endl;
}

//std::cout << sel_nuvertex_contained_ << " "
//          << sel_topo_cut_passed_ << " "
//          << sel_muoncandidate_tracklike_ << " "
//           << std::endl;

//std::cout<< "Fid vol counter " << sig_truevertex_fv << std::endl;
//std::cout<< "sig ccnc counter " << sig_ccnc << std::endl;
//std::cout<< "sig numu counter " << sig_numu << std::endl;
//std::cout<< "Sig one muon  counter " << sig_one_muon << std::endl;

//std::cout<< " muon candidate counter from selection " <<  muon_candidate_counter << std::endl;

//std::cout<<"Passed the variable " << Passed <<std::endl;
//   return Passed;
  // return false;

   //std::cout<<"This is Passed " << Passed <<std::endl;


// std::cout<<"Exiting selection" <<std::endl; 
// std::cout<<" Selection I "<< std::endl;

if ( Passed && this->is_event_mc_signal()) {

    double muon_true = TrueCandidateMuonP.Mag();

    if (CandidateMuonIndex != -1) {

        double muon_reco = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);

        // migration (truth → reco)
       // h_response_full->Fill(muon_reco, muon_true);
//	h_response_full->Fill(muon_true, muon_true);
        h_response_full->Fill(candidate_muon_mom_reco,muon_true, evt_w);
    }

    // else: inefficiency (true but not reconstructed)
    // do nothing here (handled later in unfolding)
}
 return Passed;

 }//selection


void CC1mu1piXp::finalize() {

    // NOTE: the absolute cross section (flux x N_targets normalisation, bin-width
    // division, unfolding) is produced by the dedicated Unfolder / UnfolderNuMI
    // stage via CrossSectionExtractor.  This method only writes selection
    // diagnostics: reco spectra, the efficiency/purity and per-cut histograms,
    // and the response matrix.  The former in-method "quick" cross section used a
    // hardcoded all-run flux + full-FV target count that did not match the actual
    // (Run1-only) exposure, so it was removed to avoid a misleading second answer.

    // --------------------------------------------------
// POT NORMALIZATION
// --------------------------------------------------

/*double mcPOT   = 139.69e20;
double dataPOT = 8.86e20;

double potScale = dataPOT / mcPOT;

std::cout << "MC POT     = " << mcPOT << std::endl;
std::cout << "Data POT   = " << dataPOT << std::endl;
std::cout << "POT Scale  = " << potScale << std::endl;

// Scale ALL MC histograms

h_background->Scale(potScale);

h_sel_signal->Scale(potScale);
h_all_signal->Scale(potScale);

h_eff_num->Scale(potScale);
h_eff_den->Scale(potScale);

h_mu_eff->Divide(h_mu_eff_num, h_mu_eff_den);

h_pi_eff->Divide(h_pi_eff_num, h_pi_eff_den);



h_response_full->Scale(potScale);*/

    // =========================
    // 2. Build signal histogram
    // =========================
    //
    //
    TH1D* eff_vertex =
    (TH1D*)h_mu_cut1_vertex->Clone("eff_vertex");
eff_vertex->Divide(h_mu_cut0);

TH1D* eff_topology =
    (TH1D*)h_mu_cut2_topology->Clone("eff_topology");
eff_topology->Divide(h_mu_cut0);

TH1D* eff_tracklike =
    (TH1D*)h_mu_cut3_tracklike->Clone("eff_tracklike");
eff_tracklike->Divide(h_mu_cut0);

TH1D* eff_pioncontained =
    (TH1D*)h_mu_cut4_pioncontained->Clone("eff_pioncontained");
eff_pioncontained->Divide(h_mu_cut0);

TH1D* eff_muongap =
    (TH1D*)h_mu_cut5_muongap->Clone("eff_muongap");
eff_muongap->Divide(h_mu_cut0);

TH1D* eff_piongap =
    (TH1D*)h_mu_cut6_piongap->Clone("eff_piongap");
eff_piongap->Divide(h_mu_cut0);

TH1D* eff_shower =
    (TH1D*)h_mu_cut7_shower->Clone("eff_shower");
eff_shower->Divide(h_mu_cut0);

TH1D* eff_opening =
    (TH1D*)h_mu_cut8_opening->Clone("eff_opening");
eff_opening->Divide(h_mu_cut0);

TH1D* eff_final =
    (TH1D*)h_mu_cut9_final->Clone("eff_final");
eff_final->Divide(h_mu_cut0);

h_pi_cut0->Write();
h_pi_cut1_vertex->Write();
h_pi_cut2_topology->Write();
h_pi_cut3_tracklike->Write();
h_pi_cut4_pioncontained->Write();
h_pi_cut5_muongap->Write();
h_pi_cut6_piongap->Write();
h_pi_cut7_shower->Write();
h_pi_cut8_opening->Write();
h_pi_cut9_final->Write();

h_mu_cut0->Write();
h_mu_cut1_vertex->Write();
h_mu_cut2_topology->Write();
h_mu_cut3_tracklike->Write();
h_mu_cut4_pioncontained->Write();
h_mu_cut5_muongap->Write();
h_mu_cut6_piongap->Write();
h_mu_cut7_shower->Write();
h_mu_cut8_opening->Write();
h_mu_cut9_final->Write();

// Per-cut RECO cut-flow histograms (5 obs x 10 stages)
for (int o = 0; o < 5; ++o)
  for (int c = 0; c < 10; ++c)
    if (h_cf[o][c]) h_cf[o][c]->Write();

// Cut-flow yield counters (total + signal per stage)
if (h_cutflow_tot) h_cutflow_tot->Write();
if (h_cutflow_sig) h_cutflow_sig->Write();

// Selection-diagnostic histograms
for (int k = 0; k < 2; ++k) {
  if (h_nm1_topo[k])  h_nm1_topo[k]->Write();
  if (h_nm1_oa[k])    h_nm1_oa[k]->Write();
  if (h_fin_mupid[k]) h_fin_mupid[k]->Write();
  if (h_fin_pipid[k]) h_fin_pipid[k]->Write();
  if (h_fin_mulen[k]) h_fin_mulen[k]->Write();
  if (h_fin_pilen[k]) h_fin_pilen[k]->Write();
}
if (h_bkgcat) h_bkgcat->Write();

eff_vertex->Write();
eff_topology->Write();
eff_tracklike->Write();
eff_pioncontained->Write();
eff_muongap->Write();
eff_piongap->Write();
eff_shower->Write();
eff_opening->Write();
eff_final->Write();

TH1D* eff_mu_vertex =
    (TH1D*)h_mu_cut1_vertex->Clone("eff_mu_vertex");
eff_mu_vertex->Divide(h_mu_cut0);

TH1D* eff_mu_topology =
    (TH1D*)h_mu_cut2_topology->Clone("eff_mu_topology");
eff_mu_topology->Divide(h_mu_cut1_vertex);

TH1D* eff_mu_tracklike =
    (TH1D*)h_mu_cut3_tracklike->Clone("eff_mu_tracklike");
eff_mu_tracklike->Divide(h_mu_cut2_topology);

TH1D* eff_mu_pioncontained =
    (TH1D*)h_mu_cut4_pioncontained->Clone("eff_mu_pioncontained");
eff_mu_pioncontained->Divide(h_mu_cut3_tracklike);

TH1D* eff_mu_muongap =
    (TH1D*)h_mu_cut5_muongap->Clone("eff_mu_muongap");
eff_mu_muongap->Divide(h_mu_cut4_pioncontained);

TH1D* eff_mu_piongap =
    (TH1D*)h_mu_cut6_piongap->Clone("eff_mu_piongap");
eff_mu_piongap->Divide(h_mu_cut5_muongap);

TH1D* eff_mu_shower =
    (TH1D*)h_mu_cut7_shower->Clone("eff_mu_shower");
eff_mu_shower->Divide(h_mu_cut6_piongap);

TH1D* eff_mu_opening =
    (TH1D*)h_mu_cut8_opening->Clone("eff_mu_opening");
eff_mu_opening->Divide(h_mu_cut7_shower);

TH1D* eff_mu_final =
    (TH1D*)h_mu_cut9_final->Clone("eff_mu_final");
eff_mu_final->Divide(h_mu_cut8_opening);


TH1D* eff_pi_vertex =
    (TH1D*)h_pi_cut1_vertex->Clone("eff_pi_vertex");
eff_pi_vertex->Divide(h_pi_cut0);

TH1D* eff_pi_topology =
    (TH1D*)h_pi_cut2_topology->Clone("eff_pi_topology");
eff_pi_topology->Divide(h_pi_cut0);

TH1D* eff_pi_tracklike =
    (TH1D*)h_pi_cut3_tracklike->Clone("eff_pi_tracklike");
eff_pi_tracklike->Divide(h_pi_cut0);

TH1D* eff_pi_pioncontained =
    (TH1D*)h_pi_cut4_pioncontained->Clone("eff_pi_pioncontained");
eff_pi_pioncontained->Divide(h_pi_cut0);

TH1D* eff_pi_muongap =
    (TH1D*)h_pi_cut5_muongap->Clone("eff_pi_muongap");
eff_pi_muongap->Divide(h_pi_cut0);

TH1D* eff_pi_piongap =
    (TH1D*)h_pi_cut6_piongap->Clone("eff_pi_piongap");
eff_pi_piongap->Divide(h_pi_cut0);

TH1D* eff_pi_shower =
    (TH1D*)h_pi_cut7_shower->Clone("eff_pi_shower");
eff_pi_shower->Divide(h_pi_cut0);

TH1D* eff_pi_opening =
    (TH1D*)h_pi_cut8_opening->Clone("eff_pi_opening");
eff_pi_opening->Divide(h_pi_cut0);

TH1D* eff_pi_final =
    (TH1D*)h_pi_cut9_final->Clone("eff_pi_final");
eff_pi_final->Divide(h_pi_cut0);

TH1D* eff_mu_vertex_only =
    (TH1D*)h_mu_cut1_vertex->Clone("eff_mu_vertex_only");
eff_mu_vertex_only->Divide(h_mu_cut0);

TH1D* eff_mu_topology_only =
    (TH1D*)h_mu_cut2_topology->Clone("eff_mu_topology_only");
eff_mu_topology_only->Divide(h_mu_cut1_vertex);

TH1D* eff_mu_tracklike_only =
    (TH1D*)h_mu_cut3_tracklike->Clone("eff_mu_tracklike_only");
eff_mu_tracklike_only->Divide(h_mu_cut2_topology);

TH1D* eff_mu_pioncontained_only =
    (TH1D*)h_mu_cut4_pioncontained->Clone("eff_mu_pioncontained_only");
eff_mu_pioncontained_only->Divide(h_mu_cut3_tracklike);

TH1D* eff_mu_muongap_only =
    (TH1D*)h_mu_cut5_muongap->Clone("eff_mu_muongap_only");
eff_mu_muongap_only->Divide(h_mu_cut4_pioncontained);

TH1D* eff_mu_piongap_only =
    (TH1D*)h_mu_cut6_piongap->Clone("eff_mu_piongap_only");
eff_mu_piongap_only->Divide(h_mu_cut5_muongap);

TH1D* eff_mu_shower_only =
    (TH1D*)h_mu_cut7_shower->Clone("eff_mu_shower_only");
eff_mu_shower_only->Divide(h_mu_cut6_piongap);

TH1D* eff_mu_opening_only =
    (TH1D*)h_mu_cut8_opening->Clone("eff_mu_opening_only");
eff_mu_opening_only->Divide(h_mu_cut7_shower);

TH1D* eff_mu_final_only =
    (TH1D*)h_mu_cut9_final->Clone("eff_mu_final_only");
eff_mu_final_only->Divide(h_mu_cut8_opening);

TH1D* eff_pi_vertex_only =
    (TH1D*)h_pi_cut1_vertex->Clone("eff_pi_vertex_only");
eff_pi_vertex_only->Divide(h_pi_cut0);

TH1D* eff_pi_topology_only =
    (TH1D*)h_pi_cut2_topology->Clone("eff_pi_topology_only");
eff_pi_topology_only->Divide(h_pi_cut1_vertex);

TH1D* eff_pi_tracklike_only =
    (TH1D*)h_pi_cut3_tracklike->Clone("eff_pi_tracklike_only");
eff_pi_tracklike_only->Divide(h_pi_cut2_topology);

TH1D* eff_pi_pioncontained_only =
    (TH1D*)h_pi_cut4_pioncontained->Clone("eff_pi_pioncontained_only");
eff_pi_pioncontained_only->Divide(h_pi_cut3_tracklike);

TH1D* eff_pi_muongap_only =
    (TH1D*)h_pi_cut5_muongap->Clone("eff_pi_muongap_only");
eff_pi_muongap_only->Divide(h_pi_cut4_pioncontained);

TH1D* eff_pi_piongap_only =
    (TH1D*)h_pi_cut6_piongap->Clone("eff_pi_piongap_only");
eff_pi_piongap_only->Divide(h_pi_cut5_muongap);

TH1D* eff_pi_shower_only =
    (TH1D*)h_pi_cut7_shower->Clone("eff_pi_shower_only");
eff_pi_shower_only->Divide(h_pi_cut6_piongap);

TH1D* eff_pi_opening_only =
    (TH1D*)h_pi_cut8_opening->Clone("eff_pi_opening_only");
eff_pi_opening_only->Divide(h_pi_cut7_shower);

TH1D* eff_pi_final_only =
    (TH1D*)h_pi_cut9_final->Clone("eff_pi_final_only");
eff_pi_final_only->Divide(h_pi_cut8_opening);

eff_pi_vertex->Write();
eff_pi_topology->Write();
eff_pi_tracklike->Write();
eff_pi_pioncontained->Write();
eff_pi_muongap->Write();
eff_pi_piongap->Write();
eff_pi_shower->Write();
eff_pi_opening->Write();
eff_pi_final->Write();

eff_mu_vertex_only->Write();
eff_mu_topology_only->Write();
eff_mu_tracklike_only->Write();
eff_mu_pioncontained_only->Write();
eff_mu_muongap_only->Write();
eff_mu_piongap_only->Write();
eff_mu_shower_only->Write();
eff_mu_opening_only->Write();
eff_mu_final_only->Write();

eff_pi_vertex_only->Write();
eff_pi_topology_only->Write();
eff_pi_tracklike_only->Write();
eff_pi_pioncontained_only->Write();
eff_pi_muongap_only->Write();
eff_pi_piongap_only->Write();
eff_pi_shower_only->Write();
eff_pi_opening_only->Write();
eff_pi_final_only->Write();




TCanvas* c_mu = new TCanvas(
    "c_mu",
    "Muon Cut Efficiency",
    900,
    700
);

eff_mu_vertex->SetLineColor(kBlack);
eff_mu_topology->SetLineColor(kRed);
eff_mu_tracklike->SetLineColor(kBlue);
eff_mu_pioncontained->SetLineColor(kGreen+2);
eff_mu_muongap->SetLineColor(kMagenta);
eff_mu_piongap->SetLineColor(kOrange+1);
eff_mu_shower->SetLineColor(kCyan+2);
eff_mu_opening->SetLineColor(kViolet);
eff_mu_final->SetLineColor(kGray+2);

eff_mu_vertex->SetMaximum(1.05);
eff_mu_vertex->SetMinimum(0.0);

eff_mu_vertex->Draw("hist");
eff_mu_topology->Draw("hist same");
eff_mu_tracklike->Draw("hist same");
eff_mu_pioncontained->Draw("hist same");
eff_mu_muongap->Draw("hist same");
eff_mu_piongap->Draw("hist same");
eff_mu_shower->Draw("hist same");
eff_mu_opening->Draw("hist same");
eff_mu_final->Draw("hist same");

TLegend* leg = new TLegend(0.55,0.45,0.88,0.88);

leg->AddEntry(eff_mu_vertex,"Vertex","l");
leg->AddEntry(eff_mu_topology,"Topology","l");
leg->AddEntry(eff_mu_tracklike,"Tracklike","l");
leg->AddEntry(eff_mu_pioncontained,"Pion contained","l");
leg->AddEntry(eff_mu_muongap,"Muon gap","l");
leg->AddEntry(eff_mu_piongap,"Pion gap","l");
leg->AddEntry(eff_mu_shower,"Shower veto","l");
leg->AddEntry(eff_mu_opening,"Opening angle","l");
leg->AddEntry(eff_mu_final,"Final multiplicity","l");

leg->Draw();

TCanvas* c_pi = new TCanvas(
    "c_pi",
    "Pion Cut Efficiency",
    900,
    700
);

eff_pi_vertex->SetLineColor(kBlack);
eff_pi_topology->SetLineColor(kRed);
eff_pi_tracklike->SetLineColor(kBlue);
eff_pi_pioncontained->SetLineColor(kGreen+2);
eff_pi_muongap->SetLineColor(kMagenta);
eff_pi_piongap->SetLineColor(kOrange+1);
eff_pi_shower->SetLineColor(kCyan+2);
eff_pi_opening->SetLineColor(kViolet);
eff_pi_final->SetLineColor(kGray+2);

eff_pi_vertex->SetLineWidth(2);
eff_pi_topology->SetLineWidth(2);
eff_pi_tracklike->SetLineWidth(2);
eff_pi_pioncontained->SetLineWidth(2);
eff_pi_muongap->SetLineWidth(2);
eff_pi_piongap->SetLineWidth(2);
eff_pi_shower->SetLineWidth(2);
eff_pi_opening->SetLineWidth(2);
eff_pi_final->SetLineWidth(2);

eff_pi_vertex->SetTitle(
    "Pion Momentum Cut Efficiencies; p_{#pi}^{true} [GeV/c]; Efficiency"
);

eff_pi_vertex->SetMaximum(1.05);
eff_pi_vertex->SetMinimum(0.0);

eff_pi_vertex->Draw("hist");
eff_pi_topology->Draw("hist same");
eff_pi_tracklike->Draw("hist same");
eff_pi_pioncontained->Draw("hist same");
eff_pi_muongap->Draw("hist same");
eff_pi_piongap->Draw("hist same");
eff_pi_shower->Draw("hist same");
eff_pi_opening->Draw("hist same");
eff_pi_final->Draw("hist same");

TLegend* leg_pi = new TLegend(
    0.55, 0.45,
    0.88, 0.88
);

leg_pi->AddEntry(eff_pi_vertex,"Vertex","l");
leg_pi->AddEntry(eff_pi_topology,"Topology","l");
leg_pi->AddEntry(eff_pi_tracklike,"Tracklike","l");
leg_pi->AddEntry(eff_pi_pioncontained,"Pion Contained","l");
leg_pi->AddEntry(eff_pi_muongap,"Muon Gap","l");
leg_pi->AddEntry(eff_pi_piongap,"Pion Gap","l");
leg_pi->AddEntry(eff_pi_shower,"Shower Cut","l");
leg_pi->AddEntry(eff_pi_opening,"Opening Angle","l");
leg_pi->AddEntry(eff_pi_final,"Multiplicity","l");

leg_pi->Draw();

c_mu->Write();
c_pi->Write();




    h_mu_eff->Divide(h_mu_eff_num, h_mu_eff_den);
    h_pi_eff->Divide(h_pi_eff_num, h_pi_eff_den);
    h_signal = (TH1D*)h_selected->Clone("h_signal");
    h_signal->Add(h_background, -1);

    // =========================
    // 3. Efficiency
    // =========================
    h_eff = (TH1D*)h_eff_num->Clone("h_eff");
    h_eff->SetTitle("Efficiency; p_{#mu}^{true} [GeV/c]; Efficiency");
    h_eff->Divide(h_eff_num, h_eff_den, 1.0, 1.0, "B"); // binomial errors

// =========================
    // 4. Purity
    // =========================
    h_purity = (TH1D*)h_sel_signal->Clone("h_purity");
    h_purity->SetTitle("Purity; p_{#mu}^{reco} [GeV/c]; Purity");
    h_purity->Divide(h_sel_signal, h_all_selected, 1.0, 1.0, "B");

   // h_eff->Divide(h_eff_den);

    // =========================
    // 7. Write everything (after all corrections!)
    // =========================
    h_selected->Write();
    h_background->Write();
    h_signal->Write();

    h_eff_num->Write();
    h_eff_den->Write();
    h_eff->Write();

    h_all_signal->Write();
    h_all_selected->Write();
    h_sel_signal->Write();
    h_sel_bkg->Write();
    h_purity->Write();
    h_mu_eff_num->Write();
    h_mu_eff_den->Write();
    h_mu_eff->Write();

    h_pi_eff_num->Write();
    h_pi_eff_den->Write();
    h_pi_eff->Write();
    gStyle->SetOptStat(0);

h_response_full->GetZaxis()->SetTitle("Events");
h_response_full->GetZaxis()->SetTitleOffset(1.2);
    h_response_full->Write();

    // =========================
    // 8. Debug prints
    // =========================
    std::cout << "==========================" << std::endl;
    std::cout << "Selected: " << h_selected->Integral() << std::endl;
    std::cout << "Background: " << h_background->Integral() << std::endl;
    std::cout << "Signal: " << h_signal->Integral() << std::endl;

    std::cout << "Eff numerator: " << h_eff_num->Integral() << std::endl;
    std::cout << "Eff denominator: " << h_eff_den->Integral() << std::endl;
    std::cout << "Efficiency (avg): " << h_eff->Integral() << std::endl;
    std::cout << "==========================" << std::endl;

std::cout << "sig_truevertex_fv: " << sig_truevertex_fv << std::endl;
std::cout << "sig_ccnc: " << sig_ccnc << std::endl;
std::cout << "sig_numu: " << sig_numu << std::endl;
std::cout << "sig_one_muon: " << sig_one_muon << std::endl;

std::cout << "All signal: " << h_all_signal->Integral() << std::endl;
std::cout << "All selected: " << h_all_selected->Integral() << std::endl;
std::cout << "Selected signal: " << h_sel_signal->Integral() << std::endl;
std::cout << "Selected bkg: " << h_sel_bkg->Integral() << std::endl;
std::cout << "Purity: " << h_purity->Integral() << std::endl;
double eff = 0.0;
if (h_all_signal->Integral() > 0) {
    eff = h_sel_signal->Integral() / h_all_signal->Integral();
}

double pur = 0.0;
if (h_all_selected->Integral() > 0) {
    pur = h_sel_signal->Integral() / h_all_selected->Integral();
}

std::cout << std::fixed << std::setprecision(2);
std::cout << "Efficiency: " << eff * 100.0 << " %" << std::endl;
std::cout << "Purity: " << pur * 100.0 << " %" << std::endl;

// =========================
// 9. Reco signal (selected - background)
// =========================
// Unfolding and the unfolded cross section are intentionally NOT done here —
// that is the job of the dedicated Unfolder stage.  The former in-selection
// RooUnfold SVD block was removed together with the RooUnfold dependency; only
// the reco-space signal histogram is written out for downstream use.
TH1D* h_data_signal = (TH1D*)h_selected->Clone("h_data_signal");
h_data_signal->Add(h_background, -1);
h_data_signal->SetTitle("Reco Signal Before Unfolding; p_{#mu}^{reco} [GeV/c]; Events");
h_data_signal->Write();
} //finalise

void CC1mu1piXp::define_output_branches() {

  set_branch( &sel_nslice_eq_1_, "nslice_eq_1" );
  set_branch( &sel_nuvertex_contained_, "nuvertex_contained_" );
  set_branch( &sel_topo_cut_passed_, "topo_cut_passed" );
  set_branch( &sig_truevertex_in_fv_, "sig_truevertex_in_fv" );
  set_branch( &sig_ccnc_, "sig_ccnc" );
  set_branch( &sig_is_numu_, "sig_is_numu" );
  set_branch( &sig_one_muon_above_thresh_, "sig_one_muon_above_thresh" );
  set_branch( &sig_one_proton_above_thresh_, "sig_one_proton_above_thresh" );
  set_branch( &sig_no_pions_, "sig_no_pions" );
  set_branch( &sig_no_heavy_mesons_, "sig_no_heavy_mesons" );
  set_branch( &sig_mc_n_threshold_muon, "mc_n_threshold_muon" );
  set_branch( &sig_mc_n_threshold_proton, "mc_n_threshold_proton" );
  set_branch( &sig_mc_n_threshold_pion0, "mc_n_threshold_pion0" );
  set_branch( &sig_mc_n_threshold_pionpm, "mc_n_threshold_pionpm" );
  set_branch( &sig_mc_n_heaviermeson, "mc_n_heaviermeson" );
  set_branch( &CandidateMuonIndex, "CandidateMuonIndex" );
  set_branch( &sig_one_charged_pion_, "sig_one_charged_pion" );
  set_branch( &CandidatePionIndex, "CandidatePionIndex" );
  set_branch( &sel_MCVertexInFV, "sel_MCVertexInFV" );
  set_branch( &sig_recovertex_in_fv_, "sig_recovertex_in_fv_" );
  // background-control sideband flags (signal-depleted regions)
  set_branch( &sb_cc0pi_,   "sb_cc0pi" );
  set_branch( &sb_multipi_, "sb_multipi" );
  set_branch( &sb_pi0_,     "sb_pi0" );
  set_branch( &sb_cosmic_,  "sb_cosmic" );
  set_branch( &candidate_muon_mom_mcs, "candidate_muon_mom_mcs");
  set_branch( &candidate_muon_mom_true, "candidate_muon_mom_true");
  set_branch( &candidate_muon_mom_range, "candidate_muon_mom_range");
  // The estimator that should actually be unfolded: range momentum when the
  // muon track end is contained, MCS otherwise. This is what the internal
  // efficiency/response histograms below are filled with, so the bin configs
  // must use it too -- binning on _mcs alone applies MCS to contained muons,
  // and binning on _range alone drops every uncontained muon (its _range stays
  // at the 0.0 reset value and falls below the first bin).
  set_branch( &candidate_muon_mom_reco, "candidate_muon_mom_reco");
  set_branch( &CandidateMuonTrackEndContainment, "CandidateMuonTrackEndContainment");
  set_branch( &sel_muoncandidate_tracklike_, "sel_muoncandidate_tracklike");
  set_branch( &mu_pi_opening_angle, "mu_pi_opening_angle");
  set_branch( &true_mu_pi_opening_angle, "true_mu_pi_opening_angle");
  set_branch( &candidate_pion_mom_reco, "candidate_pion_mom_reco");
  set_branch( &candidate_pion_mom_true, "candidate_pion_mom_true");
  set_branch( &candidate_muon_costh_reco, "candidate_muon_costh_reco");
  set_branch( &candidate_muon_costh_true, "candidate_muon_costh_true");
  set_branch( &candidate_pion_costh_reco, "candidate_pion_costh_reco");
  set_branch( &candidate_pion_costh_true, "candidate_pion_costh_true");






 // set_branch( &pion_number, "pion_number");
 

}

void CC1mu1piXp::reset() {

  sel_nslice_eq_1_ = false;
  sel_nuvertex_contained_ = false;
  sel_topo_cut_passed_ = false;
  sig_truevertex_in_fv_ = false;
  sig_ccnc_ = false;
  sig_is_numu_ = false;
  sig_one_muon_above_thresh_ = false;
  sig_one_proton_above_thresh_ = false;
  sig_no_pions_ = false;
  sig_no_heavy_mesons_ = false;
  sig_one_charged_pion_ = false;
  sig_no_kaons_ = false;
  sig_mc_n_threshold_muon = BOGUS_INDEX;
  sig_mc_n_threshold_proton = BOGUS_INDEX;
  sig_mc_n_threshold_pion0 = BOGUS_INDEX;
  sig_mc_n_threshold_pionpm = BOGUS_INDEX;
  sig_mc_n_heaviermeson = BOGUS_INDEX;
  sig_mc_n_kaons = BOGUS_INDEX;
  sel_muoncandidate_tracklike_ = false;
  sel_pioncandidate_tracklike_ = false;
  sb_cc0pi_ = false; sb_multipi_ = false; sb_pi0_ = false; sb_cosmic_ = false;
  //CandidateMuonIndex = BOGUS_INDEX;
  //CandidatePionIndex = BOGUS_INDEX;
  pion_number = 0;
  sel_pion_contained = false;
  muon_in_gap = false;
  pion_in_gap = false;
  sel_MCVertexInFV = false;
  sig_recovertex_in_fv_ = false;
  shower_cut = false;
  shower_index = BOGUS_INDEX;
  nPrimaryTracks = 0;
  nPrimaryShowers = 0;
  opening_angle_cut = false;
  mu_pi_opening_angle = 0.0;
  true_mu_pi_opening_angle = 0.0;
  nonproton = 0;
  tmvaOutput_mip = 0.0;
  tmvaOutput = 0.0;
  tmvaOutput_pi = 0.0;
  candidate_muon_mom_mcs = 0.0;
  candidate_muon_mom_true = 0.0;
  candidate_muon_mom_range = 0.0;
  candidate_muon_mom_reco = 0.0;
  // These two were never reset. selection() only ever assigns them true, so
  // once one event set them they stayed true for every subsequent event --
  // sel_ntrack_gt_2_ in particular gates the muon/pion candidate-finding
  // block, so a later event with too few tracks would still run it. They were
  // also read uninitialised on the very first event.
  sel_ntrack_gt_2_ = false;
  sel_nshower_eq_0_ = false;
  candidate_pion_mom_reco = -999.0;
  candidate_pion_mom_true = -999.0;
  candidate_muon_costh_reco = -999.0;
  candidate_muon_costh_true = -999.0;
  candidate_pion_costh_reco = -999.0;
  candidate_pion_costh_true = -999.0;
  truemuonindex = -1;
  CandidateMuonTrackEndContainment = false;
  CandidateMuonIndex= -1;
  CandidatePionIndex = -1;
  TrueCandidateMuonP.SetXYZ(-1.0,-1.1,-1.0);
  TrueCandidatePionP.SetXYZ(-1.0,-1.1,-1.0);
 // sig_truevertex_fv = 0;
 // sig_ccnc = 0;
 // sig_numu = 0;
 // sig_one_muon = 0;

}

void CC1mu1piXp::define_category_map() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}


