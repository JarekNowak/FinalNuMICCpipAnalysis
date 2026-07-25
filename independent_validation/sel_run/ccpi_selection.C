// ccpi_selection.C
// Event-loop selection for CCπ⁺ analysis on MicroBooNE NuMI FHC data.
//
// Reads:  MC overlay, EXT beam-off, and dirt ROOT files.
// Writes: histograms.root      — 1D cut distributions for all variables/samples
//         unfolding_inputs.root — response-matrix seed and reco spectra
//                                 consumed by ccpi_xsec.C
//
// Factorisation note:
//   This macro contains NO unfolding mathematics.  All physics quantities
//   needed downstream (h_smear_sel, h_true_gen, h_reco_*) are written to
//   unfolding_inputs.root with their Sumw2 error arrays preserved.

#include "TFile.h"
#include "FV_new.h"
#include "TVector3.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <map>
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <TMVA/Reader.h>
#include <vector>

// ─── Constants ───────────────────────────────────────────────────────────────

static const int   NCUTS      = 11;
static const int   NVARS      = 60;
static const int   NSAMPLES   = 6;

static const double TOPO_CUT        = 0.67;
static const double MUON_TRK_SCORE  = 0.8;
static const double MUON_TRK_DIST   = 4.0;   // cm
static const double MUON_TRK_LEN    = 10.0;  // cm
static const double MUON_PID        = 0.2;
static const double MUON_BDT_MIP    = -0.1;
static const double PION_TRK_LEN    = 20.0;  // cm
static const double PION_TRK_DIST   = 4.0;   // cm
static const double PION_LLR        = 0.1;
static const double PION_BDT_MIP    = -0.1;
static const double PION_BDT_PI     = -0.1;
static const double OA_CUT          = 2.6;   // rad (upper μ–π opening-angle cut;
                                             // rejects the large-angle cosmic tail)
                                             // No lower OA cut: the small-angle region
                                             // is high-purity (~58%), so a lower cut
                                             // would only cost signal.
static const double M_PROTON        = 0.938272; // GeV
static const double PI_VAL          = 3.14159265358979;

// Shared binning and exposure constants — edit in CCPiConfig.h only.
#include "../CCPiConfig.h"

// ─── Background categories (truth-based) ──────────────────────────────────────
// DIAGNOSTIC ONLY — these labels are derived from MC truth and are used solely
// to break down the cut-flow so we can see WHICH backgrounds dominate at each
// stage.  They are never used as selection cuts.  Categories are mutually
// exclusive and assigned in priority order (see the dispatch in the event loop).
enum BkgCat {
    CAT_SIGNAL = 0,   // CC νμ, 1μ 1π± 0π0 0K in FV, pμ>0.15, pπ>0.175, θμπ<2.6 rad
    CAT_OOFV,         // νμ interaction with true vertex outside the fiducial volume
    CAT_NUE_OTHER,    // wrong flavour / sign (νe, anti-νμ, …): nu_pdg != 14
    CAT_NC,           // neutral-current interaction in FV
    CAT_CC_PI0,       // CC νμ in FV with at least one π0
    CAT_CC_0PI,       // CC νμ in FV with no charged pion (protons faking the pion)
    CAT_CC_MULTIPI,   // CC νμ in FV with ≥2 charged pions
    CAT_CC_KAON,      // CC νμ in FV with a strange meson (kaon)
    CAT_CC_OTHER,     // CC νμ in FV, 1π± but failing pμ>0.15/pπ>0.175 or other edge cases
    CAT_EXT,          // beam-off (cosmic) data
    CAT_DIRT,         // dirt MC
    CAT_DATA,         // beam-on data
    NCATS
};
static const char* CatName[NCATS] = {
    "Signal","OOFV","NuE/oth","NC",
    "CC1pi0","CC0pi","CCNpi","CCK","CCoth",
    "EXT","Dirt","Data"
};

// ─── Per-event result struct (avoids recomputing inside cut/sample loops) ────

struct EventResult {
    bool   cuts[NCUTS];
    bool   signal;
    int    muon_index, pion_index, shower_index;
    int    nPrimaryShowers, nPrimaryTracks, nonproton;
    double mu_pi_opening_angle;
    double mc_muon_momentum, mc_pion_momentum, mc_opening_angle;
    double scale;
    double ppfx, wtune;
    // Cached kinematics for histogram filling
    double muon_len, muon_pid, muon_vtx_dist, muon_bdt_mip;
    double muon_range_mom, muon_mcs_mom;
    double muon_theta, muon_costheta, muon_phi;
    double muon_sx, muon_sy, muon_sz;
    double pion_len, pion_pid, pion_vtx_dist;
    double pion_bdt_mip, pion_bdt_pi;
    double pion_range_mom, pion_mcs_mom;
    double pion_theta, pion_costheta, pion_phi;
    double pion_sx, pion_sy, pion_sz;
    double mu_pi_dist;
    int    muon_hits_u, muon_hits_v, muon_hits_y;
    int    pion_hits_u, pion_hits_v, pion_hits_y;
    int    shower_hits_u, shower_hits_v, shower_hits_y;
    double shower_vtx_dist;
    double shower_theta, shower_costheta, shower_phi;
    double topo_score;
    // MC truth
    double true_mu_p, true_mu_costh, true_mu_theta, true_mu_phi;
    double true_pi_p, true_pi_costh, true_pi_theta, true_pi_phi;
    double nu_energy, Q2, W;
    double nu_mu_mc_e, nu_e_mc_e;
    int    first_pdg;
};

// ─── Helper: safe weight ─────────────────────────────────────────────────────

inline double safeWeight(double w) {
    // Clamp to the framework's [MIN_WEIGHT, MAX_WEIGHT] = [0, 30] window
    // (XSecAnalyzer UniverseMaker::safe_weight, inclusive bounds); non-finite or
    // out-of-range weights are reset to unity.
    return (std::isfinite(w) && w >= 0.0 && w <= 30.0) ? w : 1.0;
}

// ─── Helper: meson identification (PDG-2019 numbering) ───────────────────────
// Mirrors XSecAnalyzer Functions.hh::is_meson_or_antimeson so the signal
// definition can veto heavy mesons (η, ρ, …) consistently with the framework.
inline bool is_meson_or_antimeson(int pdg_code) {
    int abs_pdg = std::abs(pdg_code);
    if (abs_pdg >= 9900000) return false;          // generator-specific
    if (((abs_pdg / 1000) % 10) != 0) return false; // n_q1 (thousands) must be 0
    if (((abs_pdg / 100)  % 10) == 0) return false; // n_q2 (hundreds) must be nonzero
    if (abs_pdg >= 901 && abs_pdg <= 930) return false; // SM PDF codes
    if (abs_pdg == 110 || abs_pdg == 990) return false; // reggeon / pomeron
    if (abs_pdg == 998 || abs_pdg == 999) return false; // GEANT tracking
    if (abs_pdg == 100) return false;              // generator pseudoparticle
    return true;
}

// ─── Main ────────────────────────────────────────────────────────────────────

// pion_bragg_cut : proton-rejection threshold on the pion candidate
//                  (trk_bragg_pion >= cut).  0 disables it.  Tunable so we can
//                  compare conservative (~0.08) vs aggressive (~0.20) working points.
// cosmic_bdt_cut : cosmic-rejection threshold at the NeutrinoSlice stage — events
//                  with bdt_cosmic <= cut are rejected.  bdt_cosmic is the dedicated
//                  cosmic-vs-neutrino BDT (→1 neutrino-like, →0 cosmic-like) and is
//                  well-defined for every event that has a reconstructed slice, in
//                  all samples (MC/EXT/dirt/data).  A NON-POSITIVE cut (<= 0)
//                  DISABLES the veto entirely — important because bdt_cosmic spans
//                  [-1,1], so a literal "bdt_cosmic > 0.0" test would silently
//                  reject every event with a negative score rather than being a
//                  no-op.
//                  REPLACES the former crt_pe_cut: the CRT branches (crthitpe /
//                  crtveto) are identically zero in these NuMI ntuples, so the old
//                  CRT veto rejected nothing.  Retune with selection/phase1_roc.C
//                  (the "bdtc" cosmic-rejection row).
// DEFAULT = 0.0 (DISABLED).  An A/B run (2026-06-09) showed bdt_cosmic>0.2 nets
// negative: it removes ~14% of EXT at the slice stage, but every one of those EXT
// events is already killed by the downstream topological/containment/opening-angle
// cuts, so final-cut EXT is IDENTICAL (40.9 either way) while signal drops 0.7%.
// The topological score is doing the cosmic rejection; the residual EXT is
// irreducible by this BDT (the survivors already pass >0.2).  Parameter kept for
// future studies but left off by default.
void ccpi_selection(double pion_bragg_cut = 0.08, double cosmic_bdt_cut = 0.0) {

    std::cout << "[ccpi_selection] pion_bragg_cut = " << pion_bragg_cut
              << "   cosmic_bdt_cut = " << cosmic_bdt_cut << "\n";

    // ── Histogram names ──────────────────────────────────────────────────────

    const char* CutsName[NCUTS] = {
        "EventsInTrueFV","NeutrinoSlice","InFiducialVol","Topological",
        "MuonCandidate","ContainedPion","MuonIn3Planes","PionIn3Planes",
        "ShowerCut","OpeningAngle","Nonprotons"
    };

    const char* Variable[NVARS] = {
        "TopologicalScore","MuonTrackLength","PionTrackLength","MuPiOpeningangle",
        "MuonPID","PionPID","MuonVtxDistance","PionVtxDistance","MuonPionDistance",
        "MuonUPlaneHits","MuonYPlaneHits","MuonVPlaneHits",
        "PionUPlaneHits","PionYPlaneHits","PionVPlaneHits",
        "ShowerUPlaneHits","ShowerYPlaneHits","ShowerVPlaneHits",
        "NPrimaryShowers","NPrimaryTracks","ShowerVtxDistance",
        "PionTMVAMip","PionTMVAPi","MuonTMVAMip",
        "PionTrkMuonMom","PionMCSMuonMom","MuonTrkMCSMom","MuonMCSMuonMom",
        "MCMuonMomentum","MCPionMomentum","MCOpeningAngle",
        "Nu_muMC_Eenergy","Nu_eMC_Energy","MuToPiDistance",
        "MuonTheta","MuonCosTheta","MuonPhi",
        "PionTheta","PionCosTheta","PionPhi",
        "ShowerTheta","ShowerCosTheta","ShowerPhi",
        "MuonStartX","MuonStartY","MuonStartZ",
        "PionStartX","PionStartY","PionStartZ",
        "TrueMuonMomentum","TrueMuonCosTheta","TrueMuonTheta","TrueMuonPhi",
        "TruePionMomentum","TruePionCosTheta","TruePionTheta","TruePionPhi",
        "NuEnergy","Q2","W"
    };

    const char* Sample[NSAMPLES] = {
        "AllMC","Signal","Background","EXT","OOFV","Data"
    };

    // ── Fixed per-variable axis ranges (index matches Variable[]) ────────────
    // Every sample of a given variable MUST share identical binning, otherwise
    // the histograms cannot be stacked or divided (THStack / Data-MC ratio in
    // plot_selection.C).  The previous "100, -1, -1" auto-range gave each sample
    // its own edges — fine for unfilled overlays, broken for stacks/ratios.
    const int    VarNb[NVARS] = {
        100,100,100,100,100,100,100,100,100,
        100,100,100,100,100,100,100,100,100,
         10, 12,100,
        100,100,100,
        100,100,100,100,
        100,100,180,
        100,100,100,
        100,100,100,
        100,100,100,
        100,100,100,
        100,100,100,
        100,100,100,
        100,100,100,100,
        100,100,100,100,
        100,100,100
    };
    const double VarLo[NVARS] = {
        0,0,0,0,-1,-1,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,
        -1,-1,-1,
        0,0,0,0,
        0,0,0,
        0,0,0,
        0,-1,-3.2,
        0,-1,-3.2,
        0,-1,-3.2,
        0,-120,0,
        0,-120,0,
        0,-1,0,-3.2,
        0,-1,0,-3.2,
        0,0,0.8
    };
    const double VarHi[NVARS] = {
        1,1000,300,3.2,1,1,5,5,20,
        2000,2000,2000,1000,1000,1000,500,500,500,
        10,12,50,
        1,1,1,
        1,1.5,2,3.2,
        3.2,1.5,180,
        10,10,20,
        3.2,1,3.2,
        3.2,1,3.2,
        3.2,1,3.2,
        256,120,1040,
        256,120,1040,
        3.2,1,3.2,3.2,
        1.5,1,3.2,3.2,
        10,2,2.5
    };

    // ── Histograms ───────────────────────────────────────────────────────────
    // Allocate as flat array [cut][sample][var] then fill by pointer for speed

    TH1F* Histos[NCUTS][NSAMPLES][NVARS];
    for (int c = 0; c < NCUTS; c++) {
        for (int v = 0; v < NVARS; v++) {
            TString base = TString(Variable[v]) + "_";
            for (int s = 0; s < NSAMPLES; s++) {
                TString name = base + CutsName[c] + "_" + Sample[s];
                Histos[c][s][v] = new TH1F(name, "", VarNb[v], VarLo[v], VarHi[v]);

                if (s != 4) Histos[c][s][v]->SetLineColor(s + 1);
                Histos[c][s][v]->GetXaxis()->SetTitle(Variable[v]);
                Histos[c][s][v]->SetTitle(CutsName[c]);
            }
        }
    }

    // ── Files and POT ────────────────────────────────────────────────────────

    const std::vector<std::string> Files = {
        "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/EXT/beamoff_run1Andrun3.root",
        "/data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root"
        // Data slot (i_f==6): NuWro fake data used IN PLACE OF beam-on data for a
        // generator-closure test.  Real beam-on data was:
        //   /data/uboone/beam_on/neutrinoselection_filt_run1_beamon_beamgood.root
        ,"/data/uboone/nuwro/numi_nuwro_overlay_pion_ntuples_run1_fhc.root" // NuWro fake data
    };

    // Only Run 1 data is available (neutrinoselection_filt_run1_beamon_beamgood.root).
    // Set Run 2/4/5 data POT to zero so their MC events get zero weight.
    const double dataPoT[4] = { 3.283, 0.0, 0.0, 0.0 };
    const double MCPoT[4]   = { 23.282, 24.9337, 11.04+17.295, 63.145 };
    double totalDataPoT = 0.0, totalMCPoT = 0.0;
    double Scale[7] = {1,1,1,1,1,1,1};
    for (int i = 0; i < 4; i++) {
        Scale[i]     = dataPoT[i] / MCPoT[i];
        totalDataPoT += dataPoT[i];
        totalMCPoT   += MCPoT[i];
    }
    Scale[4] = 7809962.0 / (610496.0 + 3211097.0);  // EXT trigger scaling
    Scale[5] = totalDataPoT / 16.739;                 // Dirt
    // NuWro fake data in the data slot (i_f==6): the NuMI FHC Run1 NuWro overlay
    // has 6.6504e20 generated POT (summed from its nuselection/SubRun "pot" branch).
    // Scale it to the data exposure so its reco spectrum is normalised like beam-on
    // data would be.  ppfx_cv and weightTune CV weights are applied automatically in
    // the event loop (the branches exist in the overlay), matching the xsec_analyzer
    // fake-data CV convention.  (REAL beam-on data lacks those branches, so the slot
    // would instead get unit weight via the ppfx/weightTune reset before GetEntry.)
    Scale[6] = totalDataPoT / 6.6504;                 // NuWro fake data -> data POT

    double Selected[NCUTS][6] = {{0}};

    // Diagnostic-only accumulators (Phase 0 purity study).
    // SelectedCat[cut][category] : POT-weighted events per truth category per cut.
    // pionTruthW[|pdg|]          : truth PDG of the reco pion CANDIDATE for all MC
    //                              events at the final cut (measures, e.g., how
    //                              many "pions" are really protons).
    double SelectedCat[NCUTS][NCATS] = {{0}};
    std::map<int, double> pionTruthW;

    // ── Event-display list ───────────────────────────────────────────────────
    // Lightweight per-event record for the 2D event-display tool
    // (selection/event_display.C): every event passing the FINAL cut, except
    // beam-on data, is recorded with its parent file + entry number (for direct
    // re-access), the muon/pion candidate indices, and its truth category.  The
    // display reopens the ntuple and draws the PFPs truth-labelled by backtracked
    // PDG, so it always shows exactly the events entering the measurement.
    std::vector<int>         evl_run_v, evl_sub_v, evl_evt_v, evl_fidx_v,
                             evl_cat_v, evl_mi_v, evl_pi_v;
    std::vector<Long64_t>    evl_entry_v;
    std::vector<std::string> evl_fname_v;

    // ── Unfolding-input histograms, one set per observable ───────────────────
    // For each observable o (p_μ, p_π, cosθ_μ, cosθ_π, θ_μπ — see CCPiConfig.h):
    //   h_smear[o]   : 2D response seed (x = true, y = reco), signal MC at Cut[10]
    //   h_truegen[o] : 1D true value for ALL signal events (efficiency denominator)
    //   h_reco[o][s] : 1D reco spectra per sample s at Cut[10]
    // Sample index: 0=signal 1=bkg 2=EXT 3=dirt 4=data.
    const char* RECO_SUF[5] = { "sig", "bkg", "ext", "dirt", "data" };
    TH2D* h_smear  [N_OBS];
    TH1D* h_truegen[N_OBS];
    TH1D* h_reco   [N_OBS][5];
    for (int o = 0; o < N_OBS; o++) {
        const int      nb = OBS_NBINS[o];
        const double*  ed = OBS_EDGES(o);
        const char*    key = OBS_KEY[o];
        const char*    ttl = OBS_TITLE[o];
        h_smear[o] = new TH2D(Form("h_smear_sel_%s", key),
            Form("Response seed;True %s;Reco %s", ttl, ttl), nb, ed, nb, ed);
        h_truegen[o] = new TH1D(Form("h_true_gen_%s", key),
            Form("True signal (all generated);True %s;Events", ttl), nb, ed);
        for (int s = 0; s < 5; s++)
            h_reco[o][s] = new TH1D(Form("h_reco_%s_%s", RECO_SUF[s], key),
                Form("Reco %s;Reco %s;Events", RECO_SUF[s], ttl), nb, ed);
        h_smear[o]->Sumw2();
        h_truegen[o]->Sumw2();
        for (int s = 0; s < 5; s++) h_reco[o][s]->Sumw2();
    }

    // ── Phase-1 purity-study diagnostic histograms ───────────────────────────
    // Written to phase1_diagnostics.root; analysed offline to choose cut
    // thresholds (ROC: signal/pion-keep efficiency vs background rejection).
    //   *_pi  : pion-candidate variable when it REALLY is a π± (truth) — keep
    //   *_pr  : ... when it really is a proton — reject
    //   *_nu  : event-level cosmic variable for true SIGNAL events — keep
    //   *_ext : ... for EXT (beam-off cosmic) events — reject
    auto book = [](const char* nm, int nb, double lo, double hi) {
        TH1D* h = new TH1D(nm, "", nb, lo, hi);
        h->Sumw2(); h->SetDirectory(nullptr); return h;
    };
    // Pion-candidate PID at the final cut
    TH1D* h_chipr_pi = book("ph1_chipr_pi", 80, 0, 400);
    TH1D* h_chipr_pr = book("ph1_chipr_pr", 80, 0, 400);
    TH1D* h_chipi_pi = book("ph1_chipi_pi", 80, 0, 400);
    TH1D* h_chipi_pr = book("ph1_chipi_pr", 80, 0, 400);
    TH1D* h_chimu_pi = book("ph1_chimu_pi", 80, 0, 400);
    TH1D* h_chimu_pr = book("ph1_chimu_pr", 80, 0, 400);
    TH1D* h_pida_pi  = book("ph1_pida_pi",  60, 0,  30);
    TH1D* h_pida_pr  = book("ph1_pida_pr",  60, 0,  30);
    TH1D* h_bpion_pi = book("ph1_bpion_pi", 50, 0,   1);
    TH1D* h_bpion_pr = book("ph1_bpion_pr", 50, 0,   1);
    // Cosmic-rejection variables at the MuonCandidate stage
    TH1D* h_flash_nu  = book("ph1_flash_nu",  60,  0,  60);
    TH1D* h_flash_ext = book("ph1_flash_ext", 60,  0,  60);
    TH1D* h_crtpe_nu  = book("ph1_crtpe_nu",  60,  0, 600);
    TH1D* h_crtpe_ext = book("ph1_crtpe_ext", 60,  0, 600);
    TH1D* h_cosip_nu  = book("ph1_cosip_nu",  60,  0, 200);
    TH1D* h_cosip_ext = book("ph1_cosip_ext", 60,  0, 200);
    TH1D* h_bdtc_nu   = book("ph1_bdtc_nu",   60, -1,   1);
    TH1D* h_bdtc_ext  = book("ph1_bdtc_ext",  60, -1,   1);

    // Variables are set per-track before EvaluateMVA; readers are stateless
    // between evaluations so they're safe to reuse.

    float f_bragg_p, f_bragg_mu, f_bragg_mip, f_llr, f_score;
    float f_end_x, f_end_y, f_end_z;

    TMVA::Reader* tmvaMIP = new TMVA::Reader();
    tmvaMIP->AddVariable("trk_bragg_p_v",       &f_bragg_p);
    tmvaMIP->AddVariable("trk_bragg_mu_v",      &f_bragg_mu);
    tmvaMIP->AddVariable("trk_bragg_mip_v",     &f_bragg_mip);
    tmvaMIP->AddVariable("trk_llr_pid_score_v", &f_llr);
    tmvaMIP->AddVariable("trk_score_v",         &f_score);
    tmvaMIP->AddVariable("trk_sce_end_x_v",     &f_end_x);
    tmvaMIP->AddVariable("trk_sce_end_y_v",     &f_end_y);
    tmvaMIP->AddVariable("trk_sce_end_z_v",     &f_end_z);
    tmvaMIP->BookMVA("BDT", "booster_decision_tree/dataset_MIP_BDT_no_len/weights/TMVAClassification_BDT.weights.xml");

    float fp_bragg_p, fp_bragg_mu, fp_bragg_mip, fp_llr, fp_score;
    float fp_end_x, fp_end_y, fp_end_z;

    TMVA::Reader* tmvaPI = new TMVA::Reader();
    tmvaPI->AddVariable("trk_bragg_p_v",       &fp_bragg_p);
    tmvaPI->AddVariable("trk_bragg_mu_v",      &fp_bragg_mu);
    tmvaPI->AddVariable("trk_bragg_mip_v",     &fp_bragg_mip);
    tmvaPI->AddVariable("trk_llr_pid_score_v", &fp_llr);
    tmvaPI->AddVariable("trk_score_v",         &fp_score);
    tmvaPI->AddVariable("trk_sce_end_x_v",     &fp_end_x);
    tmvaPI->AddVariable("trk_sce_end_y_v",     &fp_end_y);
    tmvaPI->AddVariable("trk_sce_end_z_v",     &fp_end_z);
    tmvaPI->BookMVA("BDT", "booster_decision_tree/dataset_pion_BDT_no_len/weights/TMVAClassification_BDT.weights.xml");

    // ─── Branch variables (declared once, reused across files) ───────────────

    Int_t   run, sub, evt, nu_pdg, ccnc, interaction;
    Int_t   swtrig;   // software trigger flag (==1 passes); data/EXT are pre-filtered
    Float_t nu_e;   // true neutrino energy (GeV)
    Float_t true_nu_px, true_nu_py, true_nu_pz;
    Float_t reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
    Int_t   n_tracks, n_showers;
    Float_t topological_score;
    Float_t ppfx_cv, weightTune, weightSpline;

    std::vector<int>*   backtracked_pdg    = nullptr;
    std::vector<float>* backtracked_px     = nullptr;
    std::vector<float>* backtracked_py     = nullptr;
    std::vector<float>* backtracked_pz     = nullptr;
    std::vector<float>* backtracked_e      = nullptr;
    std::vector<float>* backtracked_start_x= nullptr;
    std::vector<float>* backtracked_start_y= nullptr;
    std::vector<float>* backtracked_start_z= nullptr;
    std::vector<int>*   mc_pdg             = nullptr;
    std::vector<float>* mc_E               = nullptr;
    std::vector<float>* mc_vx              = nullptr;
    std::vector<float>* mc_vy              = nullptr;
    std::vector<float>* mc_vz              = nullptr;
    std::vector<float>* mc_px              = nullptr;
    std::vector<float>* mc_py              = nullptr;
    std::vector<float>* mc_pz              = nullptr;
    std::vector<float>* trk_score_v        = nullptr;
    std::vector<float>* trk_dir_x_v        = nullptr;
    std::vector<float>* trk_dir_y_v        = nullptr;
    std::vector<float>* trk_dir_z_v        = nullptr;
    std::vector<float>* trk_sce_start_x_v  = nullptr;
    std::vector<float>* trk_sce_start_y_v  = nullptr;
    std::vector<float>* trk_sce_start_z_v  = nullptr;
    std::vector<float>* trk_sce_end_x_v    = nullptr;
    std::vector<float>* trk_sce_end_y_v    = nullptr;
    std::vector<float>* trk_sce_end_z_v    = nullptr;
    std::vector<float>* trk_len_v          = nullptr;
    std::vector<float>* trk_calo_energy_u_v= nullptr;
    std::vector<float>* trk_calo_energy_v_v= nullptr;
    std::vector<float>* trk_calo_energy_y_v= nullptr;
    std::vector<float>* trk_llr_pid_score_v= nullptr;
    std::vector<float>* trk_bragg_p_v      = nullptr;
    std::vector<float>* trk_bragg_mu_v     = nullptr;
    std::vector<float>* trk_bragg_mip_v    = nullptr;
    std::vector<float>* trk_mcs_muon_mom_v = nullptr;
    std::vector<float>* trk_range_muon_mom_v=nullptr;
    // Phase-1 purity-study candidates (proton rejection on the pion candidate)
    std::vector<float>* trk_pid_chipr_v    = nullptr;
    std::vector<float>* trk_pid_chipi_v    = nullptr;
    std::vector<float>* trk_pid_chimu_v    = nullptr;
    std::vector<float>* trk_pida_v         = nullptr;
    std::vector<float>* trk_bragg_pion_v   = nullptr;
    std::vector<float>* trk_energy_proton_v= nullptr;
    std::vector<float>* trk_energy_muon_v  = nullptr;
    // Phase-1 purity-study candidates (cosmic rejection — event-level)
    Int_t   crtveto             = 0;
    Float_t crthitpe            = 0.f;
    Float_t nu_flashmatch_score = 0.f;
    Float_t bdt_cosmic          = 0.f;
    Float_t CosmicIP            = 0.f;
    std::vector<int>*   pfnplanehits_U     = nullptr;
    std::vector<int>*   pfnplanehits_V     = nullptr;
    std::vector<int>*   pfnplanehits_Y     = nullptr;
    std::vector<unsigned int>* pfp_generation_v = nullptr;
    std::vector<float>* shr_start_x_v      = nullptr;
    std::vector<float>* shr_start_y_v      = nullptr;
    std::vector<float>* shr_start_z_v      = nullptr;
    std::vector<float>* shr_dist_v         = nullptr;
    std::vector<float>* shr_moliere_avg_v  = nullptr;
    std::vector<float>* shr_moliere_rms_v  = nullptr;
    std::vector<float>* shr_tkfit_dedx_u_v = nullptr;
    std::vector<float>* shr_tkfit_dedx_v_v = nullptr;
    std::vector<float>* shr_tkfit_dedx_y_v = nullptr;

    // ─── File loop ────────────────────────────────────────────────────────────

    for (size_t i_f = 0; i_f < Files.size(); i_f++) {

        // eff/purity cross-check: Run2/4/5 MC (i_f 1,2,3) carry zero weight
        // (their data POT is 0) and the NuWro data slot (i_f 6) does not enter
        // signal efficiency or purity, so skip them for speed. This leaves
        // Run1 MC + EXT + dirt and gives IDENTICAL eff/purity to the full run.
        if (i_f == 1 || i_f == 2 || i_f == 3 || i_f == 6) continue;

        TFile* f = TFile::Open(Files[i_f].c_str());
        if (!f || f->IsZombie()) { std::cerr << "Cannot open " << Files[i_f] << "\n"; if(f) delete f; continue; }

        TTree* t = nullptr;
        f->GetObject("nuselection/NeutrinoSelectionFilter", t);
        if (!t) { std::cerr << "No tree in " << Files[i_f] << "\n"; delete f; continue; }
        std::cout << "Processing: " << Files[i_f] << " (" << t->GetEntries() << " events)\n";

        t->SetMakeClass(1);
        t->SetBranchStatus("*", 0);

        // Enable only what we read
        const std::vector<const char*> branches = {
            "run","sub","evt","nu_pdg","ccnc","interaction","nu_e","swtrig",
            "true_nu_px","true_nu_py","true_nu_pz",
            "reco_nu_vtx_sce_x","reco_nu_vtx_sce_y","reco_nu_vtx_sce_z",
            "backtracked_pdg","backtracked_start_x","backtracked_start_y","backtracked_start_z",
            "backtracked_px","backtracked_py","backtracked_pz","backtracked_e",
            "n_tracks","n_showers","trk_score_v","mc_pdg","mc_E",
            "mc_vx","mc_vy","mc_vz","mc_px","mc_py","mc_pz",
            "trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v",
            "trk_sce_end_x_v","trk_sce_end_y_v","trk_sce_end_z_v",
            "trk_dir_x_v","trk_dir_y_v","trk_dir_z_v","trk_len_v",
            "trk_calo_energy_u_v","trk_calo_energy_v_v","trk_calo_energy_y_v",
            "trk_llr_pid_score_v","topological_score",
            "pfnplanehits_U","pfnplanehits_V","pfnplanehits_Y",
            "pfp_generation_v","trk_bragg_p_v","trk_bragg_mu_v","trk_bragg_mip_v",
            "shr_start_x_v","shr_start_y_v","shr_start_z_v","shr_dist_v",
            "backtracked_purity","shr_moliere_avg_v","shr_moliere_rms_v",
            "shr_tkfit_dedx_u_v","shr_tkfit_dedx_v_v","shr_tkfit_dedx_y_v",
            "trk_mcs_muon_mom_v","trk_range_muon_mom_v","ppfx_cv","weightTune","weightSpline",
            "trk_pid_chipr_v","trk_pid_chipi_v","trk_pid_chimu_v","trk_pida_v",
            "trk_bragg_pion_v","trk_energy_proton_v","trk_energy_muon_v",
            "crtveto","crthitpe","nu_flashmatch_score","bdt_cosmic","CosmicIP"
        };
        for (auto& b : branches) t->SetBranchStatus(b, 1);

        t->SetBranchAddress("run",                  &run);
        t->SetBranchAddress("sub",                  &sub);
        t->SetBranchAddress("evt",                  &evt);
        t->SetBranchAddress("nu_pdg",               &nu_pdg);
        t->SetBranchAddress("ccnc",                 &ccnc);
        t->SetBranchAddress("interaction",          &interaction);
        t->SetBranchAddress("nu_e",                 &nu_e);
        t->SetBranchAddress("swtrig",               &swtrig);
        t->SetBranchAddress("true_nu_px",           &true_nu_px);
        t->SetBranchAddress("true_nu_py",           &true_nu_py);
        t->SetBranchAddress("true_nu_pz",           &true_nu_pz);
        t->SetBranchAddress("reco_nu_vtx_sce_x",    &reco_nu_vtx_sce_x);
        t->SetBranchAddress("reco_nu_vtx_sce_y",    &reco_nu_vtx_sce_y);
        t->SetBranchAddress("reco_nu_vtx_sce_z",    &reco_nu_vtx_sce_z);
        t->SetBranchAddress("backtracked_pdg",      &backtracked_pdg);
        t->SetBranchAddress("backtracked_start_x",  &backtracked_start_x);
        t->SetBranchAddress("backtracked_start_y",  &backtracked_start_y);
        t->SetBranchAddress("backtracked_start_z",  &backtracked_start_z);
        t->SetBranchAddress("backtracked_px",       &backtracked_px);
        t->SetBranchAddress("backtracked_py",       &backtracked_py);
        t->SetBranchAddress("backtracked_pz",       &backtracked_pz);
        t->SetBranchAddress("backtracked_e",        &backtracked_e);
        t->SetBranchAddress("n_tracks",             &n_tracks);
        t->SetBranchAddress("n_showers",            &n_showers);
        t->SetBranchAddress("trk_score_v",          &trk_score_v);
        t->SetBranchAddress("mc_pdg",               &mc_pdg);
        t->SetBranchAddress("mc_E",                 &mc_E);
        t->SetBranchAddress("mc_vx",                &mc_vx);
        t->SetBranchAddress("mc_vy",                &mc_vy);
        t->SetBranchAddress("mc_vz",                &mc_vz);
        t->SetBranchAddress("mc_px",                &mc_px);
        t->SetBranchAddress("mc_py",                &mc_py);
        t->SetBranchAddress("mc_pz",                &mc_pz);
        t->SetBranchAddress("trk_sce_start_x_v",    &trk_sce_start_x_v);
        t->SetBranchAddress("trk_sce_start_y_v",    &trk_sce_start_y_v);
        t->SetBranchAddress("trk_sce_start_z_v",    &trk_sce_start_z_v);
        t->SetBranchAddress("trk_sce_end_x_v",      &trk_sce_end_x_v);
        t->SetBranchAddress("trk_sce_end_y_v",      &trk_sce_end_y_v);
        t->SetBranchAddress("trk_sce_end_z_v",      &trk_sce_end_z_v);
        t->SetBranchAddress("trk_dir_x_v",          &trk_dir_x_v);
        t->SetBranchAddress("trk_dir_y_v",          &trk_dir_y_v);
        t->SetBranchAddress("trk_dir_z_v",          &trk_dir_z_v);
        t->SetBranchAddress("trk_len_v",            &trk_len_v);
        t->SetBranchAddress("trk_calo_energy_u_v",  &trk_calo_energy_u_v);
        t->SetBranchAddress("trk_calo_energy_v_v",  &trk_calo_energy_v_v);
        t->SetBranchAddress("trk_calo_energy_y_v",  &trk_calo_energy_y_v);
        t->SetBranchAddress("trk_llr_pid_score_v",  &trk_llr_pid_score_v);
        t->SetBranchAddress("topological_score",    &topological_score);
        t->SetBranchAddress("pfnplanehits_U",       &pfnplanehits_U);
        t->SetBranchAddress("pfnplanehits_V",       &pfnplanehits_V);
        t->SetBranchAddress("pfnplanehits_Y",       &pfnplanehits_Y);
        t->SetBranchAddress("pfp_generation_v",     &pfp_generation_v);
        t->SetBranchAddress("trk_bragg_p_v",        &trk_bragg_p_v);
        t->SetBranchAddress("trk_bragg_mu_v",       &trk_bragg_mu_v);
        t->SetBranchAddress("trk_bragg_mip_v",      &trk_bragg_mip_v);
        t->SetBranchAddress("shr_start_x_v",        &shr_start_x_v);
        t->SetBranchAddress("shr_start_y_v",        &shr_start_y_v);
        t->SetBranchAddress("shr_start_z_v",        &shr_start_z_v);
        t->SetBranchAddress("shr_dist_v",           &shr_dist_v);
        t->SetBranchAddress("shr_moliere_avg_v",    &shr_moliere_avg_v);
        t->SetBranchAddress("shr_moliere_rms_v",    &shr_moliere_rms_v);
        t->SetBranchAddress("shr_tkfit_dedx_u_v",   &shr_tkfit_dedx_u_v);
        t->SetBranchAddress("shr_tkfit_dedx_v_v",   &shr_tkfit_dedx_v_v);
        t->SetBranchAddress("shr_tkfit_dedx_y_v",   &shr_tkfit_dedx_y_v);
        t->SetBranchAddress("trk_mcs_muon_mom_v",   &trk_mcs_muon_mom_v);
        t->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v);
        t->SetBranchAddress("ppfx_cv",              &ppfx_cv);
        t->SetBranchAddress("weightTune",           &weightTune);
        t->SetBranchAddress("weightSpline",         &weightSpline);
        // Phase-1 purity-study candidate branches
        t->SetBranchAddress("trk_pid_chipr_v",      &trk_pid_chipr_v);
        t->SetBranchAddress("trk_pid_chipi_v",      &trk_pid_chipi_v);
        t->SetBranchAddress("trk_pid_chimu_v",      &trk_pid_chimu_v);
        t->SetBranchAddress("trk_pida_v",           &trk_pida_v);
        t->SetBranchAddress("trk_bragg_pion_v",     &trk_bragg_pion_v);
        t->SetBranchAddress("trk_energy_proton_v",  &trk_energy_proton_v);
        t->SetBranchAddress("trk_energy_muon_v",    &trk_energy_muon_v);
        t->SetBranchAddress("crtveto",              &crtveto);
        t->SetBranchAddress("crthitpe",             &crthitpe);
        t->SetBranchAddress("nu_flashmatch_score",  &nu_flashmatch_score);
        t->SetBranchAddress("bdt_cosmic",           &bdt_cosmic);
        t->SetBranchAddress("CosmicIP",             &CosmicIP);

        const Long64_t nentries = t->GetEntries();

        // eff/purity cross-check: process the same FRACTION of every file so the
        // POT-weighted signal:bkg:EXT:dirt normalisation (and hence the eff/purity
        // ratios) is preserved, while cutting the I/O over these multi-GB ntuples.
        // Set SEL_FRAC=1 in the environment for the full, exact run.
        double sel_frac = 0.30;
        if (const char* ef = std::getenv("SEL_FRAC")) sel_frac = std::atof(ef);
        const Long64_t maxev = (sel_frac >= 1.0) ? nentries
                             : (Long64_t)(nentries * sel_frac);

        // ─── Event loop ───────────────────────────────────────────────────────

        for (Long64_t ientry = 0; ientry < nentries; ientry++) {
            if (ientry >= maxev) break;

            // ppfx_cv, weightTune and weightSpline are MC-only branches.  Reset to
            // 1.0 before GetEntry so that data/EXT files (which lack these branches)
            // get unit weight instead of the stale value from the last MC event.
            ppfx_cv      = 1.0f;
            weightTune   = 1.0f;
            weightSpline = 1.0f;
            // swtrig resets fail-open (==1 passes) so a file lacking the branch
            // is not wiped out by the software-trigger cut.
            swtrig     = 1;
            // Phase-1 cosmic scalars: reset to neutral sentinels so a file that
            // lacks a branch doesn't inherit the previous event's value.
            crtveto = 0; crthitpe = 0.f; CosmicIP = 9999.f;
            // bdt_cosmic resets fail-open (large value passes the > cut) so a file
            // lacking the branch is not wiped out by the cosmic veto.
            nu_flashmatch_score = 0.f; bdt_cosmic = 9999.f;

            t->GetEntry(ientry);
            if (ientry % 100000 == 0) std::cout << "  " << ientry << " / " << nentries << "\n";

            // ── Weights ───────────────────────────────────────────────────────
            // Central-value weight: clamp the PRODUCT of the CV reweights once,
            // matching the framework (UniverseMaker applies safe_weight to the
            // combined spline*tune*ppfx value), rather than clamping each factor
            // separately. Scale is POT normalisation and stays outside the clamp.
            // Includes weightSpline (GENIE spline) to match the framework CV
            // weight_TunedCentralValue_UBGenie = spline*tuned_cv*ppfx.
            const double w_cv    = safeWeight((double)ppfx_cv * (double)weightTune
                                              * (double)weightSpline);
            const double evtW    = Scale[i_f] * w_cv;

            // Skip MC files whose run period is not in the data sample.
            // (dataPoT[1..3] == 0 → Scale[1..3] == 0 for Run 2/4/5 MC.)
            if (i_f < 4 && evtW < 1e-15) continue;

            // ── MC truth pass ─────────────────────────────────────────────────
            // Single loop over mc_pdg: count particle types, store momenta.

            bool   signal              = false;
            bool   in_fv_true          = false;
            int    nmuons=0, npions=0, nprotons=0, npionszero=0, nkaons=0, nheavymesons=0;
            double mc_muon_momentum    = -999.9;
            double mc_pion_momentum    = -999.9;
            double mc_opening_angle    = -999.9;
            TVector3 mc_muon_mom3, mc_pion_mom3;

            if (!mc_pdg->empty()) {
                for (size_t i_mc = 0; i_mc < mc_pdg->size(); i_mc++) {
                    int pdg = mc_pdg->at(i_mc);
                    // Count by |pdg| so the signal definition is charge-blind:
                    // it covers νμ (μ⁻,π⁺) and ν̄μ (μ⁺,π⁻) on the same footing.
                    const int apdg = std::abs(pdg);
                    if      (apdg == 2212)                             nprotons++;
                    else if (apdg == 13) {
                        nmuons++;
                        mc_muon_mom3.SetXYZ(mc_px->at(i_mc), mc_py->at(i_mc), mc_pz->at(i_mc));
                        mc_muon_momentum = mc_muon_mom3.Mag();
                    }
                    else if (apdg == 211) {
                        npions++;
                        mc_pion_mom3.SetXYZ(mc_px->at(i_mc), mc_py->at(i_mc), mc_pz->at(i_mc));
                        mc_pion_momentum = mc_pion_mom3.Mag();
                    }
                    else if (apdg == 111)                              npionszero++;
                    else if (apdg==321 || pdg==310 || pdg==130 || pdg==311) nkaons++;
                    // Any remaining meson (η, ρ, …) is a "heavy meson" — vetoed
                    // by the signal definition to match the framework.
                    else if (is_meson_or_antimeson(pdg))               nheavymesons++;
                }

                double mc_oa_rad = -1.0;   // true μ–π opening angle (rad)
                if (nmuons == 1 && npions == 1) {
                    mc_oa_rad = mc_pion_mom3.Angle(mc_muon_mom3);
                    mc_opening_angle = mc_oa_rad * 180.0 / PI_VAL;
                }

                const TVector3 true_vtx(mc_vx->at(0), mc_vy->at(0), mc_vz->at(0));
                in_fv_true = inFV(true_vtx);

                if (in_fv_true
                    && nmuons == 1 && npionszero == 0 && npions == 1 && nkaons == 0
                    && nheavymesons == 0
                    && std::abs(nu_pdg) == 14 && ccnc==0
                    && mc_muon_momentum > 0.15 && mc_pion_momentum > 0.175
                    && mc_oa_rad < OA_CUT)
                    signal = true;
            }

            // ── Truth background category (diagnostic only — not a cut) ───────
            // Priority order: signal first, then increasingly specific neutrino
            // backgrounds.  Cosmic/dirt/data are fixed by the input file index.
            int event_cat = -1;
            if (i_f < 4) {
                if      (signal)             event_cat = CAT_SIGNAL;
                else if (!in_fv_true)        event_cat = CAT_OOFV;
                else if (std::abs(nu_pdg) != 14) event_cat = CAT_NUE_OTHER;
                else if (ccnc == 1)          event_cat = CAT_NC;
                else if (npionszero >= 1)    event_cat = CAT_CC_PI0;
                else if (npions == 0)        event_cat = CAT_CC_0PI;
                else if (npions >= 2)        event_cat = CAT_CC_MULTIPI;
                else if (nkaons >= 1)        event_cat = CAT_CC_KAON;
                else                         event_cat = CAT_CC_OTHER;
            }
            else if (i_f == 4)               event_cat = CAT_EXT;
            else if (i_f == 5)               event_cat = CAT_DIRT;
            else if (i_f == 6)               event_cat = CAT_DATA;

            // ── True observable values (signal only) ──────────────────────────
            // Computed once here and reused below as the true axis of the
            // response matrix.  cosθ is w.r.t. detector z; θ_μπ in radians.
            double true_obs[N_OBS] = {0};
            if (i_f < 4 && signal) {
                true_obs[OBS_PMU]  = mc_muon_momentum;
                true_obs[OBS_PPI]  = mc_pion_momentum;
                true_obs[OBS_CTMU] = mc_muon_mom3.CosTheta();
                true_obs[OBS_CTPI] = mc_pion_mom3.CosTheta();
                true_obs[OBS_OA]   = mc_pion_mom3.Angle(mc_muon_mom3);  // rad
            }

            // ── True-space denominator for response matrix ────────────────────
            // Fill for every MC signal event regardless of reco selection.
            // Uses the same POT×ppfx×wtune weight as the response numerator so
            // that their ratio gives a pure selection probability.
            if (i_f < 4 && signal)
                for (int o = 0; o < N_OBS; o++)
                    h_truegen[o]->Fill(true_obs[o], evtW);

            // ── Reco vertex ───────────────────────────────────────────────────
            const double rvx = reco_nu_vtx_sce_x;
            const double rvy = reco_nu_vtx_sce_y;
            const double rvz = reco_nu_vtx_sce_z;
            const TVector3 recoVtx(rvx, rvy, rvz);
            const bool in_fv_reco = inFV(recoVtx);

            // ── Topological cut ───────────────────────────────────────────────
            const bool topo_ok = (topological_score > TOPO_CUT);

            // ── Single pass over tracks: compute BDTs once, cache per track ──
            // This is the main bottleneck in the original — each track was
            // evaluated up to 4 times. Here we evaluate once and cache.

            const int ntracks = (int)trk_len_v->size();
            bool flash_cut    = (ntracks > 0);

            // Per-track cache
            std::vector<float> bdt_mip_v(ntracks, 0.f);  // MIP BDT score
            std::vector<float> bdt_pi_v (ntracks, 0.f);  // pion BDT score
            std::vector<bool>  valid_v  (ntracks, false); // passes basic quality

            int muon_index   = -1;
            int pion_index   = -1;
            int shower_index = -1;
            int nPrimaryShowers = 0;
            int nPrimaryTracks  = 0;
            int nonproton       = 0;

            // Best muon/pion candidates tracked as we scan
            float best_muon_len = -1.f;
            float best_pion_llr = -1e9f;

            for (int i = 0; i < ntracks; i++) {
                if (pfp_generation_v->at(i) != 2) continue;

                const float score  = trk_score_v->at(i);
                const float len    = trk_len_v->at(i);
                const float llr    = trk_llr_pid_score_v->at(i);
                const float bp     = trk_bragg_p_v->at(i);
                const float bmu    = trk_bragg_mu_v->at(i);
                const float bmip   = trk_bragg_mip_v->at(i);
                const float ex     = trk_sce_end_x_v->at(i);
                const float ey     = trk_sce_end_y_v->at(i);
                const float ez     = trk_sce_end_z_v->at(i);

                if (score >= 0.5f) {
                    nPrimaryTracks++;
                    if (llr > 0.1f) nonproton++;

                    // Quality gate for BDT evaluation
                    const bool qualOK = (llr > -1.f && llr < 2.f)
                                     && (bp  > 0.f  && bp  < 500.f)
                                     && (bmu > 0.f  && bmu < 500.f)
                                     && (bmip > 0.f && bmip < 500.f)
                                     && (len > 0.f  && len < 1e6f);
                    valid_v[i] = qualOK;

                    if (qualOK) {
                        // Evaluate MIP BDT once
                        f_bragg_p = bp; f_bragg_mu = bmu; f_bragg_mip = bmip;
                        f_llr = llr;    f_score = score;
                        f_end_x = ex;   f_end_y = ey;   f_end_z = ez;
                        bdt_mip_v[i] = tmvaMIP->EvaluateMVA("BDT");

                        // Evaluate pion BDT once (only for contained tracks)
                        const TVector3 ep_c(ex, ey, ez);
                        if (isContained(ep_c)) {
                            fp_bragg_p = bp; fp_bragg_mu = bmu; fp_bragg_mip = bmip;
                            fp_llr = llr;    fp_score = score;
                            fp_end_x = ex;   fp_end_y = ey;   fp_end_z = ez;
                            bdt_pi_v[i] = tmvaPI->EvaluateMVA("BDT");
                        }

                        // ── Vertex distance (computed once) ───────────────────
                        const double dsx = trk_sce_start_x_v->at(i) - rvx;
                        const double dsy = trk_sce_start_y_v->at(i) - rvy;
                        const double dsz = trk_sce_start_z_v->at(i) - rvz;
                        const double vtx_dist = std::sqrt(dsx*dsx + dsy*dsy + dsz*dsz);

                        // ── Muon candidate ────────────────────────────────────
                        if (score > MUON_TRK_SCORE
                            && vtx_dist <= MUON_TRK_DIST
                            && len > MUON_TRK_LEN
                            && llr > MUON_PID
                            && bdt_mip_v[i] >= MUON_BDT_MIP) {
                            if (muon_index == -1 || len > best_muon_len) {
                                muon_index    = i;
                                best_muon_len = len;
                            }
                        }
                    }
                } else {
                    nPrimaryShowers++;
                    shower_index = i;
                }
            }

            // ── Pion pass (skip muon index, require contained endpoint) ───────
            // Separated to guarantee muon_index is finalized before pion scan.
            // Selects the highest-LLR pion and counts all qualifying candidates
            // in the same pass (pion_count is used for the one-pion requirement).

            int pion_count = 0;
            for (int i = 0; i < ntracks; i++) {
                if (pfp_generation_v->at(i) != 2) continue;
                if (trk_score_v->at(i) < 0.5f)    continue;
                if (i == muon_index)                continue;
                if (!valid_v[i])                    continue;

                const float ex  = trk_sce_end_x_v->at(i);
                const float ey  = trk_sce_end_y_v->at(i);
                const float ez  = trk_sce_end_z_v->at(i);
                const TVector3 ep(ex, ey, ez);
                if (!isContained(ep)) continue;

                const float llr     = trk_llr_pid_score_v->at(i);
                const float len     = trk_len_v->at(i);
                const double dsx    = trk_sce_start_x_v->at(i) - rvx;
                const double dsy    = trk_sce_start_y_v->at(i) - rvy;
                const double dsz    = trk_sce_start_z_v->at(i) - rvz;
                const double vtxd   = std::sqrt(dsx*dsx + dsy*dsy + dsz*dsz);

                // Proton rejection: protons have a low pion-Bragg likelihood.
                // Default 1.0 (passes) if the branch is unavailable in a file.
                const float bragg_pi = (trk_bragg_pion_v && i < (int)trk_bragg_pion_v->size())
                                       ? trk_bragg_pion_v->at(i) : 1.0f;

                if (llr       > PION_LLR
                    && vtxd   < PION_TRK_DIST
                    && bdt_mip_v[i] > PION_BDT_MIP
                    && bdt_pi_v[i]  > PION_BDT_PI
                    && bragg_pi >= pion_bragg_cut
                    && len    > PION_TRK_LEN) {
                    pion_count++;
                    // Select highest-LLR pion (original had unreachable else branch)
                    if (pion_index == -1 || llr > best_pion_llr) {
                        pion_index    = i;
                        best_pion_llr = llr;
                    }
                }
            }

            // ── Derived booleans ──────────────────────────────────────────────
            const bool sel_has_muon = (muon_index != -1);
            bool muon_in_gap = false, pion_in_gap = false;
            bool sel_contained_pion = false;

            if (muon_index != -1)
                muon_in_gap = pfnplanehits_U->at(muon_index) >= 1
                           && pfnplanehits_V->at(muon_index) >= 1
                           && pfnplanehits_Y->at(muon_index) >= 1;

            if (pion_index != -1) {
                pion_in_gap = pfnplanehits_U->at(pion_index) > 0
                           && pfnplanehits_V->at(pion_index) > 0
                           && pfnplanehits_Y->at(pion_index) > 0;
                if (pion_count == 1) sel_contained_pion = true;
            }

            // ── Opening angle ─────────────────────────────────────────────────
            double mu_pi_oa = 0.0;
            bool   oa_ok    = false;
            if (muon_index != -1 && pion_index != -1) {
                const TVector3 muDir(trk_dir_x_v->at(muon_index),
                                     trk_dir_y_v->at(muon_index),
                                     trk_dir_z_v->at(muon_index));
                const TVector3 piDir(trk_dir_x_v->at(pion_index),
                                     trk_dir_y_v->at(pion_index),
                                     trk_dir_z_v->at(pion_index));
                mu_pi_oa = muDir.Angle(piDir);
                // Reco opening-angle cut: θμπ < 2.6 rad (upper bound only), rejecting
                // the large-angle cosmic/dirt background. The same cut is imposed on
                // the signal definition. No lower bound (small angle is high-purity).
                oa_ok    = (mu_pi_oa < OA_CUT);
            }

            // ── Shower cut ────────────────────────────────────────────────────
            bool shower_ok = (nPrimaryShowers == 0);
            if (!shower_ok && nPrimaryShowers == 1 && shower_index != -1) {
                const double dsx = trk_sce_start_x_v->at(shower_index) - rvx;
                const double dsy = trk_sce_start_y_v->at(shower_index) - rvy;
                const double dsz = trk_sce_start_z_v->at(shower_index) - rvz;
                const double sd  = std::sqrt(dsx*dsx + dsy*dsy + dsz*dsz);
                shower_ok = (pfnplanehits_Y->at(shower_index) < 50)
                         && (pfnplanehits_Y->at(shower_index) > 0)
                         && (sd > 10.0);
            }

            // ── Cut array ─────────────────────────────────────────────────────
            bool Cuts[NCUTS];
            // Cosmic veto folded into the NeutrinoSlice stage: reject slices the
            // cosmic BDT flags as cosmic-like (bdt_cosmic <= cut).  Applied to every
            // sample identically, so any residual data/MC modelling difference is
            // absorbed by the efficiency correction.  (The former CRT veto here was
            // inert — crthitpe/crtveto are unfilled in these NuMI ntuples.)
            // A non-positive cosmic_bdt_cut disables the veto (true no-op).
            // bdt_cosmic spans [-1,1], so only a positive threshold rejects
            // anything; otherwise every event passes.
            const bool cosmic_ok = (cosmic_bdt_cut <= 0.0)
                                 ? true
                                 : (bdt_cosmic > cosmic_bdt_cut);
            // Software-trigger requirement.  The data and EXT files are already
            // filtered to swtrig==1 (so this is a no-op for them), but the MC
            // overlay is NOT — only ~87% of MC passes.  Without this the MC
            // prediction is ~13% over-normalised relative to data/EXT.  Applied
            // at NeutrinoSlice (folded in, NCUTS stays 11); it stays OUT of the
            // true-signal denominator (h_true_gen), so swtrig acceptance is
            // correctly absorbed into the efficiency.
            const bool swtrig_ok = (swtrig == 1);
            Cuts[0]  = true;
            Cuts[1]  = Cuts[0] && flash_cut && cosmic_ok && swtrig_ok;
            Cuts[2]  = Cuts[1] && in_fv_reco;
            Cuts[3]  = Cuts[2] && topo_ok;
            Cuts[4]  = Cuts[3] && sel_has_muon;
            Cuts[5]  = Cuts[4] && sel_contained_pion;
            Cuts[6]  = Cuts[5] && muon_in_gap;
            Cuts[7]  = Cuts[6] && pion_in_gap;
            Cuts[8]  = Cuts[7] && shower_ok;
            Cuts[9]  = Cuts[8] && oa_ok;
            Cuts[10] = Cuts[9] && (nonproton < 4) && (nPrimaryTracks < 5);

            // ── Determine sample index ────────────────────────────────────────
            int s_specific = -1;
            if      (i_f < 4 &&  signal) s_specific = 1;
            else if (i_f < 4 && !signal) s_specific = 2;
            else if (i_f == 4)           s_specific = 3;  // EXT
            else if (i_f == 5)           s_specific = 4;  // Dirt
            else if (i_f == 6)           s_specific = 5;  // Data

            // ── Pre-compute all histogram values once ─────────────────────────
            // Stored in local variables to avoid repeated vector indexing in the
            // cut loop below.

            double muon_len=0, muon_pid=0, muon_vtx_dist=0, muon_bdt_mip=0;
            double muon_range_mom=0, muon_mcs_mom=0;
            double muon_theta=0, muon_costheta=0, muon_phi=0;
            double muon_sx=0, muon_sy=0, muon_sz=0;
            int    muon_hU=0, muon_hV=0, muon_hY=0;

            if (muon_index != -1) {
                const int mi = muon_index;
                muon_len      = trk_len_v->at(mi);
                muon_pid      = trk_llr_pid_score_v->at(mi);
                muon_bdt_mip  = bdt_mip_v[mi];
                muon_range_mom= trk_range_muon_mom_v->at(mi);
                muon_mcs_mom  = trk_mcs_muon_mom_v->at(mi);
                muon_hU       = pfnplanehits_U->at(mi);
                muon_hV       = pfnplanehits_V->at(mi);
                muon_hY       = pfnplanehits_Y->at(mi);
                muon_sx       = trk_sce_start_x_v->at(mi);
                muon_sy       = trk_sce_start_y_v->at(mi);
                muon_sz       = trk_sce_start_z_v->at(mi);
                const double dx = muon_sx - rvx, dy = muon_sy - rvy, dz = muon_sz - rvz;
                muon_vtx_dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                const TVector3 mDir(trk_dir_x_v->at(mi), trk_dir_y_v->at(mi), trk_dir_z_v->at(mi));
                muon_theta    = mDir.Theta();
                muon_costheta = mDir.CosTheta();
                muon_phi      = mDir.Phi();
            }

            double pion_len=0, pion_pid=0, pion_vtx_dist=0;
            double pion_bdt_mip=0, pion_bdt_pi=0;
            double pion_range_mom=0, pion_mcs_mom=0;
            double pion_theta=0, pion_costheta=0, pion_phi=0;
            double pion_sx=0, pion_sy=0, pion_sz=0;
            int    pion_hU=0, pion_hV=0, pion_hY=0;

            if (pion_index != -1) {
                const int pi = pion_index;
                pion_len      = trk_len_v->at(pi);
                pion_pid      = trk_llr_pid_score_v->at(pi);
                pion_bdt_mip  = bdt_mip_v[pi];
                pion_bdt_pi   = bdt_pi_v[pi];
                pion_range_mom= trk_range_muon_mom_v->at(pi);
                pion_mcs_mom  = trk_mcs_muon_mom_v->at(pi);
                pion_hU       = pfnplanehits_U->at(pi);
                pion_hV       = pfnplanehits_V->at(pi);
                pion_hY       = pfnplanehits_Y->at(pi);
                pion_sx       = trk_sce_start_x_v->at(pi);
                pion_sy       = trk_sce_start_y_v->at(pi);
                pion_sz       = trk_sce_start_z_v->at(pi);
                const double dx = pion_sx - rvx, dy = pion_sy - rvy, dz = pion_sz - rvz;
                pion_vtx_dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                const TVector3 pDir(trk_dir_x_v->at(pi), trk_dir_y_v->at(pi), trk_dir_z_v->at(pi));
                pion_theta    = pDir.Theta();
                pion_costheta = pDir.CosTheta();
                pion_phi      = pDir.Phi();
            }

            double mu_pi_dist = 0.0;
            if (muon_index != -1 && pion_index != -1) {
                const double dx = muon_sx - pion_sx;
                const double dy = muon_sy - pion_sy;
                const double dz = muon_sz - pion_sz;
                mu_pi_dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            }

            double shr_theta=0, shr_costheta=0, shr_phi=0, shr_vtxd=0;
            int    shr_hU=0, shr_hV=0, shr_hY=0;
            if (shower_index != -1) {
                const int si = shower_index;
                shr_hU = pfnplanehits_U->at(si);
                shr_hV = pfnplanehits_V->at(si);
                shr_hY = pfnplanehits_Y->at(si);
                const double dx = trk_sce_start_x_v->at(si) - rvx;
                const double dy = trk_sce_start_y_v->at(si) - rvy;
                const double dz = trk_sce_start_z_v->at(si) - rvz;
                shr_vtxd = std::sqrt(dx*dx + dy*dy + dz*dz);
                const TVector3 sDir(trk_dir_x_v->at(si), trk_dir_y_v->at(si), trk_dir_z_v->at(si));
                shr_theta    = sDir.Theta();
                shr_costheta = sDir.CosTheta();
                shr_phi      = sDir.Phi();
            }

            // MC truth kinematics (only for signal events in MC files)
            double true_mu_p=0, true_mu_ct=0, true_mu_th=0, true_mu_ph=0;
            double true_pi_p=0, true_pi_ct=0, true_pi_th=0, true_pi_ph=0;
            double nu_energy=0, Q2_val=0, W_val=0;
            bool   fill_mc_truth = (i_f < 4) && signal && (std::abs(nu_pdg) == 14);

            if (fill_mc_truth) {
                const TVector3 nu3(true_nu_px, true_nu_py, true_nu_pz);
                nu_energy = nu3.Mag();
                const TLorentzVector nu4(true_nu_px, true_nu_py, true_nu_pz, nu_energy);
                double muonE = -1;
                for (size_t i_mc = 0; i_mc < mc_pdg->size(); i_mc++) {
                    const TVector3 p3(mc_px->at(i_mc), mc_py->at(i_mc), mc_pz->at(i_mc));
                    const int apdg = std::abs(mc_pdg->at(i_mc));
                    if (apdg == 13) {
                        const TLorentzVector p4(mc_px->at(i_mc), mc_py->at(i_mc), mc_pz->at(i_mc), mc_E->at(i_mc));
                        Q2_val      = -(nu4 - p4).Mag2();
                        muonE       = mc_E->at(i_mc);
                        true_mu_p   = p3.Mag();
                        true_mu_ct  = p3.CosTheta();
                        true_mu_th  = p3.Theta();
                        true_mu_ph  = p3.Phi();
                    } else if (apdg == 211) {
                        true_pi_p   = p3.Mag();
                        true_pi_ct  = p3.CosTheta();
                        true_pi_th  = p3.Theta();
                        true_pi_ph  = p3.Phi();
                    }
                }
                if (muonE > 0) {
                    const double q0  = nu_energy - muonE;
                    double W2        = M_PROTON*M_PROTON + 2*M_PROTON*q0 - Q2_val;
                    if (W2 <= M_PROTON*M_PROTON) {
                        W2 = 2*M_PROTON*2*M_PROTON + 2*2*M_PROTON*q0 - Q2_val;
                    }
                    W_val = std::sqrt(std::max(W2, 0.0));
                }
            }

            // ── Fill histograms (cut loop, two samples per cut) ───────────────
            // Each cut that passes fires one AllMC fill and one sample-specific fill.
            // All values are already computed above — no more vector indexing here.

            for (int c = 0; c < NCUTS; c++) {
                if (!Cuts[c]) continue;

                // Diagnostic per-category cut-flow (one count per passed cut).
                if (event_cat >= 0) SelectedCat[c][event_cat] += evtW;

                // Samples to fill: the prediction stack AllMC (s=0) plus the
                // sample-specific slot.  Beam-on data (s_specific==5) must NOT
                // enter AllMC, which represents the MC+EXT+dirt prediction only.
                const bool is_data = (s_specific == 5);
                const int  fill_samples[2] = { is_data ? s_specific : 0, s_specific };
                const int  n_fill = (s_specific < 0) ? 1 : (is_data ? 1 : 2);

                for (int fi = 0; fi < n_fill; fi++) {
                    const int s = fill_samples[fi];

                    Selected[c][s] += evtW;

                    auto* H = Histos[c][s]; // alias for readability

                    H[0]->Fill(topological_score, evtW);

                    if (muon_index != -1) {
                        H[1] ->Fill(muon_len,       evtW);
                        H[4] ->Fill(muon_pid,       evtW);
                        H[6] ->Fill(muon_vtx_dist,  evtW);
                        H[9] ->Fill(muon_hU,        evtW);
                        H[10]->Fill(muon_hY,        evtW);
                        H[11]->Fill(muon_hV,        evtW);
                        H[23]->Fill(muon_bdt_mip,   evtW);
                        H[26]->Fill(muon_range_mom, evtW);
                        H[27]->Fill(muon_mcs_mom,   evtW);
                        H[34]->Fill(muon_theta,     evtW);
                        H[35]->Fill(muon_costheta,  evtW);
                        H[36]->Fill(muon_phi,       evtW);
                        H[43]->Fill(muon_sx,        evtW);
                        H[44]->Fill(muon_sy,        evtW);
                        H[45]->Fill(muon_sz,        evtW);
                    }

                    if (pion_index != -1) {
                        H[2] ->Fill(pion_len,       evtW);
                        H[5] ->Fill(pion_pid,       evtW);
                        H[7] ->Fill(pion_vtx_dist,  evtW);
                        H[12]->Fill(pion_hU,        evtW);
                        H[13]->Fill(pion_hY,        evtW);
                        H[14]->Fill(pion_hV,        evtW);
                        H[21]->Fill(pion_bdt_mip,   evtW);
                        H[22]->Fill(pion_bdt_pi,    evtW);
                        H[24]->Fill(pion_range_mom, evtW);
                        H[25]->Fill(pion_mcs_mom,   evtW);
                        H[37]->Fill(pion_theta,     evtW);
                        H[38]->Fill(pion_costheta,  evtW);
                        H[39]->Fill(pion_phi,       evtW);
                        H[46]->Fill(pion_sx,        evtW);
                        H[47]->Fill(pion_sy,        evtW);
                        H[48]->Fill(pion_sz,        evtW);
                    }

                    if (muon_index != -1 && pion_index != -1) {
                        H[3] ->Fill(mu_pi_oa,    evtW);
                        H[8] ->Fill(mu_pi_dist,  evtW);
                        H[33]->Fill(mu_pi_dist,  evtW);
                    }

                    if (shower_index != -1) {
                        H[15]->Fill(shr_hU,       evtW);
                        H[16]->Fill(shr_hY,       evtW);
                        H[17]->Fill(shr_hV,       evtW);
                        H[20]->Fill(shr_vtxd,     evtW);
                        H[40]->Fill(shr_theta,    evtW);
                        H[41]->Fill(shr_costheta, evtW);
                        H[42]->Fill(shr_phi,      evtW);
                    }

                    H[18]->Fill(nPrimaryShowers, evtW);
                    H[19]->Fill(nPrimaryTracks,  evtW);
                    H[28]->Fill(mc_muon_momentum, evtW);
                    H[29]->Fill(mc_pion_momentum, evtW);
                    H[30]->Fill(mc_opening_angle, evtW);

                    if (i_f < 4) {
                        // True neutrino energy from the nu_pdg/nu_e branches.
                        // (mc_pdg[0] is NOT reliably the neutrino — it is the muon
                        // ~60% of the time — so the old mc_E[0] logic filled only
                        // ~⅓ of events.)  |nu_pdg|==14 folds ν̄μ with νμ.
                        if (std::abs(nu_pdg) == 14) H[31]->Fill(nu_e, evtW);
                        if (std::abs(nu_pdg) == 12) H[32]->Fill(nu_e, evtW);
                    }

                    if (fill_mc_truth) {
                        H[49]->Fill(true_mu_p,  evtW);
                        H[50]->Fill(true_mu_ct, evtW);
                        H[51]->Fill(true_mu_th, evtW);
                        H[52]->Fill(true_mu_ph, evtW);
                        H[53]->Fill(true_pi_p,  evtW);
                        H[54]->Fill(true_pi_ct, evtW);
                        H[55]->Fill(true_pi_th, evtW);
                        H[56]->Fill(true_pi_ph, evtW);
                        H[57]->Fill(nu_energy,  evtW);
                        H[58]->Fill(Q2_val,     evtW);
                        H[59]->Fill(W_val,      evtW);
                    }

                } // sample loop
            } // cut loop

            // ── Event-display record: every final-cut event except beam-on data
            if (Cuts[NCUTS-1] && i_f != 6) {
                evl_run_v.push_back(run);     evl_sub_v.push_back(sub);
                evl_evt_v.push_back(evt);     evl_fidx_v.push_back((int)i_f);
                evl_entry_v.push_back(ientry); evl_cat_v.push_back(event_cat);
                evl_mi_v.push_back(muon_index); evl_pi_v.push_back(pion_index);
                evl_fname_v.push_back(Files[i_f]);
            }

            // ── Wiener-SVD inputs: fill at the final selection (Cut[10]) ─────
            // Muon momentum: range-based for contained muons (endpoint inside the
            // containment volume), MCS-based for exiting muons.  The pion is
            // required to be contained (Cut[5]) but the muon leg has no such
            // requirement; range momentum for an exiting track is the distance
            // to the detector wall, which is unrelated to the true momentum.
            // Cut[10] implies Cut[5] (contained pion) and Cut[4] (muon), so both
            // candidate indices are valid here.
            if (Cuts[NCUTS-1] && muon_index != -1 && pion_index != -1) {
                const TVector3 muon_end(trk_sce_end_x_v->at(muon_index),
                                        trk_sce_end_y_v->at(muon_index),
                                        trk_sce_end_z_v->at(muon_index));
                // Muon momentum: range-based for contained muons (endpoint inside
                // the containment volume), MCS-based for exiting muons. Range is
                // the accurate estimator for contained tracks; MCS for exiting.
                const double reco_pmu =
                    isContained(muon_end) ? muon_range_mom : muon_mcs_mom;

                // Reco observable values (parallel to true_obs above).
                // Pion momentum: trk_energy_proton_v (proton-hypothesis range KE)
                // empirically tracks true p_π far better than the muon-hypothesis
                // range momentum (which compresses/saturates near ~0.2 GeV/c).
                // Falls back to the range momentum if the branch is unavailable.
                const double reco_ppi =
                    (trk_energy_proton_v && pion_index < (int)trk_energy_proton_v->size())
                    ? trk_energy_proton_v->at(pion_index) : pion_range_mom;
                double reco_obs[N_OBS];
                reco_obs[OBS_PMU]  = reco_pmu;
                reco_obs[OBS_PPI]  = reco_ppi;
                reco_obs[OBS_CTMU] = muon_costheta;
                reco_obs[OBS_CTPI] = pion_costheta;
                reco_obs[OBS_OA]   = mu_pi_oa;          // rad

                // Sample slot: 0=sig 1=bkg 2=ext 3=dirt 4=data
                int s_reco = -1;
                if      (i_f < 4 &&  signal) s_reco = 0;
                else if (i_f < 4 && !signal) s_reco = 1;
                else if (i_f == 4)           s_reco = 2;
                else if (i_f == 5)           s_reco = 3;
                else if (i_f == 6)           s_reco = 4;

                if (s_reco >= 0)
                    for (int o = 0; o < N_OBS; o++) {
                        h_reco[o][s_reco]->Fill(reco_obs[o], evtW);
                        // Response numerator (true,reco) for selected signal only
                        if (s_reco == 0)
                            h_smear[o]->Fill(true_obs[o], reco_obs[o], evtW);
                    }
            }

            // ── Diagnostic: what is the reco PION CANDIDATE really? ───────────
            // At the final cut, bin selected MC events by the truth PDG that the
            // pion-candidate track backtracks to.  A large proton (2212) fraction
            // here is the smoking gun for proton→pion contamination.
            if (Cuts[NCUTS-1] && i_f < 4 && pion_index >= 0 && backtracked_pdg
                && pion_index < (int)backtracked_pdg->size()) {
                const int bt = std::abs(backtracked_pdg->at(pion_index));
                pionTruthW[bt] += evtW;

                // Fill PID-variable ROC histograms for the two classes we care
                // about: true π± (keep) vs true proton (reject).
                if (bt == 211 || bt == 2212) {
                    const bool isp = (bt == 2212);
                    auto fillv = [&](TH1D* hpi, TH1D* hpr, std::vector<float>* v) {
                        if (!v || pion_index >= (int)v->size()) return;
                        (isp ? hpr : hpi)->Fill(v->at(pion_index), evtW);
                    };
                    fillv(h_chipr_pi, h_chipr_pr, trk_pid_chipr_v);
                    fillv(h_chipi_pi, h_chipi_pr, trk_pid_chipi_v);
                    fillv(h_chimu_pi, h_chimu_pr, trk_pid_chimu_v);
                    fillv(h_pida_pi,  h_pida_pr,  trk_pida_v);
                    fillv(h_bpion_pi, h_bpion_pr, trk_bragg_pion_v);
                }
            }

            // ── Diagnostic: cosmic-rejection variables (event level) ──────────
            // Filled once the event has a muon candidate (Cuts[4]), where cosmic
            // contamination is still substantial.  Signal events = the "keep"
            // class; EXT (beam-off) = the "reject" class.
            if (Cuts[4]) {
                if (event_cat == CAT_SIGNAL) {
                    h_flash_nu->Fill(nu_flashmatch_score, evtW);
                    h_crtpe_nu->Fill(crthitpe,            evtW);
                    h_cosip_nu->Fill(CosmicIP,            evtW);
                    h_bdtc_nu ->Fill(bdt_cosmic,          evtW);
                } else if (i_f == 4) {
                    h_flash_ext->Fill(nu_flashmatch_score, evtW);
                    h_crtpe_ext->Fill(crthitpe,            evtW);
                    h_cosip_ext->Fill(CosmicIP,            evtW);
                    h_bdtc_ext ->Fill(bdt_cosmic,          evtW);
                }
            }

        } // event loop

        delete f;

    } // file loop

    // ── Write histograms.root — cut distributions ─────────────────────────────

    TFile* myfile = new TFile("histograms.root", "RECREATE");
    for (int c = 0; c < NCUTS; c++)
        for (int v = 0; v < NVARS; v++)
            for (int s = 0; s < NSAMPLES; s++)
                Histos[c][s][v]->Write();
    myfile->Close();
    std::cout << "  Wrote: histograms.root\n";

    // ── Write unfolding_inputs.root — consumed by ccpi_xsec.C ────────────────
    // All histograms are written with their Sumw2 error arrays intact so that
    // ccpi_xsec.C can propagate statistical uncertainties correctly.

    TFile* f_unf = new TFile("unfolding_inputs.root", "RECREATE");
    // One set of histograms per observable, names suffixed with OBS_KEY[o]
    // (e.g. h_smear_sel_pmu, h_true_gen_costhpi, h_reco_data_thmupi).
    for (int o = 0; o < N_OBS; o++) {
        h_smear[o]  ->Write();
        h_truegen[o]->Write();
        for (int s = 0; s < 5; s++) h_reco[o][s]->Write();
    }
    f_unf->Close();
    std::cout << "  Wrote: unfolding_inputs.root (" << N_OBS << " observables)\n";

    // ── Write phase1_diagnostics.root — purity-study ROC histograms ──────────
    TFile* f_ph1 = new TFile("phase1_diagnostics.root", "RECREATE");
    for (TH1D* h : {h_chipr_pi,h_chipr_pr, h_chipi_pi,h_chipi_pr,
                    h_chimu_pi,h_chimu_pr, h_pida_pi,h_pida_pr,
                    h_bpion_pi,h_bpion_pr,
                    h_flash_nu,h_flash_ext, h_crtpe_nu,h_crtpe_ext,
                    h_cosip_nu,h_cosip_ext, h_bdtc_nu,h_bdtc_ext})
        h->Write();
    f_ph1->Close();
    std::cout << "  Wrote: phase1_diagnostics.root\n";

    // ── Write event_display_list.root — consumed by selection/event_display.C ─
    {
        TFile* f_evl = new TFile("event_display_list.root", "RECREATE");
        TTree* evtree = new TTree("evdisp", "Selected events for 2D display");
        Int_t    e_run, e_sub, e_evt, e_fidx, e_cat, e_mi, e_pi;
        Long64_t e_entry;
        std::string e_fname;
        evtree->Branch("run",        &e_run);
        evtree->Branch("sub",        &e_sub);
        evtree->Branch("evt",        &e_evt);
        evtree->Branch("file_index", &e_fidx);
        evtree->Branch("entry",      &e_entry);
        evtree->Branch("category",   &e_cat);     // BkgCat enum value
        evtree->Branch("muon_index", &e_mi);
        evtree->Branch("pion_index", &e_pi);
        evtree->Branch("filename",   &e_fname);
        for (size_t k = 0; k < evl_run_v.size(); k++) {
            e_run = evl_run_v[k]; e_sub = evl_sub_v[k]; e_evt = evl_evt_v[k];
            e_fidx = evl_fidx_v[k]; e_entry = evl_entry_v[k]; e_cat = evl_cat_v[k];
            e_mi = evl_mi_v[k]; e_pi = evl_pi_v[k]; e_fname = evl_fname_v[k];
            evtree->Fill();
        }
        evtree->Write();
        f_evl->Close();
        std::cout << "  Wrote: event_display_list.root (" << evl_run_v.size()
                  << " events)\n";
    }

    // ── Print cut table ───────────────────────────────────────────────────────

    std::cout << std::setprecision(1) << std::fixed;
    std::cout << "\nCut name                       Signal      BG          EXT         Dirt        Data        Eff(%)      Pur(%)      E*P\n";
    std::cout << std::string(115, '-') << "\n";
    for (int c = 0; c < NCUTS; c++) {
        const double eff = Selected[c][1] / std::max(Selected[0][1], 1e-9) * 100.0;
        const double pur = Selected[c][1] / std::max(Selected[c][1]+Selected[c][2]+Selected[c][3]+Selected[c][4], 1e-9) * 100.0;
        std::cout << std::left << std::setw(26) << CutsName[c]
                  << std::right
                  << std::setw(12) << Selected[c][1]
                  << std::setw(12) << Selected[c][2]
                  << std::setw(12) << Selected[c][3]
                  << std::setw(12) << Selected[c][4]
                  << std::setw(12) << Selected[c][5]
                  << std::setw(12) << eff
                  << std::setw(12) << pur
                  << std::setw(12) << eff * pur
                  << "\n";
    }

    // ── Diagnostic table 1: cut-flow broken down by truth category ────────────
    // Shows the composition of every stage so we can see which backgrounds
    // dominate.  The neutrino-MC categories (Signal..CCoth) sum to the AllMC
    // prediction; EXT/Dirt/Data are shown alongside for context.
    std::cout << "\n── Truth-category cut-flow (POT-weighted events) ──\n";
    std::cout << std::left << std::setw(20) << "Cut name";
    for (int k = 0; k < NCATS; k++) std::cout << std::right << std::setw(9) << CatName[k];
    std::cout << "\n" << std::string(20 + 9*NCATS, '-') << "\n";
    for (int c = 0; c < NCUTS; c++) {
        std::cout << std::left << std::setw(20) << CutsName[c] << std::right;
        for (int k = 0; k < NCATS; k++) std::cout << std::setw(9) << SelectedCat[c][k];
        std::cout << "\n";
    }

    // ── Diagnostic table 2: truth PDG of the reco pion candidate (final cut) ──
    // A large proton (2212) row quantifies proton→pion misidentification, the
    // most likely purity sink for this topology.
    {
        double tot = 0.0;
        for (auto& kv : pionTruthW) tot += kv.second;
        std::cout << "\n── Reco pion-candidate truth composition at final cut "
                     "(MC, POT-weighted) ──\n";
        std::cout << std::left << std::setw(14) << "|PDG|"
                  << std::setw(14) << "Events"
                  << std::setw(10) << "Frac(%)" << "\n";
        std::cout << std::string(38, '-') << "\n";
        for (auto& kv : pionTruthW)
            std::cout << std::left << std::setw(14) << kv.first
                      << std::setw(14) << kv.second
                      << std::setw(10) << (tot > 0 ? 100.0*kv.second/tot : 0.0) << "\n";
        std::cout << std::left << std::setw(14) << "TOTAL"
                  << std::setw(14) << tot << "\n";
        std::cout << "  (211=π±  2212=proton  13=μ  321=K  11=e/γ  else=other)\n";
    }

    // ── Plots ─────────────────────────────────────────────────────────────────

    TString directory = "plots/";
    TCanvas* c1 = new TCanvas("c1", "", 2000, 800);
  /*  c1->Divide(5, 2);

    TLegend* l1 = new TLegend(0.1, 0.7, 0.3, 0.9);
    l1->SetHeader("Samples");
    for (int s = 0; s < NSAMPLES; s++)
        l1->AddEntry(Histos[0][s][0], Sample[s]);

    for (int i = 0; i < NVARS; i++) {
        for (int c = 0; c < NCUTS; c++) {
            c1->cd(c + 1);
            Histos[c][0][i]->Draw();
            Histos[c][1][i]->Draw("same");
            Histos[c][2][i]->Draw("same");
            Histos[c][3][i]->Draw("same");
            Histos[c][4][i]->Draw("same");
            if (c + 1 == NCUTS) l1->Draw();
        }
        c1->Print(directory + Variable[i] + ".png");
    }

    // Efficiency plots for MC truth variables
    for (int i = 28; i < 30; i++) {
        for (int c = 1; c < NCUTS; c++) {
            c1->cd(c + 1);
            Histos[c][1][i]->Divide(Histos[0][1][i]);
            Histos[c][1][i]->Draw();
        }
        c1->Print(directory + "Eff" + Variable[i] + ".png");
    }
    for (int i = 49; i < 60; i++) {
        for (int c = 1; c < NCUTS; c++) {
            c1->cd(c + 1);
            Histos[c][1][i]->Divide(Histos[0][1][i]);
            Histos[c][1][i]->Draw();
        }
        c1->Print(directory + "Eff" + Variable[i] + ".png");
    }
*/

    // ── Cut-flow plots ─────────────────────────────────────────────────────────
    // Two diagnostic figures drawn in the same style as the unfolding step plots
    // (xsec_analyzer UnfolderNuMI plot_step*.pdf) so the full chain —
    // selection cut-flow → reco spectrum → background subtraction → smearceptance
    // → efficiency → unfolded d#sigma/dx — can be laid out side-by-side:
    //   plots/cutflow_stacked.pdf  POT-weighted events surviving each cut, MC
    //                              stacked by truth category with data overlaid
    //   plots/cutflow_eff_pur.pdf  signal efficiency and purity vs cut stage
    // Both read the SelectedCat[cut][category] table built in the event loop.
    {
        gSystem->mkdir( "plots", kTRUE );
        gStyle->SetOptStat( 0 );
        gStyle->SetLegendBorderSize( 0 );

        const std::string pot_label = "3.28 #times 10^{20} POT";

        // Shared MicroBooNE NuMI + POT annotation (normalized coordinates)
        auto draw_header = []( double x, double y, const std::string& pot ) {
            TLatex l;
            l.SetNDC();
            l.SetTextAlign( 12 );
            l.SetTextSize( 0.035 );
            l.DrawLatex( x, y, "MicroBooNE NuMI" );
            l.DrawLatex( x, y - 0.045, pot.c_str() );
        };

        // Truth categories to stack (bottom → top), paired with a fill colour.
        // Data (CAT_DATA) is drawn separately as points.  The neutrino-MC
        // categories sum to the AllMC prediction; EXT and Dirt are stacked on top.
        struct CFCat { int cat; int color; };
        const std::vector<CFCat> stack_cats = {
            { CAT_SIGNAL,     kRed - 4    },
            { CAT_CC_OTHER,   kOrange + 7 },
            { CAT_CC_0PI,     kOrange - 2 },
            { CAT_CC_PI0,     kYellow - 6 },
            { CAT_CC_MULTIPI, kSpring - 6 },
            { CAT_CC_KAON,    kGreen + 2  },
            { CAT_NC,         kCyan - 3   },
            { CAT_NUE_OTHER,  kAzure + 1  },
            { CAT_OOFV,       kBlue - 7   },
            { CAT_DIRT,       kMagenta - 7},
            { CAT_EXT,        kGray + 1   },
        };

        // ── Plot 1: stacked cut-flow ──────────────────────────────────────────
        {
            TCanvas* cf = new TCanvas( "c_cutflow", "cutflow", 1400, 800 );
            cf->SetRightMargin( 0.27 );
            cf->SetLeftMargin( 0.1 );
            cf->SetBottomMargin( 0.18 );
            cf->SetLogy();

            THStack* hs = new THStack( "hs_cutflow", "" );
            TLegend* lg = new TLegend( 0.74, 0.18, 0.99, 0.9 );

            // Build one histogram per stacked category over the NCUTS stages.
            std::vector<TH1D*> cat_hists;
            for ( const auto& cf_cat : stack_cats ) {
                TH1D* h = new TH1D( Form( "h_cutflow_%d", cf_cat.cat ),
                    "", NCUTS, 0, NCUTS );
                for ( int c = 0; c < NCUTS; ++c ) {
                    h->SetBinContent( c + 1, SelectedCat[c][ cf_cat.cat ] );
                    h->GetXaxis()->SetBinLabel( c + 1, CutsName[c] );
                }
                h->SetFillColor( cf_cat.color );
                h->SetLineColor( kBlack );
                h->SetLineWidth( 1 );
                cat_hists.push_back( h );
                hs->Add( h );
            }

            // Data points (beam-on) overlaid
            TH1D* h_data = new TH1D( "h_cutflow_data", "", NCUTS, 0, NCUTS );
            for ( int c = 0; c < NCUTS; ++c ) {
                double d = SelectedCat[c][ CAT_DATA ];
                h_data->SetBinContent( c + 1, d );
                h_data->SetBinError( c + 1, std::sqrt( std::max( 0., d ) ) );
            }
            h_data->SetMarkerStyle( kFullCircle );
            h_data->SetMarkerSize( 0.8 );
            h_data->SetLineColor( kBlack );
            h_data->SetLineWidth( 2 );

            hs->Draw( "hist" );
            hs->GetXaxis()->SetTitle( "" );
            hs->GetYaxis()->SetTitle( "POT-weighted events" );
            hs->GetYaxis()->SetTitleOffset( 1.2 );
            hs->GetXaxis()->LabelsOption( "v" );
            hs->SetMinimum( 1.0 );
            h_data->Draw( "e1 same" );

            // Legend in top-to-bottom (top of stack first) order, then data
            lg->AddEntry( h_data, "NuMI data", "lep" );
            for ( auto it = cat_hists.rbegin(); it != cat_hists.rend(); ++it ) {
                int cat = stack_cats[ cat_hists.rend() - it - 1 ].cat;
                lg->AddEntry( *it, CatName[ cat ], "f" );
            }
            lg->Draw();

            draw_header( 0.12, 0.86, pot_label );
            cf->SaveAs( "plots/cutflow_stacked.pdf" );
            delete cf;
        }

        // ── Plot 2: signal efficiency and purity vs cut stage ─────────────────
        {
            TH1D* h_eff = new TH1D( "h_cutflow_eff", "", NCUTS, 0, NCUTS );
            TH1D* h_pur = new TH1D( "h_cutflow_pur", "", NCUTS, 0, NCUTS );

            const double sig0 = std::max( SelectedCat[0][ CAT_SIGNAL ], 1e-9 );
            for ( int c = 0; c < NCUTS; ++c ) {
                double sig = SelectedCat[c][ CAT_SIGNAL ];
                // Total predicted selected events = every category except data
                double tot = 0.;
                for ( int k = 0; k < NCATS; ++k )
                    if ( k != CAT_DATA ) tot += SelectedCat[c][k];

                h_eff->SetBinContent( c + 1, 100.0 * sig / sig0 );
                h_pur->SetBinContent( c + 1, tot > 0. ? 100.0 * sig / tot : 0. );
                h_eff->GetXaxis()->SetBinLabel( c + 1, CutsName[c] );
                h_pur->GetXaxis()->SetBinLabel( c + 1, CutsName[c] );
            }

            h_eff->SetLineColor( kAzure - 3 );
            h_eff->SetLineWidth( 3 );
            h_eff->SetMarkerColor( kAzure - 3 );
            h_eff->SetMarkerStyle( kFullCircle );
            h_pur->SetLineColor( kRed - 4 );
            h_pur->SetLineWidth( 3 );
            h_pur->SetMarkerColor( kRed - 4 );
            h_pur->SetMarkerStyle( kFullSquare );

            TCanvas* cep = new TCanvas( "c_cutflow_ep", "cutflow_ep", 1400, 800 );
            cep->SetRightMargin( 0.05 );
            cep->SetLeftMargin( 0.1 );
            cep->SetBottomMargin( 0.18 );

            h_eff->GetYaxis()->SetRangeUser( 0., 105. );
            h_eff->GetYaxis()->SetTitle( "Signal efficiency / purity [%]" );
            h_eff->GetYaxis()->SetTitleOffset( 1.1 );
            h_eff->GetXaxis()->LabelsOption( "v" );
            h_eff->Draw( "hist p l" );
            h_pur->Draw( "hist p l same" );

            TLegend* lg = new TLegend( 0.62, 0.78, 0.94, 0.9 );
            lg->AddEntry( h_eff, "Signal efficiency", "lp" );
            lg->AddEntry( h_pur, "Purity", "lp" );
            lg->Draw();

            draw_header( 0.13, 0.86, pot_label );
            cep->SaveAs( "plots/cutflow_eff_pur.pdf" );
            delete cep;
        }
    }

    delete tmvaMIP;
    delete tmvaPI;
    delete c1;

    std::cout << "\nDone. Run ccpi_xsec.C next to extract the cross section.\n";
}
