#pragma once

// ROOT includes
#include "TF1.h"
#include <TMVA/Reader.h>
#include <TH1D.h>
#include <TH2D.h>

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1piXp : public SelectionBase {

public:

  CC1mu1piXp();

  virtual int categorize_event( AnalysisEvent* event ) override final;
  virtual void compute_reco_observables( AnalysisEvent* event ) override final;
  virtual void compute_true_observables( AnalysisEvent* event ) override final;
  virtual void define_category_map() override final;
  virtual void define_constants() override final;
  virtual void define_output_branches() override final;
  virtual void finalize() override;
  virtual bool define_signal( AnalysisEvent* event ) override final;
  virtual void reset() override final;
  virtual bool selection( AnalysisEvent* event ) override final;

private:

  TMVA::Reader * tmvaReader;
  TMVA::Reader * tmvaReader_mu;
  TMVA::Reader * tmvaReader_pi;

  bool sel_nslice_eq_1_;
  bool sel_nshower_eq_0_;
  bool sel_ntrack_gt_2_;
  bool sel_muoncandidate_tracklike_;
  bool sel_pioncandidate_tracklike_;
  bool sel_protoncandidate_tracklike_;
  bool sel_nuvertex_contained_;
  bool sel_muoncandidate_above_p_thresh;
  bool sel_protoncandidate_above_p_thresh;
  bool sel_muoncandidate_contained;
  bool sel_protoncandidate_contained;
  bool sel_muon_momentum_quality;
  bool sel_no_flipped_tracks_;
  bool sel_proton_cand_passed_LLRCut;
  bool sel_muon_momentum_in_range;
  bool sel_muon_costheta_in_range;
  bool sel_muon_phi_in_range;
  bool sel_proton_momentum_in_range;
  bool sel_proton_costheta_in_range;
  bool sel_proton_phi_in_range;
  bool sel_topo_cut_passed_;
  bool sig_truevertex_in_fv_;
  bool sig_ccnc_;
  bool sig_is_numu_;
  bool sig_one_muon_above_thresh_;
  bool sig_one_proton_above_thresh_;
  bool sig_no_pions_;
  bool sig_no_heavy_mesons_;
  bool sig_one_charged_pion_;
  bool sig_no_kaons_;
  bool sel_pion_contained;
  bool muon_in_gap;
  bool pion_in_gap;
  bool sel_MCVertexInFV;
  bool sig_recovertex_in_fv_;
  bool shower_cut;
  bool opening_angle_cut;
  bool CandidateMuonTrackEndContainment;
  int sig_truevertex_fv = 0;
  int sig_ccnc = 0;
  int sig_numu = 0;
  int sig_one_muon = 0;
  int  muon_candidate_counter = 0;
int signal_counter = 0; 
int selected_counter = 0; 
int selected_signal_counter = 0;
  int sig_mc_n_threshold_muon;
  int sig_mc_n_threshold_proton;
  int sig_mc_n_threshold_pion0;
  int sig_mc_n_threshold_pionpm;
  int sig_mc_n_heaviermeson;
  int sig_mc_n_kaons;

  int CandidateMuonIndex;
  int CandidateProtonIndex;
  int CandidatePionIndex;
  int truemuonindex;
  int pion_number;
  int shower_index;
  int nPrimaryShowers;
  int nPrimaryTracks;
  int nonproton;
  TVector3 TrueCandidateMuonP;
  TVector3 TrueCandidatePionP;


//TMVA stuff
//
	float trk_bragg_p_v_tmva;
	float trk_bragg_mu_v_tmva;
	float trk_bragg_mip_v_tmva;
	float trk_llr_pid_score_v_tmva;
	float trk_len_v_tmva;
	float trk_sce_end_x_v_tmva;
	float trk_sce_end_y_v_tmva;
	float trk_sce_end_z_v_tmva;
	float trk_score_v_tmva;
	// float isContained_tmva;

	float trk_bragg_p_v_tmva_mu;
	float trk_bragg_mu_v_tmva_mu;
	float trk_bragg_mip_v_tmva_mu;
	float trk_llr_pid_score_v_tmva_mu;
	float trk_len_v_tmva_mu;
	float trk_sce_end_x_v_tmva_mu;
	float trk_sce_end_y_v_tmva_mu;
	float trk_sce_end_z_v_tmva_mu;
	float trk_score_v_tmva_mu;

	float trk_bragg_p_v_tmva_pi;
	float trk_bragg_mu_v_tmva_pi;
	float trk_bragg_mip_v_tmva_pi;
	float trk_llr_pid_score_v_tmva_pi;
	float trk_len_v_tmva_pi;
	float trk_sce_end_x_v_tmva_pi;
	float trk_sce_end_y_v_tmva_pi;
	float trk_sce_end_z_v_tmva_pi;
	float trk_score_v_tmva_pi;
        float tmvaOutput_mip;
        float tmvaOutput;
        float tmvaOutput_pi;
  
  


  double Reco_Pt;
  double Reco_Ptx;
  double Reco_Pty;
  double Reco_PL;
  double Reco_Pn;
  double Reco_PnPerp;
  double Reco_PnPerpx;
  double Reco_PnPerpy;
  double Reco_PnPar;
  double Reco_DeltaAlphaT;
  double Reco_DeltaAlpha3Dq;
  double Reco_DeltaAlpha3DMu;
  double BackTrack_DeltaPhiT;
  double BackTrack_DeltaPhi3D;
  double BackTrack_ECal;
  double BackTrack_EQE;
  double BackTrack_Q2;
  double BackTrack_A;
  double BackTrack_EMiss;
  double BackTrack_kMiss;
  double BackTrack_PMiss;
  double BackTrack_PMissMinus;
  double True_Pt;
  double True_Ptx;
  double True_Pty;
  double True_PL;
  double True_Pn;
  double True_PnPerp;
  double True_PnPerpx;
  double True_PnPerpy;
  double True_PnPar;
  double True_DeltaAlphaT;
  double True_DeltaAlpha3Dq;
  double True_DeltaAlpha3DMu;
  double True_DeltaPhiT;
  double True_DeltaPhi3D;
  double True_ECal;
  double True_EQE;
  double True_Q2;
  double True_A;
  double True_EMiss;
  double True_kMiss;
  double True_PMiss;
  double True_PMissMinus;
  double mu_pi_opening_angle;
  double true_mu_pi_opening_angle;
  double candidate_muon_mom_mcs;
  double candidate_muon_mom_true;
  double candidate_muon_mom_range;
  double candidate_muon_mom_reco;
  // Additional differential observables (added for the multi-observable NuWro
  // closure): reco pion momentum (proton-hypothesis range KE, matching the
  // custom pipeline), reco muon/pion cos(theta) w.r.t. detector z, and the
  // corresponding true quantities.
  double candidate_pion_mom_reco;
  double candidate_pion_mom_true;
  double candidate_muon_costh_reco;
  double candidate_muon_costh_true;
  double candidate_pion_costh_reco;
  double candidate_pion_costh_true;

TH1D* h_selected;
TH1D* h_background;
TH1D* h_signal;

TH1D* h_eff_num;   // numerator: selected AND signal
TH1D* h_eff_den;   // denominator: all true signal
TH1D* h_eff; //final eff

TH1D* h_all_signal;
TH1D* h_all_selected;
TH1D* h_sel_signal;
TH1D* h_sel_bkg;
TH1D* h_purity;
TH2D* h_response;
TH2D* h_response_full;
TH1D* h_mu_eff_num;
TH1D* h_mu_eff_den;
TH1D* h_mu_eff;

TH1D* h_pi_eff_num;
TH1D* h_pi_eff_den;
TH1D* h_pi_eff;

// ===============================
// Muon cut study
// ===============================

TH1D* h_mu_cut0;
TH1D* h_mu_cut1_vertex;
TH1D* h_mu_cut2_topology;
TH1D* h_mu_cut3_tracklike;
TH1D* h_mu_cut4_pioncontained;
TH1D* h_mu_cut5_muongap;
TH1D* h_mu_cut6_piongap;
TH1D* h_mu_cut7_shower;
TH1D* h_mu_cut8_opening;
TH1D* h_mu_cut9_final;

// ===============================
// Pion cut study
// ===============================

TH1D* h_pi_cut0;
TH1D* h_pi_cut1_vertex;
TH1D* h_pi_cut2_topology;
TH1D* h_pi_cut3_tracklike;
TH1D* h_pi_cut4_pioncontained;
TH1D* h_pi_cut5_muongap;
TH1D* h_pi_cut6_piongap;
TH1D* h_pi_cut7_shower;
TH1D* h_pi_cut8_opening;
TH1D* h_pi_cut9_final;

  // Per-cut RECO cut-flow, all 5 observables x 10 stages, ALL events (signal+bkg).
  // obs index: 0=pmu 1=ppi 2=costhmu 3=costhpi 4=thmupi ; per-run breakdown comes
  // from processing each run file separately.
  TH1D* h_cf[5][10];

  STVCalcType CalcType;
  TF1* fPP;

  //int truemuonindex;
  int trueprotonindex;
};
                  
                                                        
