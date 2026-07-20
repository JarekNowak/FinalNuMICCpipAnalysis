//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 13 11:49:21 2025 by ROOT version 6.34.08
// from TTree stv_tree/STV analysis tree
// found on file: Run1_fhc_numi.root
//////////////////////////////////////////////////////////

#ifndef stv_tree_h
#define stv_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class stv_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxCC1mu1piXp_nuvertex_contained = 1;

   // Declaration of leaf types
   Bool_t          CC1mu1piXp_Selected;
   Bool_t          CC1mu1piXp_MC_Signal;
   Int_t           CC1mu1piXp_EventCategory;
   Bool_t          CC1mu1piXp_nslice_eq_1;
   Bool_t          CC1mu1piXp_nuvertex_contained_;
   Bool_t          CC1mu1piXp_topo_cut_passed;
   Bool_t          CC1mu1piXp_sig_truevertex_in_fv;
   Bool_t          CC1mu1piXp_sig_ccnc;
   Bool_t          CC1mu1piXp_sig_is_numu;
   Bool_t          CC1mu1piXp_sig_one_muon_above_thresh;
   Bool_t          CC1mu1piXp_sig_one_proton_above_thresh;
   Bool_t          CC1mu1piXp_sig_no_pions;
   Bool_t          CC1mu1piXp_sig_no_heavy_mesons;
   Bool_t          CC1mu1piXp_CandidateMuonTrackEndContainment;
   Int_t           CC1mu1piXp_mc_n_threshold_muon;
   Int_t           CC1mu1piXp_mc_n_threshold_proton;
   Int_t           CC1mu1piXp_mc_n_threshold_pion0;
   Int_t           CC1mu1piXp_mc_n_threshold_pionpm;
   Int_t           CC1mu1piXp_mc_n_heaviermeson;
   Int_t           CC1mu1piXp_CandidateMuonIndex;
   Bool_t          CC1mu1piXp_sig_one_charged_pion;
   Bool_t          CC1mu1piXp_sel_muoncandidate_tracklike;
   Int_t           CC1mu1piXp_CandidatePionIndex;
   Bool_t          is_mc;
   Float_t         spline_weight;
   Float_t         tuned_cv_weight;
   Float_t         ppfx_cv_weight;
   Float_t         normalisation_weight;
   vector<double>  *weight_All_UBGenie;
   vector<double>  *weight_AxFFCCQEshape_UBGenie;
   vector<double>  *weight_DecayAngMEC_UBGenie;
   vector<double>  *weight_NormCCCOH_UBGenie;
   vector<double>  *weight_NormNCCOH_UBGenie;
   vector<double>  *weight_RPA_CCQE_UBGenie;
   vector<double>  *weight_RootinoFix_UBGenie;
   vector<double>  *weight_ThetaDelta2NRad_UBGenie;
   vector<double>  *weight_Theta_Delta2Npi_UBGenie;
   vector<double>  *weight_TunedCentralValue_UBGenie;
   vector<double>  *weight_VecFFCCQEshape_UBGenie;
   vector<double>  *weight_XSecShape_CCMEC_UBGenie;
   vector<double>  *weight_ppfx_all;
   vector<double>  *weight_ppfx_cv_UBPPFXCV;
   vector<double>  *weight_reint_all;
   vector<double>  *weight_splines_general_Spline;
   vector<double>  *weight_xsr_scc_Fa3_SCC;
   vector<double>  *weight_xsr_scc_Fv3_SCC;
   Float_t         nu_completeness_from_pfp;
   Float_t         nu_purity_from_pfp;
   Int_t           nslice;
   Float_t         topological_score;
   Float_t         CosmicIP;
   Float_t         contained_fraction;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           mc_nu_pdg;
   Float_t         mc_nu_vtx_x;
   Float_t         mc_nu_vtx_y;
   Float_t         mc_nu_vtx_z;
   Float_t         mc_nu_energy;
   Int_t           mc_ccnc;
   Int_t           mc_interaction;
   vector<unsigned int> *pfp_generation_v;
   vector<unsigned int> *pfp_trk_daughters_v;
   vector<unsigned int> *pfp_shr_daughters_v;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *shr_start_x_v;
   vector<float>   *shr_start_y_v;
   vector<float>   *shr_start_z_v;
   vector<float>   *shr_dist_v;
   Int_t           shr_id;
   Float_t         shr_score;
   Float_t         shr_energy_cali;
   Float_t         hits_ratio;
   Float_t         shrmoliereavg;
   Float_t         shr_distance;
   Float_t         shr_tkfit_gap10_dedx_Y;
   Float_t         shr_tkfit_2cm_dedx_Y;
   vector<float>   *trk_len_v;
   vector<float>   *trk_sce_start_x_v;
   vector<float>   *trk_sce_start_y_v;
   vector<float>   *trk_sce_start_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_sce_end_x_v;
   vector<float>   *trk_sce_end_y_v;
   vector<float>   *trk_sce_end_z_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_range_muon_mom_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_llr_pid_v;
   vector<float>   *trk_llr_pid_u_v;
   vector<float>   *trk_llr_pid_v_v;
   vector<float>   *trk_llr_pid_y_v;
   vector<float>   *trk_llr_pid_score_v;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   Double_t          candidate_muon_mom_mcs;
   Double_t          candidate_muon_mom_true;
   Double_t          candidate_muon_mom_range;
   Double_t          true_mu_pi_opening_angle;
   Double_t          mu_pi_opening_angle;


   // List of branches
   TBranch        *b_CC1mu1piXp_Selected;   //!
   TBranch        *b_CC1mu1piXp_MC_Signal;   //!
   TBranch        *b_CC1mu1piXp_EventCategory;   //!
   TBranch        *b_CC1mu1piXp_nslice_eq_1;   //!
   TBranch        *b_CC1mu1piXp_nuvertex_contained_;   //!
   TBranch        *b_CC1mu1piXp_topo_cut_passed;   //!
   TBranch        *b_CC1mu1piXp_sig_truevertex_in_fv;   //!
   TBranch        *b_CC1mu1piXp_sig_ccnc;   //!
   TBranch        *b_CC1mu1piXp_sig_is_numu;   //!
   TBranch        *b_CC1mu1piXp_sig_one_muon_above_thresh;   //!
   TBranch        *b_CC1mu1piXp_sig_one_proton_above_thresh;   //!
   TBranch        *b_CC1mu1piXp_sig_no_pions;   //!
   TBranch        *b_CC1mu1piXp_sig_no_heavy_mesons;   //!
   TBranch        *b_CC1mu1piXp_mc_n_threshold_muon;   //!
   TBranch        *b_CC1mu1piXp_mc_n_threshold_proton;   //!
   TBranch        *b_CC1mu1piXp_mc_n_threshold_pion0;   //!
   TBranch        *b_CC1mu1piXp_mc_n_threshold_pionpm;   //!
   TBranch        *b_CC1mu1piXp_mc_n_heaviermeson;   //!
   TBranch        *b_CC1mu1piXp_CandidateMuonIndex;   //!
   TBranch        *b_CC1mu1piXp_sig_one_charged_pion;   //!
   TBranch        *b_CC1mu1piXp_CandidatePionIndex;   //!
   TBranch        *b_is_mc;   //!
   TBranch        *b_spline_weight;   //!
   TBranch        *b_tuned_cv_weight;   //!
   TBranch        *b_ppfx_cv_weight;   //!
   TBranch        *b_normalisation_weight;   //!
   TBranch        *b_weight_All_UBGenie;   //!
   TBranch        *b_weight_AxFFCCQEshape_UBGenie;   //!
   TBranch        *b_weight_DecayAngMEC_UBGenie;   //!
   TBranch        *b_weight_NormCCCOH_UBGenie;   //!
   TBranch        *b_weight_NormNCCOH_UBGenie;   //!
   TBranch        *b_weight_RPA_CCQE_UBGenie;   //!
   TBranch        *b_weight_RootinoFix_UBGenie;   //!
   TBranch        *b_weight_ThetaDelta2NRad_UBGenie;   //!
   TBranch        *b_weight_Theta_Delta2Npi_UBGenie;   //!
   TBranch        *b_weight_TunedCentralValue_UBGenie;   //!
   TBranch        *b_weight_VecFFCCQEshape_UBGenie;   //!
   TBranch        *b_weight_XSecShape_CCMEC_UBGenie;   //!
   TBranch        *b_weight_ppfx_all;   //!
   TBranch        *b_weight_ppfx_cv_UBPPFXCV;   //!
   TBranch        *b_weight_reint_all;   //!
   TBranch        *b_weight_splines_general_Spline;   //!
   TBranch        *b_weight_xsr_scc_Fa3_SCC;   //!
   TBranch        *b_weight_xsr_scc_Fv3_SCC;   //!
   TBranch        *b_nu_completeness_from_pfp;   //!
   TBranch        *b_nu_purity_from_pfp;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_CosmicIP;   //!
   TBranch        *b_contained_fraction;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_mc_nu_pdg;   //!
   TBranch        *b_mc_nu_vtx_x;   //!
   TBranch        *b_mc_nu_vtx_y;   //!
   TBranch        *b_mc_nu_vtx_z;   //!
   TBranch        *b_mc_nu_energy;   //!
   TBranch        *b_mc_ccnc;   //!
   TBranch        *b_mc_interaction;   //!
   TBranch        *b_pfp_generation_v;   //!
   TBranch        *b_pfp_trk_daughters_v;   //!
   TBranch        *b_pfp_shr_daughters_v;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_shr_start_x_v;   //!
   TBranch        *b_shr_start_y_v;   //!
   TBranch        *b_shr_start_z_v;   //!
   TBranch        *b_shr_dist_v;   //!
   TBranch        *b_shr_id;   //!
   TBranch        *b_shr_score;   //!
   TBranch        *b_shr_energy_cali;   //!
   TBranch        *b_hits_ratio;   //!
   TBranch        *b_shrmoliereavg;   //!
   TBranch        *b_shr_distance;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_Y;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_Y;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_sce_start_x_v;   //!
   TBranch        *b_trk_sce_start_y_v;   //!
   TBranch        *b_trk_sce_start_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_sce_end_x_v;   //!
   TBranch        *b_trk_sce_end_y_v;   //!
   TBranch        *b_trk_sce_end_z_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_range_muon_mom_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_llr_pid_v;   //!
   TBranch        *b_trk_llr_pid_u_v;   //!
   TBranch        *b_trk_llr_pid_v_v;   //!
   TBranch        *b_trk_llr_pid_y_v;   //!
   TBranch        *b_trk_llr_pid_score_v;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *CC1mu1piXp_candidate_muon_mom_mcs;
   TBranch        *CC1mu1piXp_candidate_muon_mom_true;
   TBranch        *CC1mu1piXp_candidate_muon_mom_range;
   TBranch        *b_CC1mu1piXp_CandidateMuonTrackEndContainment;
   TBranch        *CC1mu1piXp_true_mu_pi_opening_angle;
   TBranch        *CC1mu1piXp_mu_pi_opening_angle;






   stv_tree(TTree *tree=0);
   virtual ~stv_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef stv_tree_cxx
stv_tree::stv_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Run1_all_cuts_bdt_mumom_all_branches.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Run1_all_cuts_bdt_mumom_all_branches.root");
      }
      f->GetObject("stv_tree",tree);

   }
   Init(tree);
}

stv_tree::~stv_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t stv_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t stv_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void stv_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   weight_All_UBGenie = 0;
   weight_AxFFCCQEshape_UBGenie = 0;
   weight_DecayAngMEC_UBGenie = 0;
   weight_NormCCCOH_UBGenie = 0;
   weight_NormNCCOH_UBGenie = 0;
   weight_RPA_CCQE_UBGenie = 0;
   weight_RootinoFix_UBGenie = 0;
   weight_ThetaDelta2NRad_UBGenie = 0;
   weight_Theta_Delta2Npi_UBGenie = 0;
   weight_TunedCentralValue_UBGenie = 0;
   weight_VecFFCCQEshape_UBGenie = 0;
   weight_XSecShape_CCMEC_UBGenie = 0;
   weight_ppfx_all = 0;
   weight_ppfx_cv_UBPPFXCV = 0;
   weight_reint_all = 0;
   weight_splines_general_Spline = 0;
   weight_xsr_scc_Fa3_SCC = 0;
   weight_xsr_scc_Fv3_SCC = 0;
   pfp_generation_v = 0;
   pfp_trk_daughters_v = 0;
   pfp_shr_daughters_v = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   shr_start_x_v = 0;
   shr_start_y_v = 0;
   shr_start_z_v = 0;
   shr_dist_v = 0;
   trk_len_v = 0;
   trk_sce_start_x_v = 0;
   trk_sce_start_y_v = 0;
   trk_sce_start_z_v = 0;
   trk_distance_v = 0;
   trk_sce_end_x_v = 0;
   trk_sce_end_y_v = 0;
   trk_sce_end_z_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_energy_proton_v = 0;
   trk_range_muon_mom_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_pid_chipr_v = 0;
   trk_llr_pid_v = 0;
   trk_llr_pid_u_v = 0;
   trk_llr_pid_v_v = 0;
   trk_llr_pid_y_v = 0;
   trk_llr_pid_score_v = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("CC1mu1piXp_Selected", &CC1mu1piXp_Selected, &b_CC1mu1piXp_Selected);
   fChain->SetBranchAddress("CC1mu1piXp_MC_Signal", &CC1mu1piXp_MC_Signal, &b_CC1mu1piXp_MC_Signal);
   fChain->SetBranchAddress("CC1mu1piXp_EventCategory", &CC1mu1piXp_EventCategory, &b_CC1mu1piXp_EventCategory);
   fChain->SetBranchAddress("CC1mu1piXp_nslice_eq_1", &CC1mu1piXp_nslice_eq_1, &b_CC1mu1piXp_nslice_eq_1);
   fChain->SetBranchAddress("CC1mu1piXp_nuvertex_contained_", &CC1mu1piXp_nuvertex_contained_, &b_CC1mu1piXp_nuvertex_contained_);
   fChain->SetBranchAddress("CC1mu1piXp_topo_cut_passed", &CC1mu1piXp_topo_cut_passed, &b_CC1mu1piXp_topo_cut_passed);
   fChain->SetBranchAddress("CC1mu1piXp_sig_truevertex_in_fv", &CC1mu1piXp_sig_truevertex_in_fv, &b_CC1mu1piXp_sig_truevertex_in_fv);
   fChain->SetBranchAddress("CC1mu1piXp_sig_ccnc", &CC1mu1piXp_sig_ccnc, &b_CC1mu1piXp_sig_ccnc);
   fChain->SetBranchAddress("CC1mu1piXp_sig_is_numu", &CC1mu1piXp_sig_is_numu, &b_CC1mu1piXp_sig_is_numu);
   fChain->SetBranchAddress("CC1mu1piXp_sig_one_muon_above_thresh", &CC1mu1piXp_sig_one_muon_above_thresh, &b_CC1mu1piXp_sig_one_muon_above_thresh);
   fChain->SetBranchAddress("CC1mu1piXp_sig_one_proton_above_thresh", &CC1mu1piXp_sig_one_proton_above_thresh, &b_CC1mu1piXp_sig_one_proton_above_thresh);
   fChain->SetBranchAddress("CC1mu1piXp_sig_no_pions", &CC1mu1piXp_sig_no_pions, &b_CC1mu1piXp_sig_no_pions);
   fChain->SetBranchAddress("CC1mu1piXp_sig_no_heavy_mesons", &CC1mu1piXp_sig_no_heavy_mesons, &b_CC1mu1piXp_sig_no_heavy_mesons);
   fChain->SetBranchAddress("CC1mu1piXp_mc_n_threshold_muon", &CC1mu1piXp_mc_n_threshold_muon, &b_CC1mu1piXp_mc_n_threshold_muon);
   fChain->SetBranchAddress("CC1mu1piXp_mc_n_threshold_proton", &CC1mu1piXp_mc_n_threshold_proton, &b_CC1mu1piXp_mc_n_threshold_proton);
   fChain->SetBranchAddress("CC1mu1piXp_mc_n_threshold_pion0", &CC1mu1piXp_mc_n_threshold_pion0, &b_CC1mu1piXp_mc_n_threshold_pion0);
   fChain->SetBranchAddress("CC1mu1piXp_mc_n_threshold_pionpm", &CC1mu1piXp_mc_n_threshold_pionpm, &b_CC1mu1piXp_mc_n_threshold_pionpm);
   fChain->SetBranchAddress("CC1mu1piXp_mc_n_heaviermeson", &CC1mu1piXp_mc_n_heaviermeson, &b_CC1mu1piXp_mc_n_heaviermeson);
   fChain->SetBranchAddress("CC1mu1piXp_CandidateMuonIndex", &CC1mu1piXp_CandidateMuonIndex, &b_CC1mu1piXp_CandidateMuonIndex);
   fChain->SetBranchAddress("CC1mu1piXp_sig_one_charged_pion", &CC1mu1piXp_sig_one_charged_pion, &b_CC1mu1piXp_sig_one_charged_pion);
   fChain->SetBranchAddress("CC1mu1piXp_CandidatePionIndex", &CC1mu1piXp_CandidatePionIndex, &b_CC1mu1piXp_CandidatePionIndex);
   fChain->SetBranchAddress("is_mc", &is_mc, &b_is_mc);
   fChain->SetBranchAddress("spline_weight", &spline_weight, &b_spline_weight);
   fChain->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight, &b_tuned_cv_weight);
   fChain->SetBranchAddress("ppfx_cv_weight", &ppfx_cv_weight, &b_ppfx_cv_weight);
   fChain->SetBranchAddress("normalisation_weight", &normalisation_weight, &b_normalisation_weight);
   fChain->SetBranchAddress("weight_All_UBGenie", &weight_All_UBGenie, &b_weight_All_UBGenie);
   fChain->SetBranchAddress("weight_AxFFCCQEshape_UBGenie", &weight_AxFFCCQEshape_UBGenie, &b_weight_AxFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("weight_DecayAngMEC_UBGenie", &weight_DecayAngMEC_UBGenie, &b_weight_DecayAngMEC_UBGenie);
   fChain->SetBranchAddress("weight_NormCCCOH_UBGenie", &weight_NormCCCOH_UBGenie, &b_weight_NormCCCOH_UBGenie);
   fChain->SetBranchAddress("weight_NormNCCOH_UBGenie", &weight_NormNCCOH_UBGenie, &b_weight_NormNCCOH_UBGenie);
   fChain->SetBranchAddress("weight_RPA_CCQE_UBGenie", &weight_RPA_CCQE_UBGenie, &b_weight_RPA_CCQE_UBGenie);
   fChain->SetBranchAddress("weight_RootinoFix_UBGenie", &weight_RootinoFix_UBGenie, &b_weight_RootinoFix_UBGenie);
   fChain->SetBranchAddress("weight_ThetaDelta2NRad_UBGenie", &weight_ThetaDelta2NRad_UBGenie, &b_weight_ThetaDelta2NRad_UBGenie);
   fChain->SetBranchAddress("weight_Theta_Delta2Npi_UBGenie", &weight_Theta_Delta2Npi_UBGenie, &b_weight_Theta_Delta2Npi_UBGenie);
   fChain->SetBranchAddress("weight_TunedCentralValue_UBGenie", &weight_TunedCentralValue_UBGenie, &b_weight_TunedCentralValue_UBGenie);
   fChain->SetBranchAddress("weight_VecFFCCQEshape_UBGenie", &weight_VecFFCCQEshape_UBGenie, &b_weight_VecFFCCQEshape_UBGenie);
   fChain->SetBranchAddress("weight_XSecShape_CCMEC_UBGenie", &weight_XSecShape_CCMEC_UBGenie, &b_weight_XSecShape_CCMEC_UBGenie);
   fChain->SetBranchAddress("weight_ppfx_all", &weight_ppfx_all, &b_weight_ppfx_all);
   fChain->SetBranchAddress("weight_ppfx_cv_UBPPFXCV", &weight_ppfx_cv_UBPPFXCV, &b_weight_ppfx_cv_UBPPFXCV);
   fChain->SetBranchAddress("weight_reint_all", &weight_reint_all, &b_weight_reint_all);
   fChain->SetBranchAddress("weight_splines_general_Spline", &weight_splines_general_Spline, &b_weight_splines_general_Spline);
   fChain->SetBranchAddress("weight_xsr_scc_Fa3_SCC", &weight_xsr_scc_Fa3_SCC, &b_weight_xsr_scc_Fa3_SCC);
   fChain->SetBranchAddress("weight_xsr_scc_Fv3_SCC", &weight_xsr_scc_Fv3_SCC, &b_weight_xsr_scc_Fv3_SCC);
   fChain->SetBranchAddress("nu_completeness_from_pfp", &nu_completeness_from_pfp, &b_nu_completeness_from_pfp);
   fChain->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp, &b_nu_purity_from_pfp);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("CosmicIP", &CosmicIP, &b_CosmicIP);
   fChain->SetBranchAddress("contained_fraction", &contained_fraction, &b_contained_fraction);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("mc_nu_pdg", &mc_nu_pdg, &b_mc_nu_pdg);
   fChain->SetBranchAddress("mc_nu_vtx_x", &mc_nu_vtx_x, &b_mc_nu_vtx_x);
   fChain->SetBranchAddress("mc_nu_vtx_y", &mc_nu_vtx_y, &b_mc_nu_vtx_y);
   fChain->SetBranchAddress("mc_nu_vtx_z", &mc_nu_vtx_z, &b_mc_nu_vtx_z);
   fChain->SetBranchAddress("mc_nu_energy", &mc_nu_energy, &b_mc_nu_energy);
   fChain->SetBranchAddress("mc_ccnc", &mc_ccnc, &b_mc_ccnc);
   fChain->SetBranchAddress("mc_interaction", &mc_interaction, &b_mc_interaction);
   fChain->SetBranchAddress("pfp_generation_v", &pfp_generation_v, &b_pfp_generation_v);
   fChain->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v, &b_pfp_trk_daughters_v);
   fChain->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v, &b_pfp_shr_daughters_v);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("shr_start_x_v", &shr_start_x_v, &b_shr_start_x_v);
   fChain->SetBranchAddress("shr_start_y_v", &shr_start_y_v, &b_shr_start_y_v);
   fChain->SetBranchAddress("shr_start_z_v", &shr_start_z_v, &b_shr_start_z_v);
   fChain->SetBranchAddress("shr_dist_v", &shr_dist_v, &b_shr_dist_v);
   fChain->SetBranchAddress("shr_id", &shr_id, &b_shr_id);
   fChain->SetBranchAddress("shr_score", &shr_score, &b_shr_score);
   fChain->SetBranchAddress("shr_energy_cali", &shr_energy_cali, &b_shr_energy_cali);
   fChain->SetBranchAddress("hits_ratio", &hits_ratio, &b_hits_ratio);
   fChain->SetBranchAddress("shrmoliereavg", &shrmoliereavg, &b_shrmoliereavg);
   fChain->SetBranchAddress("shr_distance", &shr_distance, &b_shr_distance);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y, &b_shr_tkfit_gap10_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y, &b_shr_tkfit_2cm_dedx_Y);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v, &b_trk_sce_start_x_v);
   fChain->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v, &b_trk_sce_start_y_v);
   fChain->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v, &b_trk_sce_start_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v, &b_trk_sce_end_x_v);
   fChain->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v, &b_trk_sce_end_y_v);
   fChain->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v, &b_trk_sce_end_z_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v, &b_trk_range_muon_mom_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_llr_pid_v", &trk_llr_pid_v, &b_trk_llr_pid_v);
   fChain->SetBranchAddress("trk_llr_pid_u_v", &trk_llr_pid_u_v, &b_trk_llr_pid_u_v);
   fChain->SetBranchAddress("trk_llr_pid_v_v", &trk_llr_pid_v_v, &b_trk_llr_pid_v_v);
   fChain->SetBranchAddress("trk_llr_pid_y_v", &trk_llr_pid_y_v, &b_trk_llr_pid_y_v);
   fChain->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v, &b_trk_llr_pid_score_v);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_mcs",&candidate_muon_mom_mcs, &CC1mu1piXp_candidate_muon_mom_mcs);
   fChain->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_true",&candidate_muon_mom_true, &CC1mu1piXp_candidate_muon_mom_true);
   fChain->SetBranchAddress("CC1mu1piXp_candidate_muon_mom_range",&candidate_muon_mom_range, &CC1mu1piXp_candidate_muon_mom_range);
   fChain->SetBranchAddress("CC1mu1piXp_CandidateMuonTrackEndContainment", &CC1mu1piXp_CandidateMuonTrackEndContainment, 
    &b_CC1mu1piXp_CandidateMuonTrackEndContainment);
   fChain->SetBranchAddress("CC1mu1piXp_true_mu_pi_opening_angle", &true_mu_pi_opening_angle, &CC1mu1piXp_true_mu_pi_opening_angle);
   fChain->SetBranchAddress("CC1mu1piXp_mu_pi_opening_angle", &mu_pi_opening_angle, &CC1mu1piXp_mu_pi_opening_angle);



   Notify();
}

bool stv_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void stv_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t stv_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef stv_tree_cxx
