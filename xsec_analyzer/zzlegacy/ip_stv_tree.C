//#define stv_tree_cxx
#define ip_stv_tree_cxx
#include "TFile.h"
#include <TH1.h>
#include "TLine.h"
#include <TH2.h>
#include "stv_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ip_stv_tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L stv_tree.C
//      root> stv_tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1D *h_mag_diff_mom = new TH1D("h_mag_diff_mom", "Magnitude Difference for Uncontained Muons; Momentum difference GeV ; N Particles",50,-5,5);
   TH1D *h_mag_diff_mom_range = new TH1D("h_mag_diff_mom_range", "Magnitude Difference for Contained Muons; Momentum difference GeV ; N Particles",50,-5,5);

   TH1D *h_uncontained_muons_only = new TH1D("h_uncontained_muons_only", "Momentum Resolution for Uncontained Muons; Momntum Resolution GeV ; N Particles",50,-5,5);
   TH1D *h_contained_muons_only = new TH1D("h_contained_muons_only", "Momentum Resolution for Contained Muons; Momntum Resolution GeV ; N Particles",100,-1.5,1.5);
   TH1D *h_contained_muon_pion = new TH1D("h_contained_muon_pion", "Momentum Resolution for Contained Muons; Momntum Resolution GeV ; N Particles",100,-1.5,1.5);
   TH1D *h_contained_muon_proton = new TH1D("h_contained_muon_proton", "Momentum Resolution for Contained Muons; Momntum Resolution GeV ; N Particles",100,-1.5,1.5);
   TH1D *h_contained_muon_other = new TH1D("h_contained_muon_other", "Momentum Resolution for Contained Muons; Momntum Resolution GeV ; N Particles",100,-1.5,1.5);


   TH1D *h_all_muons = new TH1D("h_all_muons", "Momentum Resolution for All Muons; Momntum Resolution GeV ; N Particles",50,-5,5);
   TH2F *h_res = new TH2F("h_res","Contained muon momentum; reco mu mom ; reco - true mom",100,0,2,100,-4.0,2.0);


   TH1D *h_contained_muons_trk_len = new TH1D("h_contained_muons_trk_len", "Track Length for Contained Muons; Track Length ; N Particles",50,-500,2000);
   TH1D *h_contained_muons_trk_start_x = new TH1D("h_contained_muons_trk_start_x", "Track Start x for Contained Muons; trk_sce_start_x ; N Particles",50,-100,300);
   TH1D *h_contained_muons_trk_start_y = new TH1D("h_contained_muons_trk_start_y", "Track Start y for Contained Muons; trk_sce_start_y ; N Particles",50,-350,350);
   TH1D *h_contained_muons_trk_start_z = new TH1D("h_contained_muons_trk_start_z", "Track Start z for Contained Muons; trk_sce_start_z ; N Particles",50,-500,1500);

   TH1D *h_contained_muons_trk_end_x = new TH1D("h_contained_muons_trk_end_x", "Track End x for Contained Muons; trk_sce_end_x ; N Particles",50,-100,300);
   TH1D *h_contained_muons_trk_end_y = new TH1D("h_contained_muons_trk_end_y", "Track End y for Contained Muons; trk_sce_end_y ; N Particles",50,-350,350);
   TH1D *h_contained_muons_trk_end_z = new TH1D("h_contained_muons_trk_end_z", "Track End z for Contained Muons; trk_sce_end_z ; N Particles",50,-500,1500);



   TH1D *h_contained_muons_trk_len_biased = new TH1D("h_contained_muons_trk_len_biased", "Track Length for Contained Muons; Track Length ; N Particles",50,-500,2000);
   TH1D *h_contained_muons_trk_start_x_biased = new TH1D("h_contained_muons_trk_start_x_biased", "Track Start x for Contained Muons; trk_sce_start_x ; N Particles",50,-100,300);
   TH1D *h_contained_muons_trk_start_y_biased = new TH1D("h_contained_muons_trk_start_y_biased", "Track Start y for Contained Muons; trk_sce_start_y ; N Particles",50,-350,350);
   TH1D *h_contained_muons_trk_start_z_biased = new TH1D("h_contained_muons_trk_start_z_biased", "Track Start z for Contained Muons; trk_sce_start_z ; N Particles",50,-500,1500);

   TH1D *h_contained_muons_trk_end_x_biased = new TH1D("h_contained_muons_trk_end_x_biased", "Track End x for Contained Muons; trk_sce_end_x ; N Particles",50,-100,300);
   TH1D *h_contained_muons_trk_end_y_biased = new TH1D("h_contained_muons_trk_end_y_biased", "Track End y for Contained Muons; trk_sce_end_y ; N Particles",50,-350,350);
   TH1D *h_contained_muons_trk_end_z_biased = new TH1D("h_contained_muons_trk_end_z_biased", "Track End z for Contained Muons; trk_sce_end_z ; N Particles",50,-500,1500);




 



   double selected_events = 0.0;  //these are selected signal events
   double mc_signal_events = 0.0;
   double not_selected_events = 0.0; //these are the rejected/background events
   double sig_selected_events = 0.0;
   double bag_selected_events = 0.0;
   double muon_muon = 0.0;
   double muon_pion = 0.0;
   double muon_proton = 0.0;
   double muon_other = 0.0;
   double muon_total = 0.0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;



      if(CC1mu1piXp_Selected == true){

	h_mag_diff_mom->Fill((candidate_muon_mom_mcs) - (candidate_muon_mom_true));

//	for (size_t i = 0; i < trk_len_v->size(); i++) {
//
   	if((CC1mu1piXp_CandidateMuonIndex != -1) && (CC1mu1piXp_CandidateMuonTrackEndContainment)){

	muon_total++;

	 if((abs(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex)) == 13)) {

                muon_muon++;
	}

	 else if((abs(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex)) == 211))  {

                muon_pion++;
	}

	
	else  if(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex) == 2212)  {

                muon_proton++;

	}

	else {

		muon_other++;
	}
	}

	if((CC1mu1piXp_CandidateMuonIndex != -1) && (CC1mu1piXp_CandidateMuonTrackEndContainment)){

	 if(abs(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex)) == 13){
       
	
	
       // if(CC1mu1piXp_CandidateMuonTrackEndContainment){
        //	h_mag_diff_mom_range->Fill((candidate_muon_mom_range) - (candidate_muon_mom_true)); 	
                h_contained_muons_only->Fill(((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
		h_all_muons->Fill(((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
		h_res->Fill((candidate_muon_mom_range),((candidate_muon_mom_range) - (candidate_muon_mom_true)));

		if((((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true)) > -0.25){

			h_contained_muons_trk_len->Fill(trk_len_v->at(CC1mu1piXp_CandidateMuonIndex)); 
			h_contained_muons_trk_start_x->Fill(trk_sce_start_x_v->at(CC1mu1piXp_CandidateMuonIndex));
			h_contained_muons_trk_start_y->Fill(trk_sce_start_y_v->at(CC1mu1piXp_CandidateMuonIndex));
			h_contained_muons_trk_start_z->Fill(trk_sce_start_z_v->at(CC1mu1piXp_CandidateMuonIndex));
			h_contained_muons_trk_end_x->Fill(trk_sce_end_x_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_end_y->Fill(trk_sce_end_y_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_end_z->Fill(trk_sce_end_z_v->at(CC1mu1piXp_CandidateMuonIndex));

		}

		if((((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true)) < -0.25){

                        h_contained_muons_trk_len_biased->Fill(trk_len_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_start_x_biased->Fill(trk_sce_start_x_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_start_y_biased->Fill(trk_sce_start_y_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_start_z_biased->Fill(trk_sce_start_z_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_end_x_biased->Fill(trk_sce_end_x_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_end_y_biased->Fill(trk_sce_end_y_v->at(CC1mu1piXp_CandidateMuonIndex));
                        h_contained_muons_trk_end_z_biased->Fill(trk_sce_end_z_v->at(CC1mu1piXp_CandidateMuonIndex));

                }




	}

	else if(abs(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex)) == 211){
		 h_contained_muon_pion->Fill(((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
	}

	 else if(backtracked_pdg->at(CC1mu1piXp_CandidateMuonIndex) == 2212){
                 h_contained_muon_proton->Fill(((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
        }
          else {
                 h_contained_muon_other->Fill(((candidate_muon_mom_range) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
        }


	}

	/*else{

		h_uncontained_muons_only->Fill(((candidate_muon_mom_mcs) - (candidate_muon_mom_true))/(candidate_muon_mom_true));
		h_all_muons->Fill(((candidate_muon_mom_mcs) - (candidate_muon_mom_true))/(candidate_muon_mom_true));

	}
	}*/
//	}

        selected_events++;


      		if(CC1mu1piXp_MC_Signal == true){

                	sig_selected_events++;

        	}

		//else{

		//	bag_selected_events++;
		//}
                                          
		else{


        	not_selected_events++;

       		 }

       }

      if(CC1mu1piXp_MC_Signal == true){

        mc_signal_events++;


      }

       if (Cut(ientry) < 0) continue;
   }

         double purity = 0.0;
         double efficiency = 0.0;
         purity = double(sig_selected_events)/double(sig_selected_events + not_selected_events);
         efficiency = double(sig_selected_events)/double(mc_signal_events);
         std::cout<<"This is the purity " <<purity<<std::endl;
         std::cout<<"This is the efficiency " <<efficiency<<std::endl;
         std::cout<<"selected events " << selected_events<<std::endl;
         std::cout<<"not selected events "<< not_selected_events<<std::endl;
         std::cout<<"total signal events "<< mc_signal_events<<std::endl;
         std::cout<<"Selected signal events "<< sig_selected_events<<std::endl;
	 std::cout<<"These are the number of entries "<< nentries <<std::endl;

	 std::cout<<"Muon is a muon "<< muon_muon <<" " <<  " Percentage " << ((muon_muon)/(muon_total))*100 <<std::endl;
	 std::cout<<"Muon is a pion "<< muon_pion <<" " <<  " Percentage " << ((muon_pion)/(muon_total))*100  << std::endl;
	 std::cout<<"Muon is a proton "<< muon_proton <<" " <<  " Percentage " << ((muon_proton)/(muon_total))*100 << std::endl;
   	 std::cout<<"Muon is something else "<< muon_other <<" " <<  " Percentage " << ((muon_other)/(muon_total))*100 << std::endl;
         std::cout<<"Muon total "<< muon_total <<std::endl;

	 std::cout<< std::setprecision(1)<<std::fixed;

	 //std::cout<<"\t \t \t Selected signal  \t" << sig_selected_events << " \t selected bckg \t " << not_selected_events <<
        //"\t   efficiency " << efficiency << " \t purity \t " << purity;

auto a = new TCanvas();
h_mag_diff_mom->SetLineWidth(4);
h_mag_diff_mom->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_mag_diff_mom->SetLineColor(kCyan);
//h_mag_diff_mom->Draw("HIST");

auto b = new TCanvas();
h_mag_diff_mom_range->SetLineWidth(4);
h_mag_diff_mom_range->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_mag_diff_mom_range->SetLineColor(kMagenta);
//h_mag_diff_mom_range->Draw("HIST");

auto c = new TCanvas();
h_uncontained_muons_only->SetLineWidth(4);
h_uncontained_muons_only->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_uncontained_muons_only->SetLineColor(kCyan);
h_uncontained_muons_only->Draw("HIST");


auto d = new TCanvas();
h_contained_muons_only->SetLineWidth(4);
h_contained_muons_only->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_only->SetLineColor(kMagenta);
h_contained_muon_pion->SetLineWidth(4);
h_contained_muon_pion->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muon_pion->SetLineColor(kCyan);
h_contained_muon_proton->SetLineWidth(4);
h_contained_muon_proton->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muon_proton->SetLineColor(kSpring);
h_contained_muon_other->SetLineWidth(4);
h_contained_muon_other->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muon_other->SetLineColor(kBlack);
h_contained_muons_only->Draw("HIST");
h_contained_muon_pion->Draw("HIST, same");
h_contained_muon_proton->Draw("HIST, same");
h_contained_muon_other->Draw("HIST, same");




auto e = new TCanvas();
h_all_muons->SetLineWidth(4);
h_all_muons->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_all_muons->SetLineColor(kViolet+9);
h_all_muons->Draw("HIST");

auto f = new TCanvas();
gStyle->SetOptStat("nemri");
h_res->Draw("Colz");

auto g = new TCanvas();
h_contained_muons_trk_len->SetLineWidth(4);
h_contained_muons_trk_len->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_len->SetLineColor(kCyan);
h_contained_muons_trk_len->Draw("HIST");

auto h = new TCanvas();
h_contained_muons_trk_start_x->SetLineWidth(4);
h_contained_muons_trk_start_x->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_x->SetLineColor(kRed);
h_contained_muons_trk_start_x->Draw("HIST");


auto i = new TCanvas();
h_contained_muons_trk_start_y->SetLineWidth(4);
h_contained_muons_trk_start_y->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_y->SetLineColor(kRed);
h_contained_muons_trk_start_y->Draw("HIST");

auto j = new TCanvas();
h_contained_muons_trk_start_z->SetLineWidth(4);
h_contained_muons_trk_start_z->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_z->SetLineColor(kRed);
h_contained_muons_trk_start_z->Draw("HIST");


auto k = new TCanvas();
h_contained_muons_trk_end_x->SetLineWidth(4);
h_contained_muons_trk_end_x->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_x->SetLineColor(kBlue);
h_contained_muons_trk_end_x->Draw("HIST");


auto l = new TCanvas();
h_contained_muons_trk_end_y->SetLineWidth(4);
h_contained_muons_trk_end_y->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_y->SetLineColor(kBlue);
h_contained_muons_trk_end_y->Draw("HIST");

auto m = new TCanvas();
h_contained_muons_trk_end_z->SetLineWidth(4);
h_contained_muons_trk_end_z->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_z->SetLineColor(kBlue);
h_contained_muons_trk_end_z->Draw("HIST");



auto n = new TCanvas();
h_contained_muons_trk_len_biased->SetLineWidth(4);
h_contained_muons_trk_len_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_len_biased->SetLineColor(kMagenta);
h_contained_muons_trk_len_biased->Draw("HIST");


auto o = new TCanvas();
h_contained_muons_trk_start_x_biased->SetLineWidth(4);
h_contained_muons_trk_start_x_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_x_biased->SetLineColor(kSpring);
h_contained_muons_trk_start_x_biased->Draw("HIST");


auto p = new TCanvas();
h_contained_muons_trk_start_y_biased->SetLineWidth(4);
h_contained_muons_trk_start_y_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_y_biased->SetLineColor(kSpring);
h_contained_muons_trk_start_y_biased->Draw("HIST");

auto q = new TCanvas();
h_contained_muons_trk_start_z_biased->SetLineWidth(4);
h_contained_muons_trk_start_z_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_start_z_biased->SetLineColor(kSpring);
h_contained_muons_trk_start_z_biased->Draw("HIST");

auto r = new TCanvas();
h_contained_muons_trk_end_x_biased->SetLineWidth(4);
h_contained_muons_trk_end_x_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_x_biased->SetLineColor(kOrange+1);
h_contained_muons_trk_end_x_biased->Draw("HIST");


auto s = new TCanvas();
h_contained_muons_trk_end_y_biased->SetLineWidth(4);
h_contained_muons_trk_end_y_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_y_biased->SetLineColor(kOrange+1);
h_contained_muons_trk_end_y_biased->Draw("HIST");

auto t = new TCanvas();
h_contained_muons_trk_end_z_biased->SetLineWidth(4);
h_contained_muons_trk_end_z_biased->Scale((2.0*pow(10,20))/(2.33*pow(10,21)));
h_contained_muons_trk_end_z_biased->SetLineColor(kOrange+1);
h_contained_muons_trk_end_z_biased->Draw("HIST");


















 }                                                    
