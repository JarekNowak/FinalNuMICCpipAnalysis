// Standard library includes
#include <iomanip>
#include <cctype>
#include <iostream>
#include <sstream>

#include <algorithm>
#include <cmath>
#include <exception>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "XSecAnalyzer/PGFPlotsDumpUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"


struct GeneratorInfo {
  std::string path;
  int lineColor;
  int lineStyle;
  int lineWidth;
  std::string name;
  float scaling;
};

// Helper function that dumps a lot of the results to simple text files.
// The events_to_xsec_factor is a constant that converts expected true event
// counts to a total cross section (10^{-39} cm^2 / Ar) via multiplication.
void dump_overall_results( const UnfoldedMeasurement& result,
  const std::map< std::string, std::unique_ptr<TMatrixD> >& unf_cov_matrix_map,
  double events_to_xsec_factor,
  const std::map< std::string, std::unique_ptr< PredictedTrueEvents > >& pred_map )
{
  // Dump the unfolded flux-averaged total cross sections (by converting
  // the units on the unfolded signal event counts)
  TMatrixD unf_signal = *result.unfolded_signal_;
  unf_signal *= events_to_xsec_factor;
  dump_text_column_vector( "unfold_output/vec_table_unfolded_signal.txt", unf_signal );

  // Dump similar tables for each of the theoretical predictions (and the fake
  // data truth if applicable). Note that this function expects that the
  // additional smearing matrix A_C has not been applied to these predictions.
  for ( const auto& gen_pair : pred_map ) {
    std::string gen_short_name = gen_pair.second->name();
    TMatrixD temp_gen = gen_pair.second->get_prediction();
    temp_gen *= events_to_xsec_factor;
    dump_text_column_vector( "unfold_output/vec_table_" + gen_short_name + ".txt",
      temp_gen );
  }

  // No unit conversions are necessary for the unfolding, error propagation,
  // and additional smearing matrices since they are dimensionless
  dump_text_matrix( "unfold_output/mat_table_unfolding.txt", *result.unfolding_matrix_ );
  dump_text_matrix( "unfold_output/mat_table_err_prop.txt", *result.err_prop_matrix_ );
  dump_text_matrix( "unfold_output/mat_table_add_smear.txt", *result.add_smear_matrix_ );

  // Convert units on the covariance matrices one-by-one and dump them
  for ( const auto& cov_pair : unf_cov_matrix_map ) {
    const auto& name = cov_pair.first;
    TMatrixD temp_cov_matrix = *cov_pair.second;
    // Note that we need to square the unit conversion factor for the
    // covariance matrix elements
    temp_cov_matrix *= std::pow( events_to_xsec_factor, 2 );

    std::cout << name << ": " << temp_cov_matrix[0][0] << std::endl;

    dump_text_matrix( "unfold_output/mat_table_cov_" + name + ".txt", temp_cov_matrix );
  }
}

std::string toLatexScientific(double value) {
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2) << value;
  std::string str = stream.str();
  size_t pos = str.find('e');
  if (pos != std::string::npos) {
      str.replace(pos, 2, " #times 10^{");
      str += "}";
  }
  return str;
}

TH2D TMatrixDToTH2D(const TMatrixD & mat, const char* name, const char* title, double xlow, double xup, double ylow, double yup) {
  int nX = mat.GetNrows();
  int nY = mat.GetNcols();
  TH2D hist(name, title, nX, xlow, xup, nY, ylow, yup);
  for (int i = 0; i < nX; ++i) {
      for (int j = 0; j < nY; ++j) {
          hist.SetBinContent(i+1, j+1, mat[i][j]);
      }
  }
  return hist;
}

TH1D* get_generator_hist(const TString& filePath, const unsigned int sl_idx, const float scaling = 1.f )
{
    // Open the file
    TFile* file = new TFile(filePath, "readonly");

    // Check if the file was successfully opened
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return nullptr;
    }

    // Define the plot names
    std::vector<TString> plotNames = {
        "TrueElectronEnergyPlot",
        "TrueElectronCosBetaPlot",
        "TruePionCosBetaPlot",
        "TrueElectronPionOpeningAnglePlot",
        "TrueTotalPlot"
    };

    // special case for combined histogram, all bins
    if (sl_idx == 5) {

      // retrive all histograms
      std::vector<TH1D*> hists_set;

      for (int idx = 0; idx < plotNames.size(); idx++) {
        TH1D* hist = (TH1D*)file->Get(plotNames[idx]);
        hists_set.push_back(hist);
      }

      // construct new histogram from these
      TH1D* hist = new TH1D("", "", 21, 0, 21);

      for (int h_idx = 0; h_idx < hists_set.size(); h_idx++) {

          int n_bins = hists_set[h_idx]->GetNbinsX();

          for (int bin_idx = 1; bin_idx < n_bins + 1; bin_idx++) {  // root counts from 1

            double bin_content = hists_set[h_idx]->GetBinContent(bin_idx);
            double bin_error = hists_set[h_idx]->GetBinError(bin_idx);

            // not general, just testing
            int combined_bin_idx = h_idx*5 + bin_idx;

            // set bin
            hist->SetBinContent(combined_bin_idx, bin_content);
            hist->SetBinError(combined_bin_idx, bin_error);

          }
      }

      return hist;
    }
    else {

      // Check if the index is valid
      if (sl_idx >= plotNames.size()) {
          std::cerr << "Invalid slice index: " << sl_idx << std::endl;
          return nullptr;
      }

      // Get the histogram from the file
      TH1D* hist = (TH1D*)file->Get(plotNames[sl_idx]);

      // Check if the histogram was successfully retrieved
      if (!hist) {
          std::cerr << "Failed to retrieve histogram: " << plotNames[sl_idx] << std::endl;
          return nullptr;
      }

      // Scale the histogram
      hist->Scale(scaling);

      return hist;
    }

    // fail safe
    std::cerr << "Failed to retrieve generator prediction, sl_idx = " << sl_idx << std::endl;
    return nullptr;
}

void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist)
{
    // Copy the histogram contents into a column vector
    int num_bins = mat->GetNcols();
    TMatrixD hist_mat(num_bins, 1);
    for (int r = 0; r < num_bins; ++r)
    {
        hist_mat(r, 0) = hist->GetBinContent(r + 1);
    }

    // Multiply the column vector by the input matrix
    // TODO: add error handling here related to matrix dimensions
    TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult,
                                  hist_mat);

    // Update the input histogram contents with the new values
    for (int r = 0; r < num_bins; ++r)
    {
        double val = hist_mat_transformed(r, 0);
        hist->SetBinContent(r + 1, val);
    }
}

// ---------------------------------------------------------------------------
// Diagnostic plots of the intermediate analysis steps that feed the unfolding,
// drawn in the same style as the unfolded-slice plots (large canvas, no stats
// box, right-hand legend, MicroBooNE NuMI label) so the stages can be laid out
// side-by-side with the final cross section:
//   plot_step1_reco_spectrum.pdf     selected reco spectrum: data vs MC, with
//                                    the MC signal + background stacked
//   plot_step2_bkgd_subtraction.pdf  background-subtracted data vs CV MC signal
//   plot_step3_smearing_matrix.pdf   CV smearceptance (response) matrix,
//                                    P(reco | true) including efficiency
//   plot_step4_efficiency.pdf        selection efficiency vs true bin number
// Everything is shown in the native flat reco/true bin-number space (every
// analysis slice concatenated), matching the binning of the unfolding inputs.
// These complement plot_slice_*.pdf (the final unfolded d#sigma/dx).
// ---------------------------------------------------------------------------
void make_analysis_step_plots( const SystematicsCalculator& syst,
  double total_pot )
{
  gStyle->SetOptStat( 0 );
  gStyle->SetLegendBorderSize( 0 );

  const std::string pot_label = toLatexScientific( total_pot ) + " POT";

  // Draws the shared MicroBooNE NuMI + POT annotation in the top-left of the
  // current pad (normalized coordinates).
  auto draw_header = [&]( double x, double y ) {
    TLatex l;
    l.SetNDC();
    l.SetTextAlign( 12 );
    l.SetTextSize( 0.04 );
    l.DrawLatex( x, y, "MicroBooNE NuMI Data" );
    l.DrawLatex( x, y - 0.05, pot_label.c_str() );
  };

  // Reco-space and true-space central-value quantities used in the steps below.
  MeasuredEvents me = syst.get_measured_events();           // data, bkgd, cov
  auto sig_mc  = syst.get_cv_ordinary_reco_signal();        // CV MC reco signal
  auto bkgd_mc = syst.get_cv_ordinary_reco_bkgd();          // CV reco background
  auto smear   = syst.get_cv_smearceptance_matrix();        // rows reco, cols true

  const int nreco = me.reco_signal_->GetNrows();
  const int ntrue = smear->GetNcols();

  const double right_margin = 0.30;

  // ── Step 1: selected reco spectrum, data vs stacked MC (signal + bkgd) ──────
  {
    auto* h_data = new TH1D( "h_step_data", "", nreco, 0, nreco );
    auto* h_sig  = new TH1D( "h_step_sig",  "", nreco, 0, nreco );
    auto* h_bkg  = new TH1D( "h_step_bkg",  "", nreco, 0, nreco );

    for ( int r = 0; r < nreco; ++r ) {
      // Total measured data = background-subtracted data + subtracted background
      double data_tot = ( *me.reco_signal_ )( r, 0 ) + ( *me.reco_bkgd_ )( r, 0 );
      h_data->SetBinContent( r + 1, data_tot );
      h_data->SetBinError( r + 1, std::sqrt( std::max( 0., data_tot ) ) );
      h_sig->SetBinContent( r + 1, ( *sig_mc )( r, 0 ) );
      h_bkg->SetBinContent( r + 1, ( *bkgd_mc )( r, 0 ) );
    }

    h_bkg->SetFillColor( TColor::GetColor("#E69F00") ); // Okabe-Ito orange (colorblind-safe)
    h_bkg->SetLineColor( kGray + 2 );
    h_sig->SetFillColor( TColor::GetColor("#0072B2") ); // Okabe-Ito blue
    h_sig->SetLineColor( kAzure - 3 );
    h_data->SetMarkerStyle( kFullCircle );
    h_data->SetMarkerSize( 0.7 );
    h_data->SetLineColor( kBlack );
    h_data->SetLineWidth( 2 );

    auto* hs = new THStack( "hs_step1", "" );
    hs->Add( h_bkg );   // background on the bottom of the stack
    hs->Add( h_sig );

    auto* c = new TCanvas( "c_step1", "step1", 1400, 800 );
    c->SetRightMargin( right_margin );
    c->SetLeftMargin( 0.1 );
    c->SetBottomMargin( 0.12 );

    double ymax = std::max( h_data->GetMaximum(),
      h_sig->GetMaximum() + h_bkg->GetMaximum() );
    hs->SetMaximum( ymax * 1.4 );
    hs->Draw( "hist" );
    hs->GetXaxis()->SetTitle( "Reco bin number" );
    hs->GetYaxis()->SetTitle( "Selected events" );
    h_data->Draw( "e1 same" );

    auto* lg = new TLegend( 1 - right_margin + 0.02, 0.55, 0.98, 0.88 );
    lg->AddEntry( h_data, "NuMI data", "lep" );
    lg->AddEntry( h_sig,  "MC signal (CC#pi^{#pm})", "f" );
    lg->AddEntry( h_bkg,  "MC + EXT background", "f" );
    lg->Draw();

    draw_header( 0.13, 0.85 );
    c->SaveAs( "unfold_output/plot_step1_reco_spectrum.pdf" );
  }

  // ── Step 2: background-subtracted data vs CV MC signal ──────────────────────
  {
    auto* h_bsub = new TH1D( "h_step_bsub", "", nreco, 0, nreco );
    auto* h_sig2 = new TH1D( "h_step_sig2", "", nreco, 0, nreco );

    for ( int r = 0; r < nreco; ++r ) {
      h_bsub->SetBinContent( r + 1, ( *me.reco_signal_ )( r, 0 ) );
      double var = ( *me.cov_matrix_ )( r, r );
      h_bsub->SetBinError( r + 1, var > 0. ? std::sqrt( var ) : 0. );
      h_sig2->SetBinContent( r + 1, ( *sig_mc )( r, 0 ) );
    }

    h_sig2->SetLineColor( kAzure - 3 );
    h_sig2->SetLineWidth( 3 );
    h_sig2->SetLineStyle( 5 );
    h_bsub->SetMarkerStyle( kFullCircle );
    h_bsub->SetMarkerSize( 0.7 );
    h_bsub->SetLineColor( kBlack );
    h_bsub->SetLineWidth( 2 );

    auto* c = new TCanvas( "c_step2", "step2", 1400, 800 );
    c->SetRightMargin( right_margin );
    c->SetLeftMargin( 0.1 );
    c->SetBottomMargin( 0.12 );

    double ymax = std::max( h_bsub->GetMaximum(), h_sig2->GetMaximum() );
    h_bsub->GetYaxis()->SetRangeUser( 0., ymax * 1.4 );
    h_bsub->GetXaxis()->SetTitle( "Reco bin number" );
    h_bsub->GetYaxis()->SetTitle( "Background-subtracted events" );
    h_bsub->Draw( "e1" );
    h_sig2->Draw( "hist same" );
    h_bsub->Draw( "e1 same" );  // data on top

    auto* lg = new TLegend( 1 - right_margin + 0.02, 0.6, 0.98, 0.88 );
    lg->AddEntry( h_bsub, "Bkgd-subtracted data", "lep" );
    lg->AddEntry( h_sig2, "MC signal (CC#pi^{#pm})", "l" );
    lg->Draw();

    draw_header( 0.13, 0.85 );
    c->SaveAs( "unfold_output/plot_step2_bkgd_subtraction.pdf" );
  }

  // ── Step 3: CV smearceptance (response) matrix, P(reco | true) ──────────────
  {
    // x-axis = true bin, y-axis = reco bin (migration convention)
    auto* h_sm = new TH2D( "h_step_smear", "", ntrue, 0, ntrue, nreco, 0, nreco );
    for ( int r = 0; r < nreco; ++r )
      for ( int t = 0; t < ntrue; ++t )
        h_sm->SetBinContent( t + 1, r + 1, ( *smear )( r, t ) );

    gStyle->SetPalette( kBird );

    auto* c = new TCanvas( "c_step3", "step3", 1200, 1000 );
    c->SetRightMargin( 0.15 );
    c->SetLeftMargin( 0.1 );
    c->SetBottomMargin( 0.1 );
    h_sm->Draw( "colz" );
    h_sm->GetXaxis()->SetTitle( "True bin number" );
    h_sm->GetYaxis()->SetTitle( "Reco bin number" );
    h_sm->GetZaxis()->SetTitle( "P(reco | true) #times efficiency" );

    draw_header( 0.12, 0.93 );
    c->SaveAs( "unfold_output/plot_step3_smearing_matrix.pdf" );
  }

  // ── Step 4: selection efficiency vs true bin (column sums of smearceptance) ──
  {
    auto* h_eff = new TH1D( "h_step_eff", "", ntrue, 0, ntrue );
    for ( int t = 0; t < ntrue; ++t ) {
      double eff = 0.;
      for ( int r = 0; r < nreco; ++r ) eff += ( *smear )( r, t );
      h_eff->SetBinContent( t + 1, eff );
    }

    h_eff->SetLineColor( kAzure - 3 );
    h_eff->SetLineWidth( 3 );
    h_eff->SetFillColorAlpha( kAzure - 7, 0.35 );

    auto* c = new TCanvas( "c_step4", "step4", 1400, 800 );
    c->SetRightMargin( 0.05 );
    c->SetLeftMargin( 0.1 );
    c->SetBottomMargin( 0.12 );
    h_eff->GetYaxis()->SetRangeUser( 0., std::min( 1., h_eff->GetMaximum() * 1.4 ) );
    h_eff->GetXaxis()->SetTitle( "True bin number" );
    h_eff->GetYaxis()->SetTitle( "Selection efficiency" );
    h_eff->Draw( "hist" );

    draw_header( 0.13, 0.85 );
    c->SaveAs( "unfold_output/plot_step4_efficiency.pdf" );
  }
}

void UnfolderNuMI(std::string XSEC_Config, std::string SLICE_Config, std::string OutputDirectory, std::string OutputFileName) {

  // set to using fake data
  bool using_fake_data = false;
  bool total_only = false;

  std::cout << "\nRunning Unfolder.C with options:" << std::endl;
  std::cout << "\tXSEC_Config: " << XSEC_Config << std::endl;
  std::cout << "\tSLICE_Config: " << SLICE_Config << std::endl;
  std::cout << "\tOutputDirectory: " << OutputDirectory << std::endl;
  std::cout << "\tOutputFileName: " << OutputFileName << std::endl;
  std::cout << "\n" << std::endl;

  // Use a CrossSectionExtractor object to handle the systematics and unfolding
  auto extr = std::make_unique< CrossSectionExtractor >( XSEC_Config );

  // Enable the fake-data closure comparison automatically when the "data"
  // sample carries MC truth (i.e. it is a fake-data overlay). This lets the
  // unfolded result be compared bin-by-bin to the fake data's own truth, which
  // is the closure test proper. Was hardcoded false, so the truth curve never
  // appeared even for fake-data runs.
  using_fake_data = ( extr->get_syst().fake_data_universe() != nullptr );
  std::cout << "using_fake_data = " << using_fake_data << std::endl;

  // get unfolded results
  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  const auto& pred_map = extr->get_prediction_map();
  double total_pot = extr->get_data_pot();

  double A_C_total = 1;
  if (total_only) {
    const TMatrixD &A_C_temp =  *xsec.result_.add_smear_matrix_;
    A_C_total = A_C_temp(0,0);
    std::cout << "Total Only Mode: A_C matrix element = " << A_C_total << std::endl;
  }

  // Truth Generator Plots
  const std::string genPath = "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/generatorFiles/";

  // Format: path, lineColor, lineStyle, lineWidth, name, scaling
  std::vector<GeneratorInfo> generators = {
    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //{genPath + "FlatTreeAnalyzerOutput_NuWro_Combined.root", kTeal+3, 1, 3, "NuWro 21.09.2", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_NEUT_Combined.root", kOrange+8, 1, 3, "NEUT 5.4.0.1", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_GiBUU_Combined.root", kMagenta-3, 1, 3, "GiBUU 2023", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_GiBUU2025_Combined.root", kMagenta-3, 1, 3, "GiBUU 2025", 1.f},
    //{genPath + "FlatTreeAnalyzerOutput_Genie_Combined.root", kRed+1, 1, 3, "GENIE 3.4.2 AR23", 1.f},
  };

  // Was hardcoded to 2: a slice config with more than two slices silently
  // dropped the rest, and one with a single slice threw from slices_.at(1).
  // Unfolder.C already loops over sb.slices_.size().
  for (int sl_idx = 0; sl_idx < static_cast<int>( sb.slices_.size() ); sl_idx++) {

    const auto& slice = sb.slices_.at( sl_idx );

    // Make a histogram showing the unfolded true event counts in the current slice
    SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(
      (*xsec.result_.unfolded_signal_), slice, xsec.result_.cov_matrix_.get() );

    // Temporary copies of the unfolded true event count slices with
    // different covariance matrices
    std::map< std::string, std::unique_ptr<SliceHistogram> > sh_cov_map;
    for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
      const auto& uc_name = uc_pair.first;
      const auto& uc_matrix = uc_pair.second;

      auto& uc_ptr = sh_cov_map[ uc_name ];
      uc_ptr.reset(
        SliceHistogram::make_slice_histogram( (*xsec.result_.unfolded_signal_),
          slice, uc_matrix.get() )
      );
    }

    // Also use the GENIE CV model to do the same
    auto genie_cv_it = pred_map.find("MicroBooNETune");
    if ( genie_cv_it == pred_map.end() ) {
      throw std::runtime_error( "UnfolderNuMI: no prediction named"
        " \"MicroBooNETune\" in the xsec config. The key is the quoted"
        " description field, so it must match exactly." );
    }
    // Apply the additional smearing matrix A_C to every truth-space prediction
    // before comparing to the unfolded data. Wiener-SVD unfolding returns
    // A_C * true_signal, so a model prediction must be smeared by A_C to live
    // in the same space; otherwise sharp truth features (e.g. the forward muon
    // peak) are compared against a regularised measurement and the chi^2 is
    // meaningless. This applies uniformly to the MicroBooNE Tune and to every
    // file-based generator prediction (NuWro, ...).
    const TMatrixD& A_C_smear = *xsec.result_.add_smear_matrix_;
    auto apply_ac = [&A_C_smear]( const TMatrixD& truth ) {
      return TMatrixD( A_C_smear, TMatrixD::kMult, truth );
    };

    TMatrixD genie_cv_truth = apply_ac( genie_cv_it->second->get_prediction() );

    for ( const auto& gen_pair : pred_map ) {
      std::cout << "Key: " << gen_pair.first << std::endl;
    }

    SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(
      genie_cv_truth, slice, nullptr );

    // If present, also use the truth information from the fake data to do the same
    SliceHistogram* slice_truth = nullptr;
    // The fake-data truth curve is only drawn if the config declares a
    // "Fakedata" prediction (Prediction ... univ FakeData). using_fake_data is
    // auto-detected from the sample, so it can be true without that line; in
    // that case just skip the truth curve rather than aborting.
    auto fake_data_truth_it = pred_map.find("Fakedata");
    if ( using_fake_data && fake_data_truth_it == pred_map.end() ) {
      std::cout << "[UnfolderNuMI] fake-data sample detected but no \"Fakedata\""
        " prediction in the config; skipping the truth-comparison curve.\n";
    }
    if ( using_fake_data && fake_data_truth_it != pred_map.end() ) {
      TMatrixD fake_data_truth = fake_data_truth_it->second->get_prediction();
      TMatrixD fake_data_truth_cov(fake_data_truth.GetNrows(), fake_data_truth.GetNrows());

      for (int b = 0; b < fake_data_truth.GetNrows(); ++b)
      {
        // Was fake_data_truth(0,0): every diagonal element got bin 0's content,
        // so every truth bin was assigned the same variance and the closure
        // chi^2 below was wrong for all fake-data studies.
        fake_data_truth_cov(b,b) = fake_data_truth(b,0);
      }

      slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth,
        slice, &fake_data_truth_cov );
    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
    auto& slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map[ "unfolded data" ] = slice_unf;
    if ( slice_truth ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

    // Add every other prediction in the config (e.g. file-based generator
    // predictions like NuWro, GiBUU, NEUT) to the overlay, treated exactly like
    // the MicroBooNE Tune: same trans_mat conversion below, same chi^2 vs the
    // unfolded data. "MicroBooNE Tune" and "Fakedata" are already handled.
    for ( const auto& pred_pair : pred_map ) {
      const std::string& desc = pred_pair.first;
      // Note: the config parser strips whitespace from the quoted description
      // (it reads char-by-char with `>>`, which skips spaces), so the CV key is
      // "MicroBooNETune", not "MicroBooNE Tune". It is already in slice_gen_map
      // as slice_cv; skip it (and the fake-data truth) to avoid a duplicate.
      if ( desc == "MicroBooNETune" || desc == "Fakedata" ) continue;
      if ( slice_gen_map.count( desc ) ) continue;
      TMatrixD pred_truth = apply_ac( pred_pair.second->get_prediction() );
      slice_gen_map[ desc ] = SliceHistogram::make_slice_histogram(
        pred_truth, slice, nullptr );
    }

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string name_latex;
    std::string diff_xsec_units_denom;
    std::string diff_xsec_denom_latex;
    std::string diff_xsec_units_denom_latex;
    double other_var_width = 1.;
    for ( const auto& ov_spec : slice.other_vars_ ) {
      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb.slice_vars_.at( ov_spec.var_index_ );
      if ( high != low && std::abs(high - low) < BIG_DOUBLE ) {
        ++var_count;
        other_var_width *= ( high - low );
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;
        const std::string& temp_units = var_spec.units_;
        if ( !temp_units.empty() ) {
          diff_xsec_units_denom += " / " + temp_units;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    for ( size_t av_idx : slice.active_var_indices_ ) {
      const auto& var_spec = sb.slice_vars_.at( av_idx );
      const std::string& temp_name = var_spec.name_;
      if ( temp_name != "true bin number" ) {
        var_count += slice.active_var_indices_.size();
        name_latex += var_spec.name_;
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;

        if ( !var_spec.units_.empty() ) {
          diff_xsec_units_denom += " / " + var_spec.units_;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    // NOTE: This currently assumes that each slice is a 1D histogram
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      double width = slice_unf->hist_->GetBinWidth( b + 1 );
      width *= other_var_width;
      trans_mat( b, b ) = ( 1 / (conv_factor*width));
    }

    if (total_only) {
      const TMatrixD &A_C_temp =  *xsec.result_.add_smear_matrix_;
      double A_C_total = A_C_temp(0,0);
      trans_mat(0,0) = trans_mat(0,0) * 1/A_C_total;
      std::cout << "Ac total element: " << A_C_total << std::endl;
    }

    std::string slice_y_title;
    std::string slice_y_latex_title;
    if ( var_count > 0 && sl_idx != 4 && !total_only) {
      slice_y_title += "d";
      slice_y_latex_title += "{$d";
      if ( var_count > 1 ) {
        slice_y_title += "^{" + std::to_string( var_count ) + "}";
        slice_y_latex_title += "^{" + std::to_string( var_count ) + "}";
      }
      slice_y_title += "#sigma/" + diff_xsec_denom;
      slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;
    }
    else {
      slice_y_title += "#sigma";
      slice_y_latex_title += "\\sigma";
    }
    // Units are 10^-38 cm^2 per argon nucleus (see conversion_factor()); the
    // labels previously read 10^-39 / nucleon, from before the normalisation
    // was fixed to the BNB-note convention.
    slice_y_title += " [10^{-38} cm^{2}" + diff_xsec_units_denom + " / Ar]";
    slice_y_latex_title += "\\text{ }(10^{-38}\\text{ cm}^{2}"
      + diff_xsec_units_denom_latex + " / \\mathrm{Ar})$}";

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for ( auto& pair : slice_gen_map ) {
      auto* slice_h = pair.second;
      slice_h->transform( trans_mat );
      slice_h->hist_->GetYaxis()->SetTitle( slice_y_title.c_str() );
      slice_h->hist_->GetYaxis()->SetTitleSize( 0.055 );
      slice_h->hist_->GetYaxis()->SetTitleOffset( 0.75 );

      slice_h->hist_->GetXaxis()->SetTitleSize( 0.055 );
      slice_h->hist_->GetXaxis()->SetTitleOffset( 1 );

      slice_h->hist_->GetXaxis()->SetLabelSize( 0.05 );
      slice_h->hist_->GetYaxis()->SetLabelSize( 0.05 );
    }

    // Also transform all of the unfolded data slice histograms which have
    // specific covariance matrices
    for ( auto& sh_cov_pair : sh_cov_map ) {
      auto& slice_h = sh_cov_pair.second;
      slice_h->transform( trans_mat );
    }

    // SYSTDUMP: integrated cross section + bin-averaged fractional uncertainty per
    // systematic category (in cross-section space), for the physical slice only.
    if ( sl_idx == 0 && !sh_cov_map.empty() ) {
      TH1* href = sh_cov_map.begin()->second->hist_.get();
      int nb = href->GetNbinsX();
      double sig = 0.;
      for ( int b = 1; b <= nb; ++b ) sig += href->GetBinContent(b) * href->GetBinWidth(b);
      std::cout << "[SYSTDUMP] sigma_int " << sig << std::endl;
      for ( auto& p : sh_cov_map ) {
        TH1* hc = p.second->hist_.get();
        double s = 0.; int m = 0;
        for ( int b = 1; b <= nb; ++b ) {
          double v = hc->GetBinContent(b);
          if ( v > 0. ) { s += hc->GetBinError(b) / v; ++m; }
        }
        std::cout << "[SYSTDUMP] " << p.first << ' ' << ( m ? 100.*s/m : 0. ) << std::endl;
      }
    }

    // Dump the final differential-cross-section curves for the physical slice
    // (sl_idx 0) to a sidecar ROOT file so a clean, native multi-observable
    // closure montage can be drawn without rasterising the PDF canvases.
    if ( sl_idx == 0 ) {
      TDirectory* save_dir = gDirectory;
      std::string ch_name = OutputDirectory + "closure_hists_" + OutputFileName;
      TFile ch_file( ch_name.c_str(), "RECREATE" );
      auto write_clone = [&]( const char* key, const char* outname ) {
        auto it = slice_gen_map.find( key );
        if ( it != slice_gen_map.end() && it->second && it->second->hist_ ) {
          TH1* h = static_cast<TH1*>( it->second->hist_->Clone( outname ) );
          h->SetDirectory( &ch_file );
          h->Write();
        }
      };
      write_clone( "unfolded data", "h_unfolded_nuwro" );
      write_clone( "MicroBooNE Tune", "h_genie_tune" );
      write_clone( "truth", "h_fakedata_truth" );
      // Dump every extra file-based generator prediction (A_C-smeared, xsec
      // units) so its shape can be compared to the data outside the plot.
      for ( const auto& gp : slice_gen_map ) {
        if ( gp.first == "unfolded data" || gp.first == "MicroBooNE Tune"
          || gp.first == "truth" ) continue;
        std::string clean = gp.first;
        for ( char& ch : clean ) if ( !std::isalnum((unsigned char)ch) ) ch = '_';
        write_clone( gp.first.c_str(), ("h_gen_" + clean).c_str() );
      }
      ch_file.Close();
      if ( save_dir ) save_dir->cd();
      std::cout << "[closure-dump] wrote " << ch_name << std::endl;
    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
    std::cout << '\n';
    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram* other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if ( name == "unfolded data" ) {
        chi2_map[ name ] = SliceHistogram::Chi2Result();
        continue;
      }
      // Compare all other distributions to the unfolded data
      else {
        other = slice_gen_map.at( "unfolded data" );
      }

      // Store the chi^2 results in the map.  The unfolded covariance can be
      // numerically singular (e.g. Wiener-SVD regularization collapses it into a
      // lower-rank space, condition number ~1e13), in which case get_chi2()'s
      // matrix inversion throws "Matrix inversion failed".  This chi^2 is only a
      // diagnostic legend annotation, so fall back to a sentinel result and keep
      // producing the plots instead of aborting the whole binary.
      SliceHistogram::Chi2Result cr;
      try {
        cr = slice_h->get_chi2( *other );
      } catch ( const std::exception& e ) {
        const int nb = slice_h->hist_->GetNbinsX();
        cr = SliceHistogram::Chi2Result( std::nan(""), nb, nb, -1. );
        std::cerr << "[UnfolderNuMI] chi2 for '" << name << "' failed ("
                  << e.what() << "); using sentinel (singular covariance).\n";
      }
      const auto& chi2_result = chi2_map[ name ] = cr;

      std::cout << name << ": \u03C7\u00b2 = "
        << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
      if ( chi2_result.num_bins_ > 1 ) std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
    }

    TCanvas *c1 = new TCanvas(std::to_string(sl_idx).c_str(), std::to_string(sl_idx).c_str(), 2100, 900);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetCanvasPreferGL(0);

    bool noRatioPlot = false;   // draw a data/MC ratio panel below the main plot
    const auto rightMargin = 0.33;

    TPad* pad1 = new TPad(("pad1 slice "+std::to_string(sl_idx)).c_str(), "", 0.0, noRatioPlot ? 0.05 : 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.125);
    pad1->SetTopMargin(0.1);
    pad1->SetRightMargin(rightMargin);
    pad1->Draw();
    pad1->cd();

    slice_unf->hist_->SetLineColor( kBlack );
    slice_unf->hist_->SetLineWidth( 3 );
    slice_unf->hist_->SetMarkerStyle( kFullCircle );
    slice_unf->hist_->SetMarkerSize( 0.7 );
    slice_unf->hist_->SetStats( false );

    slice_unf->hist_->SetTitle("");

    if (total_only) {
      slice_unf->hist_->GetXaxis()->SetLabelSize(0);
      slice_unf->hist_->GetXaxis()->SetTickLength(0);
    }

    double ymax = -DBL_MAX;
    slice_unf->hist_->GetYaxis()->SetRangeUser( 0., slice_unf->hist_->GetMaximum()*1.5 );
    slice_unf->hist_->Draw( "e" );

    slice_cv->hist_->SetStats( false );
    slice_cv->hist_->SetLineColor( kAzure - 7 );
    slice_cv->hist_->SetLineWidth( 3 );
    slice_cv->hist_->SetLineStyle( 5 );
    slice_cv->hist_->Draw( "hist same" );

    if ( slice_truth ) {
      slice_truth->hist_->SetStats( false );
      slice_truth->hist_->SetLineColor( kOrange );
      slice_truth->hist_->SetLineWidth( 5 );
      slice_truth->hist_->Draw( "hist same" );
    }

    // Draw the extra file-based generator predictions (NuWro, GiBUU, NEUT, ...)
    // added to slice_gen_map above, each in a distinct colour/style.
    const int gen_colors[] = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"), TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00") }; // Okabe-Ito
    const int gen_styles[] = { 2, 7, 9, 3 };
    int gen_i = 0;
    for ( auto& pair : slice_gen_map ) {
      const std::string& nm = pair.first;
      if ( nm == "unfolded data" || nm == "MicroBooNE Tune" || nm == "truth" )
        continue;
      TH1* gh = pair.second->hist_.get();
      gh->SetStats( false );
      gh->SetLineColor( gen_colors[ gen_i % 4 ] );
      gh->SetLineStyle( gen_styles[ gen_i % 4 ] );
      gh->SetLineWidth( 3 );
      gh->Draw( "hist same" );
      ++gen_i;
    }

    // Print Values
    for (int i = 1; i <= slice_unf->hist_->GetNbinsX(); ++i) {
        double bin_content = slice_unf->hist_->GetBinContent(i);
        double bin_error = slice_unf->hist_->GetBinError(i);
        std::cout << "Bin: " << i << ", Val = " << bin_content << ", Unc = " << bin_error << std::endl;
    }

    // Draw generator predictions
    // Ac matrix
    const TMatrixD &A_C =  *xsec.result_.add_smear_matrix_;
    // Find bins for this slice
    size_t start = std::numeric_limits<size_t>::max();
    size_t stop = 0;
    for (const auto& entry : slice.bin_map_) {
        const auto& set = entry.second;
        if (set.size() != 1) {
            throw std::runtime_error("Error: set in bin_map_ has more or less than 1 entry");
        }
        start = std::min(start, *set.begin());
        stop = std::max(stop, *set.rbegin());
    }

    TMatrixD ac_hist_slice(stop - start + 1, stop - start + 1);
    for (int i = start; i <= stop; i++) {
        for (int j = start; j <= stop; j++) {
            ac_hist_slice(i - start, j - start) = A_C.operator()(i, j);
        }
    }

    TLegend *lg = new TLegend(1 - rightMargin + 0.02, 0.09, 1 - 0.02, 0.93);
    // With ~7 entries (data + up to 4 standalone generators + tune + fake-data
    // truth), most of them two-line #splitline{name}{chi2}, TLegend's auto text
    // sizing assumes few entries and makes the two lines of each entry overlap.
    // Pin an explicit small text size and a tighter symbol margin so all entries
    // fit without overlapping.
    lg->SetTextSize( 0.026 );
    lg->SetMargin( 0.16 );
    lg->SetEntrySeparation( 0.4 );

    // Add in generator predictions
    // loop through generators
    for(const auto& generator : generators) {
      int h_idx = sl_idx;
      if (total_only) h_idx = 4;
      const auto gen_hist = get_generator_hist(generator.path, h_idx, generator.scaling);
      if (gen_hist) {

        // Multiple by AC matrix
        if (!total_only) multiply_1d_hist_by_matrix(&ac_hist_slice, gen_hist);

        // Normalise by bin width
        for (int i = 1; i <= gen_hist->GetNbinsX(); ++i) {
            double bin_content = gen_hist->GetBinContent(i);
            double bin_error = gen_hist->GetBinError(i);
            double bin_width = slice_unf->hist_->GetXaxis()->GetBinWidth(i); // Get the bin width from the slice_unf histogram
            gen_hist->SetBinContent(i, bin_content / bin_width);
            gen_hist->SetBinError(i, bin_error / bin_width);
        }

        gen_hist->SetLineColor(generator.lineColor); // Set the line color
        gen_hist->SetLineWidth(generator.lineWidth); // Set the line width
        gen_hist->SetLineStyle(generator.lineStyle); // Set the line style
        gen_hist->Draw( "hist same" );

        // Calculate a Chi2 and p-value
        SliceHistogram *gen_slice_h = SliceHistogram::slice_histogram_from_histogram(*gen_hist);
        SliceHistogram::Chi2Result chi2_result;
        try {
          chi2_result = gen_slice_h->get_chi2( *slice_unf );
        } catch ( const std::exception& e ) {
          const int nb = gen_slice_h->hist_->GetNbinsX();
          chi2_result = SliceHistogram::Chi2Result( std::nan(""), nb, nb, -1. );
          std::cerr << "[UnfolderNuMI] generator chi2 for '" << generator.name
                    << "' failed (" << e.what() << "); using sentinel.\n";
        }
        //std::cout << chi2_result.chi2_ << ", p-value = " << chi2_result.p_value_ << std::endl;

        std::ostringstream oss;
        oss << "#splitline{" << generator.name << "}{"
            << "#chi^{2} = " << (chi2_result.chi2_>= 0.01 && chi2_result.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin" << (chi2_result.num_bins_ > 1 ? "s" : "") << "}";
        std::string label = oss.str();

        lg->AddEntry(gen_hist, label.c_str(), "l");
      }
    }

    for ( const auto& pair : slice_gen_map ) {
      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      const auto& chi2_result = chi2_map.at( name );

      std::string name_clean = name;

      if (name_clean == "truth") name_clean = "Fake-data truth";
      if (name_clean == "MicroBooNE Tune") name_clean = "GENIE 3.0.6 G18 #muB"; // _10a_02_11a
      //if (label == "unfolded data") label = "Unfolded Fake Data";
      if (name_clean == "unfolded data") name_clean = "Unfolded Data";

      std::ostringstream oss;

      if (name != "unfolded data") {
        oss << "#splitline{" << name_clean << "}{"
            << "#chi^{2} = " << (chi2_result.chi2_>= 0.01 && chi2_result.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin" << (chi2_result.num_bins_ > 1 ? "s" : "") << "}";
      }
      else {
        oss << name_clean;
      }
      std::string label = oss.str();

      if (name == "unfolded data") lg->AddEntry( slice_h->hist_.get(), label.c_str(), "lep" );
      else lg->AddEntry( slice_h->hist_.get(), label.c_str(), "l" );
    }

    // redraw data points to put on top
    slice_unf->hist_->Draw( "e same" );

    // Draw Legend
    lg->Draw( "same" );

     // Create the label text with the POT value
    TLatex label;
    label.SetTextAlign(12); // Set text alignment (left-aligned)
    label.SetNDC(); // Set position in normalized coordinates
    std::string labelText1( "MicroBooNE NuMI Data" );
    // POT was hardcoded to 2.2e20; use the actual data POT for this run
    // (total_pot from get_data_pot()) so the label matches the sample.
    std::string labelText2( toLatexScientific( total_pot ) + " POT" );
    label.SetTextSize(0.045);

    if (sl_idx == 0 && !total_only) {
      label.DrawLatex(0.4, 0.85, labelText1.c_str() );
      label.DrawLatex(0.4, 0.80, labelText2.c_str() );
    }
    else {
      label.DrawLatex(0.135, 0.85, labelText1.c_str() );
      label.DrawLatex(0.135, 0.80, labelText2.c_str() );
    }

    // ---- Ratio panel: data / MC for each generator ------------------------
    // Below the main cross-section plot, show data/prediction per bin for each
    // generator (GENIE tune, NuWro, ...), with the data's fractional
    // uncertainty as the error bar and a reference line at 1.
    if ( !noRatioPlot ) {
      // Suppress the main pad's x-axis labels/title (the ratio pad carries them)
      slice_unf->hist_->GetXaxis()->SetLabelSize( 0 );
      slice_unf->hist_->GetXaxis()->SetTitleSize( 0 );

      c1->cd();
      TPad* pad2 = new TPad( ("pad2_" + std::to_string(sl_idx)).c_str(), "",
        0.0, 0.0, 1.0, 0.3 );
      pad2->SetTopMargin( 0.03 );
      pad2->SetBottomMargin( 0.32 );
      pad2->SetRightMargin( rightMargin );
      pad2->SetGridy();
      pad2->Draw();
      pad2->cd();

      const TH1* h_data = slice_unf->hist_.get();
      bool first = true;
      const int rc[] = { TColor::GetColor("#0072B2"), TColor::GetColor("#009E73"), TColor::GetColor("#CC79A7"), TColor::GetColor("#D55E00"), TColor::GetColor("#E69F00") }; // Okabe-Ito
      const int rs[] = { 5, 2, 7, 9, 3 };
      int ri = 0;
      for ( const auto& gp : slice_gen_map ) {
        const std::string& nm = gp.first;
        if ( nm == "unfolded data" || nm == "truth" ) continue;
        const TH1* h_gen = gp.second->hist_.get();
        TH1* h_ratio = static_cast<TH1*>( h_data->Clone(
          ("ratio_" + nm + "_" + std::to_string(sl_idx)).c_str() ) );
        h_ratio->SetDirectory( nullptr );
        for ( int b = 1; b <= h_ratio->GetNbinsX(); ++b ) {
          double d = h_data->GetBinContent( b );
          double de = h_data->GetBinError( b );
          double g = h_gen->GetBinContent( b );
          if ( g != 0. ) {
            h_ratio->SetBinContent( b, d / g );
            h_ratio->SetBinError( b, de / g );   // data uncertainty only
          } else {
            h_ratio->SetBinContent( b, 0. );
            h_ratio->SetBinError( b, 0. );
          }
        }
        int col = ( nm == "MicroBooNE Tune" ) ? rc[0] : rc[1 + (ri % 4)];
        int sty = ( nm == "MicroBooNE Tune" ) ? rs[0] : rs[1 + (ri % 4)];
        if ( nm != "MicroBooNE Tune" ) ++ri;
        h_ratio->SetLineColor( col );
        h_ratio->SetMarkerColor( col );
        h_ratio->SetLineStyle( sty );
        h_ratio->SetLineWidth( 2 );
        h_ratio->SetMarkerStyle( 20 );
        h_ratio->SetMarkerSize( 0.7 );
        if ( first ) {
          h_ratio->SetStats( false );
          h_ratio->SetTitle( "" );
          h_ratio->GetYaxis()->SetTitle( "Data / MC" );
          h_ratio->GetYaxis()->SetRangeUser( 0.0, 3.0 );
          h_ratio->GetYaxis()->SetNdivisions( 505 );
          h_ratio->GetYaxis()->SetTitleSize( 0.12 );
          h_ratio->GetYaxis()->SetTitleOffset( 0.35 );
          h_ratio->GetYaxis()->SetLabelSize( 0.10 );
          h_ratio->GetXaxis()->SetTitleSize( 0.13 );
          h_ratio->GetXaxis()->SetTitleOffset( 1.0 );
          h_ratio->GetXaxis()->SetLabelSize( 0.11 );
          h_ratio->Draw( "e1" );
          first = false;
        } else {
          h_ratio->Draw( "e1 same" );
        }
      }
      // Reference line at 1
      double xlo = h_data->GetXaxis()->GetXmin();
      double xhi = h_data->GetXaxis()->GetXmax();
      TLine* l1 = new TLine( xlo, 1.0, xhi, 1.0 );
      l1->SetLineColor( kBlack );
      l1->SetLineStyle( 2 );
      l1->Draw( "same" );
    }

    // write to file. Name by the slice's active variable, not just the slice
    // index: the index alone (plot_slice_0.pdf ...) is reused across every
    // observable, so running several observables in sequence (e.g. via
    // run_full_chain.sh) silently overwrites all but the last. Fall back to the
    // index if the variable name is empty.
    std::string var_tag;
    if ( !slice.active_var_indices_.empty() ) {
      var_tag = sb.slice_vars_.at( slice.active_var_indices_.front() ).name_;
      // strip characters that are awkward in a filename (spaces, ROOT latex)
      std::string cleaned;
      for ( char ch : var_tag ) {
        if ( std::isalnum( static_cast<unsigned char>(ch) ) ) cleaned += ch;
      }
      var_tag = cleaned;
    }
    if ( var_tag.empty() ) var_tag = "slice";
    std::string plot_name = "unfold_output/plot_" + var_tag + "_"
      + std::to_string(sl_idx) + ".pdf";
    c1->SaveAs(plot_name.c_str());

  }

  // create plot of A_C matrix
  const Int_t n = 2;
  Double_t bins[n+1] = {0, 7, 8};
  const Char_t *labels[n] = {"E_{e}", "Total"};

  // Convert TMatrixD to TH2D
  TMatrixD temp_ac = *xsec.result_.add_smear_matrix_;
  TH2D h_A_C = TMatrixDToTH2D(temp_ac, "h_A_C", "Regularization Matrix", 0, temp_ac.GetNcols(), 0, temp_ac.GetNrows());

  gStyle->SetPalette(kBird);

  TCanvas *c_ac = new TCanvas("c_ac","A_C Matrix",200,10,1920,1080);
  c_ac->SetRightMargin(0.15);
  c_ac->SetTopMargin(0.125);
  h_A_C.SetStats(0); // Disable the statistics box
  //h_A_C.GetZaxis()->SetRangeUser(-0.5, 1.5); // Set the z range
  h_A_C.Draw("colz");
  h_A_C.GetXaxis()->SetTitle("Bin Number");
  h_A_C.GetYaxis()->SetTitle("Bin Number");
  h_A_C.GetZaxis()->SetTitle("Regularization");
  c_ac->Update();

  // Draw vertical and horizontal lines at the bin edges
  for (Int_t i = 1; i < n; i++) {
      TLine *vline = new TLine(bins[i], 0, bins[i], h_A_C.GetNbinsY());
      vline->SetLineColor(kBlack);
      vline->Draw();

      TLine *hline = new TLine(0, bins[i], h_A_C.GetNbinsX(), bins[i]);
      hline->SetLineColor(kBlack);
      hline->Draw();
  }

  for (Int_t i = 1; i <= n; i++) {
      // Draw white dotted lines from bins[i-1] to bins[i]
      TLine *vline_dotted1 = new TLine(bins[i], bins[i-1], bins[i], bins[i]);
      vline_dotted1->SetLineColor(kWhite);
      vline_dotted1->SetLineStyle(2); // Set line style to dotted
      vline_dotted1->Draw();

      TLine *hline_dotted1 = new TLine(bins[i-1], bins[i], bins[i], bins[i]);
      hline_dotted1->SetLineColor(kWhite);
      hline_dotted1->SetLineStyle(2); // Set line style to dotted
      hline_dotted1->Draw();

      if(i<n)
      {
          TLine *vline_dotted2 = new TLine(bins[i], bins[i+1], bins[i], bins[i]);
          vline_dotted2->SetLineColor(kWhite);
          vline_dotted2->SetLineStyle(2); // Set line style to dotted
          vline_dotted2->Draw();

          TLine *hline_dotted2 = new TLine(bins[i+1], bins[i], bins[i], bins[i]);
          hline_dotted2->SetLineColor(kWhite);
          hline_dotted2->SetLineStyle(2); // Set line style to dotted
          hline_dotted2->Draw();
      }
  }

  // Add labels in the middle of the intervals
  for (Int_t i = 0; i < n; i++) {
      Double_t midPoint = (bins[i] + bins[i+1]) / 2.0;
      TLatex *text = new TLatex(midPoint, 1.03*h_A_C.GetNbinsY(), labels[i]);
      text->SetTextSize(0.03); // Set text size to something smaller
      text->SetTextAlign(22); // Center alignment
      text->Draw();
      // delete text;
  }

  // add labels in each bin
  for (int i = 0; i < h_A_C.GetNbinsX(); i++) {
      for (int j = 0; j < h_A_C.GetNbinsY(); j++) {

          double bin_content = h_A_C.GetBinContent(i+1, j+1);
          if (bin_content == 0) continue;

          TLatex* latex = new TLatex(h_A_C.GetXaxis()->GetBinCenter(i+1), h_A_C.GetYaxis()->GetBinCenter(j+1), Form("%.3f",bin_content));
          latex->SetTextFont(42);
          latex->SetTextSize(0.02);
          latex->SetTextAlign(22);
          latex->Draw();
      }
  }

  h_A_C.Write();

  c_ac->SaveAs("unfold_output/plot_regularization_matrix.pdf");

  // Diagnostic plots of the intermediate analysis steps (reco spectrum,
  // background subtraction, smearceptance, efficiency), in the same style so
  // they sit side-by-side with the final unfolded slices above.
  make_analysis_step_plots( extr->get_syst(), total_pot );

}

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: Unfolder.C XSEC_Config"
      << " SLICE_Config OUTPUT_FILE\n";
    return 1;
  }

  std::string XSEC_Config( argv[1] );
  std::string SLICE_Config( argv[2] );
  std::string OutputFile( argv[3] );

  //Take the output directory from the file handed as the expected output
  //Only used for dumping to text or plot, if that option is requested in the hardcoded options at start of file
  std::string OutputDirectory = OutputFile.substr(0, OutputFile.find_last_of("/") + 1);
  std::string OutputFileName = OutputFile.substr(OutputFile.find_last_of("/") + 1);

  UnfolderNuMI(XSEC_Config, SLICE_Config, OutputDirectory, OutputFileName);

  // See the equivalent call in Unfolder.C: re-report any systematic category
  // that evaluated to zero everywhere, as the last thing printed. Warning
  // only -- the exit status stays 0 so `set -e` drivers are not aborted.
  MCC9SystematicsCalculator::report_zeroed_categories();

  return 0;
}
