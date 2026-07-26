// flux_breakdown.C — decompose the PPFX flux systematic to see how much is a
// genuine (correlated, normalisation) flux uncertainty vs per-bin scatter that
// the unfolding amplifies (into which the finite-MC response statistics leak).
// Reports, per observable:
//   - flux at RECO level (before unfolding): integrated (correlated) + bin-avg (diag)
//   - flux at UNFOLDED level: integrated + bin-avg  (the bin-avg is the ~20-26% quoted)
//   - MC+data stats bin-avg (unfolded) for reference
//   flux_breakdown XSEC_Config OBSNAME
#include <iostream>
#include <cmath>
#include "XSecAnalyzer/CrossSectionExtractor.hh"

static double sumall(const TMatrixD& c){ double s=0; int n=c.GetNrows();
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)s+=c(i,j); return s; }

int main( int argc, char* argv[] ) {
  if ( argc != 3 ) { std::cout << "Usage: flux_breakdown XSEC_Config OBSNAME\n"; return 1; }
  CrossSectionExtractor extr( argv[1] );
  auto xsec = extr.get_unfolded_events();
  const auto& syst = extr.get_syst();

  // ---- RECO level (pre-unfolding) ----
  auto reco_covs = syst.get_covariances();
  auto reco_sig  = syst.get_cv_ordinary_reco_signal();
  auto reco_bkg  = syst.get_cv_ordinary_reco_bkgd();
  int nr = reco_sig->GetNrows();
  auto flux_reco = reco_covs->at("flux").get_matrix();
  double reco_tot=0; for(int i=0;i<nr;i++) reco_tot += (*reco_sig)(i,0)+(*reco_bkg)(i,0);
  double flux_reco_int = std::sqrt(std::max(0.0,sumall(*flux_reco)))/reco_tot;
  double flux_reco_diag=0; int nrd=0;
  for(int i=0;i<nr;i++){ double p=(*reco_sig)(i,0)+(*reco_bkg)(i,0);
    if(p>0){ flux_reco_diag += std::sqrt(std::max(0.0,(*flux_reco)(i,i)))/p; nrd++; } }
  flux_reco_diag = nrd? flux_reco_diag/nrd : 0;

  // ---- UNFOLDED level ----
  const TMatrixD& unf = *xsec.result_.unfolded_signal_;
  int nt = unf.GetNrows();
  double unf_tot=0; for(int i=0;i<nt;i++) unf_tot += unf(i,0);
  auto& m = xsec.unfolded_cov_matrix_map_;
  const TMatrixD& flux_unf = *m.at("flux");
  double flux_unf_int  = std::sqrt(std::max(0.0,sumall(flux_unf)))/unf_tot;
  double flux_unf_diag=0; for(int i=0;i<nt;i++) if(unf(i,0)>0) flux_unf_diag += std::sqrt(std::max(0.0,flux_unf(i,i)))/unf(i,0);
  flux_unf_diag /= nt;
  // MC+data stats (unfolded, bin-avg) for reference
  double stat_diag=0;
  for(int i=0;i<nt;i++) if(unf(i,0)>0){ double v=0;
    if(m.count("SimulationStats")) v+=(*m.at("SimulationStats"))(i,i);
    if(m.count("DataStats")) v+=(*m.at("DataStats"))(i,i);
    stat_diag += std::sqrt(std::max(0.0,v))/unf(i,0); }
  stat_diag /= nt;

  printf("%-8s | flux RECO: int(norm)=%4.1f%% bin-avg=%4.1f%% | flux UNFOLDED: int(norm)=%4.1f%% bin-avg=%4.1f%% | amp(diag)=%.2fx | stat(unf)=%4.1f%%\n",
    argv[2], 100*flux_reco_int, 100*flux_reco_diag, 100*flux_unf_int, 100*flux_unf_diag,
    flux_reco_diag>0? flux_unf_diag/flux_reco_diag:0, 100*stat_diag);
  return 0;
}
