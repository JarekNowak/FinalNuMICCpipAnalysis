// WienerSVD.h
// Wiener-SVD unfolding for CCπ⁺ differential cross section extraction.
//
// Algorithm reference:
//   Tang, Zhao et al., "A new method to unfold physics distributions"
//   arXiv:1710.01242 (2017).
//
// Usage overview:
//   1. Fill h_smear_sel and h_true_gen during the main event loop
//      (see ccpi_selection.C).
//   2. Call RunWienerSVD() to get WienerSVDResult.
//   3. Call ComputeXSec() to convert unfolded counts to dσ/dp_μ.
//   4. Call PlotWienerSVDDiagnostics() for QA plots.
//
// Coordinate convention for h_smear_sel (TH2D):
//   x-axis = true  muon momentum [GeV/c]   (column index j)
//   y-axis = reco  muon momentum [GeV/c]   (row    index i)
//
//   Fill as: h_smear_sel->Fill(mc_muon_momentum, reco_pmu, weight)
//                                ^x = true          ^y = reco
//
// The response matrix is constructed as:
//   A[i][j] = h_smear_sel(j+1, i+1) / h_true_gen(j+1)
// This encodes selection efficiency × migration probability in one matrix,
// so the unfolded spectrum directly represents true interaction rates
// (no separate efficiency division required).
//
// ─── Bug-fix changelog ───────────────────────────────────────────────────────
//  FIX-1  RunWienerSVD: guard against empty h_smear_sel / h_true_gen; these
//         caused all Wiener filter coefficients to be zero (silent all-zero output).
//  FIX-2  RunWienerSVD: C_inv.Invert() now checks determinant and aborts on
//         singular matrix instead of silently producing garbage.
//  FIX-3  RunWienerSVD: added axis-swap diagnostic on the built response matrix;
//         nr == nt made the bin-count check trivially pass for transposed fills.
//  FIX-4  RunWienerSVD: zero-efficiency true bins now reported explicitly.
//  FIX-5  ComputeXSec: guard against norm_base == 0 (unset POT/flux constants).
//  FIX-6  RunWienerSVD: explicitly initialise res.success = false and pre-size
//         all result TVectorD/TMatrixD fields at function entry.  Cling
//         (ROOT's interpreter) sometimes ignored the in-class default
//         initializer `bool success = false;`, leaving the flag with stack
//         garbage.  If an early return then fired, the caller saw success=true
//         and wrote a file whose TH1D bins came from out-of-bounds TVectorD
//         reads → all NaN content and a blank dσ/dp plot.
// ─────────────────────────────────────────────────────────────────────────────

#ifndef WIENER_SVD_H
#define WIENER_SVD_H

#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

// ─────────────────────────────────────────────────────────────────────────────
// Physics constants for cross-section normalisation
// ─────────────────────────────────────────────────────────────────────────────

namespace XSecConst {

    // Liquid argon properties
    static const double rho_LAr = 1.3836;          // g/cm³
    static const double A_Ar    = 39.948;           // g/mol
    static const double N_Avo   = 6.02214076e23;    // Avogadro's number

    // ── Fiducial volume boundaries ────────────────────────────────────────────
    // These MUST match the inFV() function defined in FV_new.h.  Keep them in
    // sync: the unfolded rate is normalised per argon nucleus in this volume,
    // so any mismatch with the signal-definition FV biases the cross section.
    //
    static const double FV_x_lo  =   10.0;          // cm  (FV_new.h: FVxmin)
    static const double FV_x_hi  =  246.0;          // cm  (FV_new.h: FVxmax)
    static const double FV_y_lo  = -101.0;          // cm  (FV_new.h: FVymin)
    static const double FV_y_hi  =  101.0;          // cm  (FV_new.h: FVymax)
    static const double FV_z_lo  =   10.0;          // cm  (FV_new.h: FVzmin)
    static const double FV_z_hi  =  986.0;          // cm  (FV_new.h: FVzmax)

    // FV_new.h also removes a dead region in z that lies inside [FV_z_lo, FV_z_hi]
    // (the full slab spans the x–y face), so its volume must be subtracted.
    static const double FV_dead_z_lo = 675.1;       // cm  (FV_new.h: deadzmin)
    static const double FV_dead_z_hi = 775.1;       // cm  (FV_new.h: deadzmax)

    inline double FV_volume_cm3() {
        const double dz_total = FV_z_hi - FV_z_lo;
        const double dz_dead  = FV_dead_z_hi - FV_dead_z_lo;
        return (FV_x_hi - FV_x_lo) * (FV_y_hi - FV_y_lo) * (dz_total - dz_dead);
    }

    // Number of argon nuclei in the fiducial volume
    inline double N_targets_Ar() {
        return rho_LAr * FV_volume_cm3() * N_Avo / A_Ar;
    }

    // ── NuMI FHC integrated νμ flux ───────────────────────────────────────────
    // Units: νμ / cm² per **single POT**, integrated over all neutrino energies.
    //
    // *** IMPORTANT ***  Replace this placeholder with the value obtained by
    // integrating your PPFX flux histogram (ppfx_flux_numu.root or similar)
    // over neutrino energy:
    //   Φ = ∫ dΦ/dE_ν dE_ν   [ν_μ/cm²/POT]
    //
    // Representative MicroBooNE NuMI FHC estimate: ~1.12e-9 νμ/cm²/POT
    // (see, e.g., MicroBooNE Public Note 1087).
    //
    // NOTE (unit convention fix):  Previously this constant was called
    // FLUX_PER_1E20POT and was meant to be in νμ/cm² per 10²⁰ POT, but the
    // placeholder value 1.12e-9 is in fact a *per-POT* number.  The
    // multiplication with TOTAL_DATA_POT (which was in units of 10²⁰)
    // therefore produced a Phi_total that was 10²⁰ too small, and the
    // resulting dσ/dp came out ~10²⁰ too large.  Both constants are now in
    // raw-POT units to remove the ambiguity.
    // (νμ + ν̄μ)/cm²/POT: TH1::Integral("width") of h_numu + h_numubar from
    // uboone_numi_flux_histograms.root (200 bins, 0–20 GeV, 0.1 GeV width).
    // νμ = 1.476120e-9, ν̄μ = 8.963927e-10; combined for charge-blind |nu_pdg|==14 signal.
    static const double FLUX_PER_POT  = 2.372512e-9; // (νμ+ν̄μ)/cm²/POT, ν̄μ fraction 37.8%
    static const bool   FLUX_CONFIRMED = true;

} // namespace XSecConst

// ─────────────────────────────────────────────────────────────────────────────
// Data structure returned by RunWienerSVD
// ─────────────────────────────────────────────────────────────────────────────

struct WienerSVDResult {
    int      n_true;               // number of true bins
    int      n_reco;               // number of reco bins
    TVectorD x_unfolded;           // unfolded true spectrum (data-POT-equivalent counts)
    TVectorD x_unfolded_err;       // sqrt of diagonal of cov_stat
    TMatrixD cov_stat;             // statistical covariance of x_unfolded
    TMatrixD response;             // response matrix A (n_reco × n_true)
    TVectorD efficiency;           // selection efficiency per true bin (= colsum of A)
    TVectorD singular_values;      // SVD singular values σ_k (descending)
    TVectorD wiener_filter;        // Wiener coefficients w_k ∈ [0,1]
    TVectorD data_minus_bkg;       // background-subtracted reco data (input to unfold)
    TVectorD x_bias;               // estimated regularisation bias per true bin
    // "Additional smearing" matrix A_C = C⁻¹ V diag(w) Vᵀ C
    TMatrixD Ac;
    bool     success = false;

    WienerSVDResult() = default;
    WienerSVDResult(const WienerSVDResult&) = default;

    // ROOT's TVectorD/TMatrixD::operator= requires matching sizes, so the
    // compiler-generated copy assignment silently fails (prints "not compatible"
    // and no-ops) when assigning a sized result into a default-constructed
    // WienerSVDResult (whose members are all size 0).  This explicit operator
    // calls ResizeTo() before each assignment to handle the 0→N case.
    WienerSVDResult& operator=(const WienerSVDResult& o) {
        if (this == &o) return *this;
        n_true  = o.n_true;
        n_reco  = o.n_reco;
        success = o.success;
        x_unfolded    .ResizeTo(o.x_unfolded    .GetNrows());               x_unfolded     = o.x_unfolded;
        x_unfolded_err.ResizeTo(o.x_unfolded_err.GetNrows());               x_unfolded_err = o.x_unfolded_err;
        efficiency    .ResizeTo(o.efficiency    .GetNrows());               efficiency     = o.efficiency;
        singular_values.ResizeTo(o.singular_values.GetNrows());             singular_values = o.singular_values;
        wiener_filter .ResizeTo(o.wiener_filter .GetNrows());               wiener_filter  = o.wiener_filter;
        data_minus_bkg.ResizeTo(o.data_minus_bkg.GetNrows());               data_minus_bkg = o.data_minus_bkg;
        x_bias        .ResizeTo(o.x_bias        .GetNrows());               x_bias         = o.x_bias;
        cov_stat  .ResizeTo(o.cov_stat  .GetNrows(), o.cov_stat  .GetNcols()); cov_stat   = o.cov_stat;
        response  .ResizeTo(o.response  .GetNrows(), o.response  .GetNcols()); response   = o.response;
        Ac        .ResizeTo(o.Ac        .GetNrows(), o.Ac        .GetNcols()); Ac         = o.Ac;
        return *this;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// Build a first- or second-difference regularisation matrix of size n×n.
// "identity"  → C = I   (no extra regularisation, recommended default)
// "smooth1"   → first-difference (lower-triangular)
// "smooth2"   → second-difference (penalises curvature)
inline TMatrixD BuildCMatrix(int n, const std::string& type = "identity") {
    TMatrixD C(n, n);
    C.Zero();
    if (type == "identity") {
        for (int i = 0; i < n; i++) C[i][i] = 1.0;
    } else if (type == "smooth1") {
        C[0][0] = 1.0;
        for (int i = 1; i < n; i++) { C[i][i] = 1.0; C[i][i-1] = -1.0; }
    } else if (type == "smooth2") {
        C[0][0] = 1.0;
        C[n-1][n-1] = 1.0;
        for (int i = 1; i < n-1; i++) {
            C[i][i-1] = -1.0;
            C[i][i]   =  2.0;
            C[i][i+1] = -1.0;
        }
    } else {
        std::cerr << "[BuildCMatrix] Unknown type '" << type
                  << "'. Falling back to identity.\n";
        for (int i = 0; i < n; i++) C[i][i] = 1.0;
    }
    return C;
}

// ─────────────────────────────────────────────────────────────────────────────
// RunWienerSVD — the core unfolding
// ─────────────────────────────────────────────────────────────────────────────
//
// Parameters:
//   h_smear_sel   2D histogram (x=true p_μ, y=reco p_μ), signal MC passing all
//                 cuts, filled with Scale×ppfx×wtune weights.
//                 Fill as: h_smear_sel->Fill(mc_muon_momentum, reco_pmu, weight)
//   h_true_gen    1D true p_μ histogram for ALL signal-definition events
//                 (before any reco cuts), same weighting as h_smear_sel.
//                 This is the efficiency denominator.
//   h_data_reco   1D reco p_μ distribution from data (or MC pseudo-data).
//                 Used as the "measurement" b.
//   h_bkg_reco    1D reco p_μ prediction for all backgrounds (NC + EXT + dirt).
//                 Subtracted from h_data_reco before unfolding.
//   ctype         Regularisation matrix type: "identity"|"smooth1"|"smooth2"
//   verbose       Print diagnostics to stdout.
//
WienerSVDResult RunWienerSVD(
    TH2D*              h_smear_sel,
    TH1D*              h_true_gen,
    TH1D*              h_data_reco,
    TH1D*              h_bkg_reco,
    const std::string& ctype   = "identity",
    bool               verbose = true)
{
    WienerSVDResult res;

    // ── FIX-6: Defensive init (ROOT Cling interpreter sometimes ignores the
    //   in-class default initializer `bool success = false;`, leaving the flag
    //   with stack garbage that evaluates truthy.  An early `return res;`
    //   below would then look successful to the caller, which would write an
    //   "all NaN" output file because the result TVectorDs were never sized.
    //   Explicitly zero everything up-front so partial results stay obvious.
    res.success = false;
    res.n_true  = 0;
    res.n_reco  = 0;

    const int nt = h_true_gen->GetNbinsX();
    const int nr = h_data_reco->GetNbinsX();
    res.n_true = nt;
    res.n_reco = nr;

    // Pre-size all output vectors/matrices to nt/nr and zero them, so that even
    // if we abort early the caller's bin-by-bin write loops access valid memory
    // (zero values, not NaN from out-of-bounds reads).
    res.x_unfolded     .ResizeTo(nt);     res.x_unfolded     .Zero();
    res.x_unfolded_err .ResizeTo(nt);     res.x_unfolded_err .Zero();
    res.efficiency     .ResizeTo(nt);     res.efficiency     .Zero();
    res.x_bias         .ResizeTo(nt);     res.x_bias         .Zero();
    res.cov_stat       .ResizeTo(nt, nt); res.cov_stat       .Zero();
    res.response       .ResizeTo(nr, nt); res.response       .Zero();
    res.Ac             .ResizeTo(nt, nt); res.Ac             .Zero();
    res.data_minus_bkg .ResizeTo(nr);     res.data_minus_bkg .Zero();
    // singular_values and wiener_filter are sized to nk = min(nr,nt) later, once
    // the SVD has actually run; pre-size them here to nt as a safe placeholder.
    res.singular_values.ResizeTo(nt);     res.singular_values.Zero();
    res.wiener_filter  .ResizeTo(nt);     res.wiener_filter  .Zero();

    // ── Sanity: bin counts ────────────────────────────────────────────────────
    // Note: when nt == nr the dimension check alone cannot detect an axis swap
    // (Fill(reco, true) instead of Fill(true, reco)).  A content-level diagonal
    // dominance check is performed after building A (see FIX-3 below).
    if (h_smear_sel->GetNbinsX() != nt || h_smear_sel->GetNbinsY() != nr) {
        std::cerr << "[WienerSVD] ERROR: h_smear_sel bin count mismatch: "
                  << "GetNbinsX()=" << h_smear_sel->GetNbinsX()
                  << " (want nt=" << nt << "), "
                  << "GetNbinsY()=" << h_smear_sel->GetNbinsY()
                  << " (want nr=" << nr << ").\n"
                  << "  Check that h_smear_sel is filled as Fill(true_p, reco_p).\n";
        return res;
    }
    if (h_bkg_reco->GetNbinsX() != nr) {
        std::cerr << "[WienerSVD] ERROR: h_bkg_reco has " << h_bkg_reco->GetNbinsX()
                  << " bins but h_data_reco has " << nr << ".\n";
        return res;
    }

    // ── FIX-1: Guard against empty input histograms ───────────────────────────
    // If h_true_gen is empty, z_prior = 0 → all d_svd[k] = 0 → all w_k = 0
    // → x_hat = 0.  The old code returned success = true with a silent zero
    // result.  Now we abort with a clear error so the user can diagnose the
    // upstream problem (missing file paths, failed selection, etc.).
    {
        // Count in-range integral manually to avoid TH1::Integral() missing
        // entries in the case where MU_EDGES is set too narrowly (see ccpi_xsec.C).
        double smear_sum = 0;
        for (int j = 0; j < nt; j++)
            for (int i = 0; i < nr; i++)
                smear_sum += h_smear_sel->GetBinContent(j+1, i+1);

        // Also include the underflow/overflow cross-terms as a diagnostic.
        const double smear_full = h_smear_sel->Integral(
            0, h_smear_sel->GetNbinsX()+1,
            0, h_smear_sel->GetNbinsY()+1);

        if (smear_sum < 1.0) {
            std::cerr << "[WienerSVD] ERROR: h_smear_sel is empty "
                         "(in-range integral = " << smear_sum << ").\n";
            if (smear_full > smear_sum + 0.5)
                std::cerr << "  " << (smear_full - smear_sum)
                          << " events are in overflow/underflow — "
                          << "MU_EDGES may not cover the momentum range.\n";
            std::cerr << "  Check file paths, selection cuts, and signal definition.\n";
            return res;
        }

        const double gen_sum  = h_true_gen->Integral();
        const double gen_full = h_true_gen->Integral(0, h_true_gen->GetNbinsX()+1);
        if (gen_sum < 1.0) {
            std::cerr << "[WienerSVD] ERROR: h_true_gen is empty "
                         "(in-range integral = " << gen_sum << ").\n";
            if (gen_full > gen_sum + 0.5)
                std::cerr << "  " << (gen_full - gen_sum)
                          << " events are in overflow/underflow — "
                          << "MU_EDGES may not cover the true momentum range.\n";
            return res;
        }
    }

    if (verbose) {
        std::cout << "\n╔══════════════════════════════════════╗\n"
                  << "║      Wiener-SVD Unfolding (CCπ⁺)     ║\n"
                  << "╚══════════════════════════════════════╝\n";
        std::cout << "  True bins : " << nt << "   Reco bins : " << nr << "\n";
        std::cout << "  C-matrix  : " << ctype << "\n";
    }

    // ── Step 1: Response matrix A (nr × nt) ──────────────────────────────────
    // A[i][j] = h_smear_sel(x=j+1, y=i+1) / h_true_gen(j+1)
    //   x-axis of h_smear_sel = TRUE momentum  (must match Fill convention)
    //   y-axis of h_smear_sel = RECO momentum
    // Column j sums to the total selection efficiency ε_j.

    TMatrixD A(nr, nt);
    TVectorD eff(nt);
    A.Zero();

    for (int j = 0; j < nt; j++) {
        const double n_gen_j = h_true_gen->GetBinContent(j+1);
        if (n_gen_j < 1e-12) { eff[j] = 0.0; continue; }
        double col_sum = 0;
        for (int i = 0; i < nr; i++) {
            const double n_ij = h_smear_sel->GetBinContent(j+1, i+1);
            A[i][j]  = n_ij / n_gen_j;
            col_sum += A[i][j];
        }
        eff[j] = col_sum;
    }
    res.response   = A;
    res.efficiency = eff;

    // ── FIX-4: Report zero-efficiency true bins ───────────────────────────────
    // A zero-efficiency column produces a zero singular value, so the unfolded
    // content for that bin is zero regardless of the data.  Report explicitly
    // so the user knows those bins have no MC coverage.
    {
        int n_zero_eff = 0;
        for (int j = 0; j < nt; j++) if (eff[j] < 1e-6) n_zero_eff++;
        if (n_zero_eff > 0) {
            std::cerr << "[WienerSVD] WARNING: " << n_zero_eff << " true bin(s) have "
                         "zero selection efficiency (empty A column).\n"
                      << "  Those bins will be zero in the unfolded result regardless "
                         "of the data.\n"
                      << "  Consider widening MU_EDGES or merging low-statistics bins.\n";
        }
    }

    // ── FIX-3: Axis-swap diagnostic ───────────────────────────────────────────
    // When nt == nr the dimension check above silently passes even if
    // h_smear_sel was filled as Fill(reco, true) instead of Fill(true, reco),
    // transposing A.  A well-behaved detector response matrix should be
    // diagonally dominant; flag if the off-diagonal weight exceeds 5× diagonal.
    {
        double diag_sum = 0, off_sum = 0;
        for (int j = 0; j < nt; j++)
            for (int i = 0; i < nr; i++)
                (i == j ? diag_sum : off_sum) += A[i][j];
        if (diag_sum > 0 && off_sum > 5.0 * diag_sum) {
            std::cerr << "[WienerSVD] WARNING: response matrix is strongly "
                         "off-diagonal (diag=" << diag_sum
                      << ", off=" << off_sum << ").\n"
                      << "  Probable cause: h_smear_sel filled as "
                         "Fill(reco_p, true_p) instead of Fill(true_p, reco_p).\n";
        }
    }

    if (verbose) {
        std::cout << "\n  Efficiency per true bin (col-sum of A):\n  ";
        for (int j = 0; j < nt; j++) printf("  [%d] %.3f", j, eff[j]);
        std::cout << "\n";
    }

    // ── Step 2: Build C (additional smearing / regularisation) ───────────────
    // For C = I the algorithm reduces to standard Wiener-SVD.
    // Work in the transformed space z = C × x_true.

    TMatrixD C     = BuildCMatrix(nt, ctype);
    TMatrixD C_inv = C;

    // ── FIX-2: Check C-matrix inversion ──────────────────────────────────────
    // TMatrixD::Invert() returns the determinant but the old code discarded it,
    // silently producing C_inv = I (ROOT fall-back) for a singular C.
    // For "smooth2" with many bins the matrix can become ill-conditioned.
    {
        Double_t det = 0.0;
        C_inv.Invert(&det);
        if (std::abs(det) < 1e-12) {
            std::cerr << "[WienerSVD] ERROR: C-matrix is singular or near-singular "
                         "(det = " << det << ") for ctype = '" << ctype << "'.\n"
                      << "  Try 'identity' or reduce the number of bins.\n";
            return res;
        }
        if (verbose)
            printf("  C-matrix det             : %.6e\n", det);
    }

    // Transformed response: Ã = A × C⁻¹
    TMatrixD A_tilde(nr, nt);
    A_tilde = A * C_inv;

    // ── Step 3: Background-subtracted data b = data - background ─────────────

    TVectorD b(nr), b_err2(nr);
    for (int i = 0; i < nr; i++) {
        const double d  = h_data_reco->GetBinContent(i+1);
        const double bg = h_bkg_reco ->GetBinContent(i+1);
        b[i]     = d - bg;
        // Statistical error of b: σ²(data) + σ²(background)
        b_err2[i] = std::pow(h_data_reco->GetBinError(i+1), 2)
                  + std::pow(h_bkg_reco ->GetBinError(i+1), 2);
        if (b_err2[i] < 1e-12) b_err2[i] = std::max(1.0, std::abs(b[i]));
    }
    res.data_minus_bkg = b;

    // ── Step 4: SVD decomposition of Ã = U Σ Vᵀ ──────────────────────────────

    TDecompSVD svd(A_tilde);
    if (!svd.Decompose()) {
        std::cerr << "[WienerSVD] ERROR: SVD decomposition of A_tilde failed.\n"
                  << "  This typically means the response matrix contains NaN/Inf.\n"
                  << "  Check event weights (ppfx_cv, weightTune) for pathological values.\n";
        return res;
    }

    TMatrixD U   = svd.GetU();    // nr × nr  (left singular vectors as columns)
    TMatrixD V   = svd.GetV();    // nt × nt  (right singular vectors as columns)
    TVectorD sig = svd.GetSig();  // min(nr,nt) singular values, descending order
    const int nk = sig.GetNrows();
    res.singular_values = sig;

    if (verbose) {
        std::cout << "\n  Singular values:\n  ";
        for (int k = 0; k < nk; k++) printf("  σ%d = %.5f", k, sig[k]);
        std::cout << "\n";
    }

    // ── Step 5: Signal power spectrum in SVD basis ────────────────────────────
    // Expected signal in the transformed space: z_prior = C × x_gen
    // where x_gen = MC truth (all-generated) used as Wiener-filter prior.

    TVectorD x_prior(nt);
    for (int j = 0; j < nt; j++) x_prior[j] = h_true_gen->GetBinContent(j+1);

    // z_prior = C × x_prior
    TVectorD z_prior(nt);
    for (int j = 0; j < nt; j++) {
        double zj = 0;
        for (int l = 0; l < nt; l++) zj += C[j][l] * x_prior[l];
        z_prior[j] = zj;
    }

    // d = Vᵀ × z_prior  (signal amplitude in SVD modes)
    TVectorD d_svd(nk);
    for (int k = 0; k < nk; k++) {
        double dk = 0;
        for (int j = 0; j < nt; j++) dk += V[j][k] * z_prior[j];
        d_svd[k] = dk;
    }

    // ── Step 6: Noise variance in SVD basis ───────────────────────────────────
    // For diagonal noise covariance V_b:
    //   n_k = [Uᵀ V_b U]_kk = Σ_i  U[i][k]² × σ²_b[i]

    TVectorD noise(nk);
    for (int k = 0; k < nk; k++) {
        double nk_val = 0;
        for (int i = 0; i < nr; i++)
            nk_val += U[i][k] * U[i][k] * b_err2[i];
        noise[k] = nk_val;
    }

    // ── Step 7: Wiener filter coefficients ────────────────────────────────────
    // Signal power in SVD mode k: (σ_k × d_k)²
    // Wiener filter: w_k = signal / (signal + noise)

    TVectorD wf(nk);
    for (int k = 0; k < nk; k++) {
        if (sig[k] < 1e-10) { wf[k] = 0; continue; }
        const double signal_pwr = sig[k] * sig[k] * d_svd[k] * d_svd[k];
        const double noise_pwr  = noise[k];
        wf[k] = (signal_pwr + noise_pwr > 0)
                ? signal_pwr / (signal_pwr + noise_pwr)
                : 0.0;
    }
    res.wiener_filter = wf;

    if (verbose) {
        std::cout << "\n  Wiener filter coefficients:\n  ";
        for (int k = 0; k < nk; k++) printf("  w%d = %.4f", k, wf[k]);
        std::cout << "\n";
        // Count effectively suppressed modes
        int n_suppressed = 0;
        for (int k = 0; k < nk; k++) if (wf[k] < 0.05) n_suppressed++;
        if (n_suppressed == nk)
            std::cerr << "[WienerSVD] WARNING: all " << nk
                      << " Wiener coefficients are near zero — "
                         "x_hat will be approximately zero.\n"
                      << "  Check that h_true_gen is filled with the same events "
                         "and binning as h_smear_sel.\n";
    }

    // ── Step 8: Unfold ────────────────────────────────────────────────────────
    // ẑ = V × diag(w_k / σ_k) × Uᵀ × b    (unfolded in transformed space)
    // x̂ = C⁻¹ × ẑ                          (back to physical space)

    // g = Uᵀ × b
    TVectorD g(nk);
    for (int k = 0; k < nk; k++) {
        double gk = 0;
        for (int i = 0; i < nr; i++) gk += U[i][k] * b[i];
        g[k] = gk;
    }

    // g_filtered[k] = w_k / σ_k × g[k]
    TVectorD g_filt(nk);
    g_filt.Zero();
    for (int k = 0; k < nk; k++) {
        if (sig[k] < 1e-10) continue;
        g_filt[k] = wf[k] / sig[k] * g[k];
    }

    // ẑ = V × g_filt
    TVectorD z_hat(nt);
    z_hat.Zero();
    for (int j = 0; j < nt; j++) {
        double zh = 0;
        for (int k = 0; k < nk; k++) zh += V[j][k] * g_filt[k];
        z_hat[j] = zh;
    }

    // x̂ = C⁻¹ × ẑ
    TVectorD x_hat(nt);
    x_hat.Zero();
    for (int j = 0; j < nt; j++) {
        double xh = 0;
        for (int l = 0; l < nt; l++) xh += C_inv[j][l] * z_hat[l];
        x_hat[j] = xh;
    }

    // ── Step 9: Statistical covariance propagation ────────────────────────────
    // Unfolding matrix: M = C⁻¹ × V × diag(w/σ) × Uᵀ   (nt × nr)
    // Cov(x̂) = M × V_b × Mᵀ   (for diagonal V_b)

    TMatrixD M(nt, nr);
    for (int j = 0; j < nt; j++) {
        for (int i = 0; i < nr; i++) {
            double mji = 0;
            for (int k = 0; k < nk; k++) {
                if (sig[k] < 1e-10) continue;
                double Vz = 0;
                for (int l = 0; l < nt; l++) Vz += C_inv[j][l] * V[l][k];
                mji += Vz * (wf[k] / sig[k]) * U[i][k];
            }
            M[j][i] = mji;
        }
    }

    TMatrixD Cov(nt, nt);
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double cv = 0;
            for (int i = 0; i < nr; i++) cv += M[j1][i] * b_err2[i] * M[j2][i];
            Cov[j1][j2] = cv;
        }

    TVectorD x_err(nt);
    for (int j = 0; j < nt; j++)
        x_err[j] = (Cov[j][j] > 0) ? std::sqrt(Cov[j][j]) : 0.0;

    res.x_unfolded     = x_hat;
    res.x_unfolded_err = x_err;
    res.cov_stat       = Cov;

    // ── Step 10: Additional smearing matrix and bias estimate ─────────────────
    // A_c = C⁻¹ × V × diag(w) × Vᵀ × C
    // Bias = (A_c - I) × x_prior
    // The bias shows how much Wiener regularisation pulls the result toward the
    // prior.  It is NOT subtracted — report it as a systematic.

    TMatrixD VWVt(nt, nt);
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double v = 0;
            for (int k = 0; k < nk; k++) v += V[j1][k] * wf[k] * V[j2][k];
            VWVt[j1][j2] = v;
        }
    TMatrixD Ac = C_inv * VWVt * C;
    res.Ac = Ac;

    TVectorD bias(nt);
    bias.Zero();
    for (int j = 0; j < nt; j++) {
        double bj = -x_prior[j];   // –I × x_prior
        for (int l = 0; l < nt; l++) bj += Ac[j][l] * x_prior[l];
        bias[j] = bj;
    }
    res.x_bias = bias;

    if (verbose) {
        std::cout << "\n  Unfolded spectrum (counts, data-POT equivalent):\n";
        std::cout << "  " << std::left
                  << std::setw(8) << "Bin"
                  << std::setw(12) << "Unfolded"
                  << std::setw(12) << "StatErr"
                  << std::setw(12) << "Efficiency"
                  << std::setw(12) << "Bias\n";
        for (int j = 0; j < nt; j++)
            printf("  [%2d]   %9.2f   %9.2f   %9.4f   %9.2f\n",
                   j, x_hat[j], x_err[j], eff[j], bias[j]);
    }

    res.success = true;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// RunWienerSVD_FW — framework-consistent Wiener-SVD
// ─────────────────────────────────────────────────────────────────────────────
//
// Faithful re-implementation of xsec_analyzer's WienerSVDUnfolder
// (src/utils/WienerSVDUnfolder.cxx), which follows the BNL Wiener-SVD paper
// (Qian et al.) exactly:
//   • whiten by the data covariance via a Cholesky factor Q of Cov^{-1};
//   • regularise INSIDE the SVD by decomposing R·C^{-1} (not A itself);
//   • build the Wiener filter from Eq. (3.24) using the MC prior.
// This differs from RunWienerSVD (which decomposes A=USVᵀ directly with an
// internal Poisson noise model). Use this method to get a differential result
// consistent with the framework Unfolder/UnfolderNuMI. Pass ctype="smooth2" to
// match the framework default ("Unfold WienerSVD 1 second-deriv").
//
// Same inputs / WienerSVDResult interface as RunWienerSVD, so ComputeXSec and
// PlotWienerSVDDiagnostics work unchanged. The data covariance is the diagonal
// statistical covariance var_i = (data_err_i)^2 + (bkg_err_i)^2, matching the
// inputs used for the StandaloneUnfolding cross-check.
//
// Regularisation matrices follow WienerSVDUnfolder::set_reg_matrix:
//   "identity"            → C = I
//   "smooth1"/first-deriv → C(r,r)=-1, C(r,r+1)=+1
//   "smooth2"/second-deriv→ tridiagonal (1, diag, 1) with diag=1e-6-2 (+1 ends)
inline TMatrixD BuildCMatrix_FW(int n, const std::string& type) {
    TMatrixD C(n, n);
    C.Zero();
    if (type == "identity") {
        C.UnitMatrix();
    } else if (type == "smooth1" || type == "first-deriv") {
        for (int r = 0; r < n; r++) { C[r][r] = -1.0; if (r < n-1) C[r][r+1] = 1.0; }
    } else { // "smooth2" / "second-deriv" (framework default)
        const double inv_epsilon = 1e-6;
        for (int r = 0; r < n; r++) {
            if (r < n-1) { C[r+1][r] = 1.0; C[r][r+1] = 1.0; }
            double diag = inv_epsilon - 2.0;
            if (r == 0 || r == n-1) diag += 1.0;
            C[r][r] = diag;
        }
    }
    return C;
}

WienerSVDResult RunWienerSVD_FW(
    TH2D*              h_smear_sel,
    TH1D*              h_true_gen,
    TH1D*              h_data_reco,
    TH1D*              h_bkg_reco,
    const std::string& ctype      = "smooth2",
    bool               verbose    = true,
    const TMatrixD*    ext_covmat = nullptr)  // full reco covariance (nr×nr); if
                                              // null, fall back to stat-only diag
{
    WienerSVDResult res;
    res.success = false;
    res.n_true  = 0;
    res.n_reco  = 0;

    const int nt = h_true_gen->GetNbinsX();
    const int nr = h_data_reco->GetNbinsX();
    res.n_true = nt;
    res.n_reco = nr;

    res.x_unfolded     .ResizeTo(nt);     res.x_unfolded     .Zero();
    res.x_unfolded_err .ResizeTo(nt);     res.x_unfolded_err .Zero();
    res.efficiency     .ResizeTo(nt);     res.efficiency     .Zero();
    res.x_bias         .ResizeTo(nt);     res.x_bias         .Zero();
    res.cov_stat       .ResizeTo(nt, nt); res.cov_stat       .Zero();
    res.response       .ResizeTo(nr, nt); res.response       .Zero();
    res.Ac             .ResizeTo(nt, nt); res.Ac             .Zero();
    res.data_minus_bkg .ResizeTo(nr);     res.data_minus_bkg .Zero();
    res.singular_values.ResizeTo(nt);     res.singular_values.Zero();
    res.wiener_filter  .ResizeTo(nt);     res.wiener_filter  .Zero();

    if (h_smear_sel->GetNbinsX() != nt || h_smear_sel->GetNbinsY() != nr) {
        std::cerr << "[WienerSVD_FW] ERROR: h_smear_sel bin mismatch (X="
                  << h_smear_sel->GetNbinsX() << " want " << nt << ", Y="
                  << h_smear_sel->GetNbinsY() << " want " << nr << ").\n";
        return res;
    }
    if (h_bkg_reco->GetNbinsX() != nr) {
        std::cerr << "[WienerSVD_FW] ERROR: h_bkg_reco/h_data_reco bin mismatch.\n";
        return res;
    }
    // Wiener-SVD assumes (#reco bins) >= (#true bins).
    if (nr < nt) {
        std::cerr << "[WienerSVD_FW] ERROR: nr (" << nr << ") < nt (" << nt
                  << "); Wiener-SVD requires nr >= nt.\n";
        return res;
    }

    if (verbose) {
        std::cout << "\n╔══════════════════════════════════════╗\n"
                  << "║  Wiener-SVD Unfolding (framework-cons) ║\n"
                  << "╚══════════════════════════════════════╝\n";
        std::cout << "  True bins : " << nt << "   Reco bins : " << nr << "\n";
        std::cout << "  C-matrix  : " << ctype << " (BNL formulation)\n";
    }

    // ── Smearceptance A (nr×nt), efficiency, prior, data−bkg, covariance ─────
    TMatrixD A(nr, nt); A.Zero();
    TVectorD eff(nt);
    for (int j = 0; j < nt; j++) {
        const double gen = h_true_gen->GetBinContent(j+1);
        double col = 0.0;
        if (gen > 1e-12)
            for (int i = 0; i < nr; i++) {
                A[i][j] = h_smear_sel->GetBinContent(j+1, i+1) / gen;
                col += A[i][j];
            }
        eff[j] = col;
    }
    res.response   = A;
    res.efficiency = eff;

    TMatrixD prior(nt, 1);
    for (int j = 0; j < nt; j++) prior[j][0] = h_true_gen->GetBinContent(j+1);

    TMatrixD data_signal(nr, 1);
    TMatrixD covmat(nr, nr); covmat.Zero();
    for (int i = 0; i < nr; i++) {
        const double d  = h_data_reco->GetBinContent(i+1);
        const double b  = h_bkg_reco ->GetBinContent(i+1);
        data_signal[i][0]   = d - b;
        res.data_minus_bkg[i] = d - b;
        const double ed = h_data_reco->GetBinError(i+1);
        const double eb = h_bkg_reco ->GetBinError(i+1);
        double var = ed*ed + eb*eb;
        covmat[i][i] = (var > 0.0) ? var : 1.0;   // keep invertible
    }

    // External full covariance (framework systematic+correlated cov, already
    // scaled to this POT). Replaces the stat-only diagonal above. A stat-only
    // covariance makes the Wiener filter under-regularise the heavily-migrated
    // response and biases the unfolded integral high (~26%); the framework's
    // full covariance is the methodologically correct regularisation metric.
    if (ext_covmat) {
        if (ext_covmat->GetNrows() == nr && ext_covmat->GetNcols() == nr) {
            covmat = *ext_covmat;
            if (verbose) std::cout << "  Covariance: EXTERNAL full (nr=" << nr
                                   << ", diag0=" << covmat[0][0] << ")\n";
        } else {
            std::cerr << "[WienerSVD_FW] WARNING: ext_covmat is "
                      << ext_covmat->GetNrows() << "×" << ext_covmat->GetNcols()
                      << " but nr=" << nr << "; using stat-only covariance.\n";
        }
    } else if (verbose) {
        std::cout << "  Covariance: stat-only diagonal (no external cov supplied)\n";
    }

    // ── Whitening: Q = Cholesky factor of Cov^{-1}  (mirrors WienerSVDUnfolder) ─
    TMatrixD inv_cov(covmat);
    double det_cov = 0.0; inv_cov.Invert(&det_cov);
    if (det_cov == 0.0) { std::cerr << "[WienerSVD_FW] ERROR: singular covariance.\n"; return res; }
    TDecompChol chol(inv_cov);
    if (!chol.Decompose()) { std::cerr << "[WienerSVD_FW] ERROR: Cholesky failed.\n"; return res; }
    TMatrixD Q(chol.GetU());                 // upper-triangular factor (as framework)

    TMatrixD R = Q * A;                       // nr×nt

    // ── Regularisation matrix C and its inverse ───────────────────────────────
    TMatrixD C = BuildCMatrix_FW(nt, ctype);
    TMatrixD Cinv(C);
    double det_C = 0.0; Cinv.Invert(&det_C);
    if (det_C == 0.0) { std::cerr << "[WienerSVD_FW] ERROR: singular C matrix.\n"; return res; }

    // ── SVD of R·C^{-1} = U_C · D_C · V_C^T ──────────────────────────────────
    TMatrixD RCinv = R * Cinv;                // nr×nt
    TDecompSVD svd(RCinv);
    if (!svd.Decompose()) { std::cerr << "[WienerSVD_FW] ERROR: SVD failed.\n"; return res; }
    TMatrixD U_C = svd.GetU();                // nr×nr
    TVectorD D_C = svd.GetSig();              // length nt (descending)
    TMatrixD V_C = svd.GetV();                // nt×nt
    TMatrixD U_C_tr(TMatrixD::kTransposed, U_C);
    TMatrixD V_C_tr(TMatrixD::kTransposed, V_C);

    res.singular_values.ResizeTo(D_C.GetNrows());
    res.singular_values = D_C;

    // D_C^T as nt×nr diagonal
    TMatrixD D_C_tr(nt, nr); D_C_tr.Zero();
    for (int t = 0; t < nt; t++) D_C_tr[t][t] = D_C[t];

    // ── Wiener filter (Eq. 3.24): numer = (D_C · [V_C^T C x_prior])² ──────────
    TMatrixD numer_vec = V_C_tr * C * prior;  // nt×1
    for (int e = 0; e < nt; e++) {
        const double el = numer_vec[e][0];
        const double dC = D_C[e];
        numer_vec[e][0] = dC*dC*el*el;
    }
    TMatrixD W_C(nt, nt);      W_C.Zero();
    TVectorD wfilt(nt);        wfilt.Zero();
    for (int t = 0; t < nt; t++) {
        double numer = numer_vec[t][0];
        double denom = numer + 1.0;
        double w = (denom != 0.0) ? numer/denom : 0.0;
        W_C[t][t] = w;
        wfilt[t]  = w;
    }
    res.wiener_filter.ResizeTo(nt);
    res.wiener_filter = wfilt;

    TMatrixD W_C_tilde(W_C);
    for (int t = 0; t < nt; t++) {
        const double dC = D_C[t];
        if (dC*dC > 0.0) W_C_tilde[t][t] /= dC*dC;
        else             W_C_tilde[t][t]  = 0.0;
    }

    // ── Additional smearing matrix A_C and unfolding matrix R_tot ────────────
    TMatrixD A_C   = Cinv * V_C * W_C * V_C_tr * C;                  // nt×nt
    TMatrixD R_tot = Cinv * V_C * W_C_tilde * D_C_tr * U_C_tr * Q;   // nt×nr
    res.Ac = A_C;

    // ── Unfold and propagate covariance ───────────────────────────────────────
    TMatrixD unf = R_tot * data_signal;                             // nt×1
    TMatrixD R_tot_tr(TMatrixD::kTransposed, R_tot);
    TMatrixD cov = R_tot * covmat * R_tot_tr;                       // nt×nt

    for (int t = 0; t < nt; t++) {
        res.x_unfolded[t]     = unf[t][0];
        res.x_unfolded_err[t] = (cov[t][t] > 0.0) ? std::sqrt(cov[t][t]) : 0.0;
    }
    res.cov_stat = cov;

    // Regularisation bias = (A_C − I)·x_prior (reported, not subtracted)
    TVectorD bias(nt); bias.Zero();
    for (int j = 0; j < nt; j++) {
        double bj = -prior[j][0];
        for (int l = 0; l < nt; l++) bj += A_C[j][l] * prior[l][0];
        bias[j] = bj;
    }
    res.x_bias = bias;

    if (verbose) {
        std::cout << "  " << std::left
                  << std::setw(8) << "Bin" << std::setw(12) << "Unfolded"
                  << std::setw(12) << "StatErr" << std::setw(12) << "Eff"
                  << std::setw(12) << "Wiener_w\n";
        for (int j = 0; j < nt; j++)
            printf("  [%2d]   %9.2f   %9.2f   %9.4f   %9.4f\n",
                   j, res.x_unfolded[j], res.x_unfolded_err[j], eff[j], wfilt[j]);
    }

    res.success = true;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// RunBinByBin — efficiency-correction-only (no migration unfolding)
// ─────────────────────────────────────────────────────────────────────────────
//
// The simplest possible correction: divide background-subtracted reco counts
// by the total column-sum efficiency per bin.  Migration between true and reco
// bins is entirely ignored — reco bin i is assumed to correspond one-to-one to
// true bin i.  This method is fast and model-independent but biased when
// detector migrations are significant.  Use it as a cross-check or as the
// zeroth iteration before running IBU.
//
// Requires: identical true and reco binning (nt == nr).
// Bias estimate: (A[j][j]/eff[j] − 1) × x_true_mc[j], quantifying how much
// the diagonal fraction of the response deviates from unity.
//
WienerSVDResult RunBinByBin(
    TH2D*  h_smear_sel,
    TH1D*  h_true_gen,
    TH1D*  h_data_reco,
    TH1D*  h_bkg_reco,
    bool   verbose = true)
{
    WienerSVDResult res;
    res.success = false;
    res.n_true  = 0;
    res.n_reco  = 0;

    const int nt = h_true_gen->GetNbinsX();
    const int nr = h_data_reco->GetNbinsX();
    res.n_true = nt;
    res.n_reco = nr;

    res.x_unfolded     .ResizeTo(nt);     res.x_unfolded     .Zero();
    res.x_unfolded_err .ResizeTo(nt);     res.x_unfolded_err .Zero();
    res.efficiency     .ResizeTo(nt);     res.efficiency     .Zero();
    res.x_bias         .ResizeTo(nt);     res.x_bias         .Zero();
    res.cov_stat       .ResizeTo(nt, nt); res.cov_stat       .Zero();
    res.response       .ResizeTo(nr, nt); res.response       .Zero();
    res.Ac             .ResizeTo(nt, nt); res.Ac             .Zero();
    res.data_minus_bkg .ResizeTo(nr);     res.data_minus_bkg .Zero();
    res.singular_values.ResizeTo(nt);     res.singular_values.Zero();
    res.wiener_filter  .ResizeTo(nt);     res.wiener_filter  .Zero();

    if (nt != nr) {
        std::cerr << "[BinByBin] ERROR: true bins (" << nt
                  << ") != reco bins (" << nr
                  << "). Bin-by-bin correction requires identical true/reco binning.\n";
        return res;
    }

    // ── Response matrix ───────────────────────────────────────────────────────
    TMatrixD A(nr, nt);
    A.Zero();
    TVectorD eff(nt);
    for (int j = 0; j < nt; j++) {
        const double n_gen_j = h_true_gen->GetBinContent(j+1);
        if (n_gen_j < 1e-12) { eff[j] = 0.0; continue; }
        for (int i = 0; i < nr; i++) {
            A[i][j] = h_smear_sel->GetBinContent(j+1, i+1) / n_gen_j;
            eff[j] += A[i][j];
        }
    }
    res.response   = A;
    res.efficiency = eff;

    // ── Background-subtracted data ────────────────────────────────────────────
    TVectorD b(nr), b_err2(nr);
    for (int i = 0; i < nr; i++) {
        b[i]      = h_data_reco->GetBinContent(i+1) - h_bkg_reco->GetBinContent(i+1);
        b_err2[i] = std::pow(h_data_reco->GetBinError(i+1), 2)
                  + std::pow(h_bkg_reco ->GetBinError(i+1), 2);
        if (b_err2[i] < 1e-12) b_err2[i] = std::max(1.0, std::abs(b[i]));
    }
    res.data_minus_bkg = b;

    // ── Bin-by-bin correction ─────────────────────────────────────────────────
    // x_hat[j] = (data[j] - bkg[j]) / eff[j]
    // Ac is diagonal: Ac[j][j] = A[j][j] / eff[j]  (fraction of diagonal response)
    for (int j = 0; j < nt; j++) {
        if (eff[j] < 1e-10) continue;
        res.x_unfolded[j]     = b[j] / eff[j];
        res.cov_stat[j][j]    = b_err2[j] / (eff[j] * eff[j]);
        res.x_unfolded_err[j] = std::sqrt(res.cov_stat[j][j]);
        res.Ac[j][j]          = A[j][j] / eff[j];
        // Bias from neglected off-diagonal migration: (A[j][j]/eff[j] - 1) * x_mc[j]
        res.x_bias[j] = (res.Ac[j][j] - 1.0) * h_true_gen->GetBinContent(j+1);
    }

    if (verbose) {
        std::cout << "\n╔════════════════════════════════════════╗\n"
                  << "║   Bin-by-Bin Correction (CCπ⁺)         ║\n"
                  << "╚════════════════════════════════════════╝\n";
        std::cout << "  True/Reco bins : " << nt << "\n";
        std::cout << "\n  Corrected spectrum (counts):\n";
        std::cout << "  " << std::left
                  << std::setw(8)  << "Bin"
                  << std::setw(12) << "Corrected"
                  << std::setw(12) << "StatErr"
                  << std::setw(12) << "Efficiency"
                  << std::setw(12) << "Bias\n";
        for (int j = 0; j < nt; j++)
            printf("  [%2d]   %9.2f   %9.2f   %9.4f   %9.2f\n",
                   j, res.x_unfolded[j], res.x_unfolded_err[j], eff[j], res.x_bias[j]);
    }

    res.success = true;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// RunIBU — D'Agostini iterative Bayesian unfolding
// ─────────────────────────────────────────────────────────────────────────────
//
// Implements the algorithm from: G. D'Agostini, NIM A362 (1995) 487.
// Also known as "Bayes unfolding" (RooUnfold Bayes method).
//
// Each iteration applies the Bayes posterior:
//   w[i][j] = A[i][j] × x_prior[j] / (Σ_j' A[i][j'] × x_prior[j'])
//   x_new[j] = Σ_i b[i] × w[i][j] / ε[j]
// and updates x_prior ← x_new.
//
// The initial prior is the MC truth spectrum h_true_gen.
// Statistical covariance is propagated linearly through the final-iteration
// unfolding matrix M (standard approximation; underestimates for n_iter > 1).
//
// Parameters:
//   n_iter   Number of iterations (4 is a common starting point; scan for
//            convergence using the bias estimate).
//
WienerSVDResult RunIBU(
    TH2D*  h_smear_sel,
    TH1D*  h_true_gen,
    TH1D*  h_data_reco,
    TH1D*  h_bkg_reco,
    int    n_iter  = 4,
    bool   verbose = true)
{
    WienerSVDResult res;
    res.success = false;
    res.n_true  = 0;
    res.n_reco  = 0;

    const int nt = h_true_gen->GetNbinsX();
    const int nr = h_data_reco->GetNbinsX();
    res.n_true = nt;
    res.n_reco = nr;

    res.x_unfolded     .ResizeTo(nt);     res.x_unfolded     .Zero();
    res.x_unfolded_err .ResizeTo(nt);     res.x_unfolded_err .Zero();
    res.efficiency     .ResizeTo(nt);     res.efficiency     .Zero();
    res.x_bias         .ResizeTo(nt);     res.x_bias         .Zero();
    res.cov_stat       .ResizeTo(nt, nt); res.cov_stat       .Zero();
    res.response       .ResizeTo(nr, nt); res.response       .Zero();
    res.Ac             .ResizeTo(nt, nt); res.Ac             .Zero();
    res.data_minus_bkg .ResizeTo(nr);     res.data_minus_bkg .Zero();
    res.singular_values.ResizeTo(nt);     res.singular_values.Zero();
    res.wiener_filter  .ResizeTo(nt);     res.wiener_filter  .Zero();

    if (n_iter < 1) {
        std::cerr << "[IBU] ERROR: n_iter must be >= 1 (got " << n_iter << ").\n";
        return res;
    }

    // ── Response matrix ───────────────────────────────────────────────────────
    TMatrixD A(nr, nt);
    A.Zero();
    TVectorD eff(nt);
    for (int j = 0; j < nt; j++) {
        const double n_gen_j = h_true_gen->GetBinContent(j+1);
        if (n_gen_j < 1e-12) { eff[j] = 0.0; continue; }
        for (int i = 0; i < nr; i++) {
            A[i][j] = h_smear_sel->GetBinContent(j+1, i+1) / n_gen_j;
            eff[j] += A[i][j];
        }
    }
    res.response   = A;
    res.efficiency = eff;

    // ── Background-subtracted data ────────────────────────────────────────────
    TVectorD b(nr), b_err2(nr);
    for (int i = 0; i < nr; i++) {
        b[i]      = h_data_reco->GetBinContent(i+1) - h_bkg_reco->GetBinContent(i+1);
        b_err2[i] = std::pow(h_data_reco->GetBinError(i+1), 2)
                  + std::pow(h_bkg_reco ->GetBinError(i+1), 2);
        if (b_err2[i] < 1e-12) b_err2[i] = std::max(1.0, std::abs(b[i]));
    }
    res.data_minus_bkg = b;

    // ── D'Agostini iterations ─────────────────────────────────────────────────
    // Initial prior: MC truth spectrum (unnormalized)
    TVectorD x_prior(nt);
    for (int j = 0; j < nt; j++) x_prior[j] = h_true_gen->GetBinContent(j+1);

    if (verbose) {
        std::cout << "\n╔════════════════════════════════════════╗\n"
                  << "║   Iterative Bayesian Unfolding (IBU)   ║\n"
                  << "╚════════════════════════════════════════╝\n";
        std::cout << "  Iterations : " << n_iter
                  << "   True bins : " << nt
                  << "   Reco bins : " << nr << "\n";
    }

    // Jacobian J[j][l] = dx_N[j]/db[l], propagated through all iterations via
    // the chain rule.  The initial prior x_0 = x_MC is independent of b, so J_0 = 0.
    // At each step: J_{n+1} = M_n + T_n × J_n
    // where T_n[j][k] = Σ_i b[i] × (M_n[j][i]/x_n[j]) × (δ_{jk} − q_n[i][k])
    //       q_n[i][k] = A[i][k] × x_n[k] / expct_n[i]
    // For n_iter=1, J = M (identical to the previous linear approximation).
    TMatrixD J(nt, nr);
    J.Zero();

    TMatrixD M_unfold(nt, nr);

    for (int iter = 0; iter < n_iter; iter++) {
        M_unfold.Zero();
        TVectorD x_new(nt);
        x_new.Zero();

        // Cache per-reco-bin quantities needed for the Jacobian update.
        TVectorD expct_v(nr);
        expct_v.Zero();
        TMatrixD q_mat(nr, nt);   // q_mat[i][k] = A[i][k]*x_prior[k]/expct_i
        q_mat.Zero();

        for (int i = 0; i < nr; i++) {
            double expct = 0;
            for (int j = 0; j < nt; j++) expct += A[i][j] * x_prior[j];
            expct_v[i] = expct;
            if (expct < 1e-12) continue;
            for (int k = 0; k < nt; k++)
                q_mat[i][k] = A[i][k] * x_prior[k] / expct;

            for (int j = 0; j < nt; j++) {
                if (eff[j] < 1e-10) continue;
                const double w_ij = A[i][j] * x_prior[j] / (expct * eff[j]);
                x_new[j]         += b[i] * w_ij;
                M_unfold[j][i]    = w_ij;
            }
        }

        // Transfer matrix T[j][k] encodes d(M_n × b)[j] / d(x_prior)[k].
        // J_new = M_unfold + T × J  (C++ evaluates RHS before assigning)
        TMatrixD T(nt, nt);
        T.Zero();
        for (int j = 0; j < nt; j++) {
            if (x_prior[j] < 1e-12) continue;
            for (int k = 0; k < nt; k++) {
                double sum = 0;
                for (int i = 0; i < nr; i++) {
                    if (expct_v[i] < 1e-12) continue;
                    const double delta = (j == k) ? 1.0 : 0.0;
                    sum += b[i] * (M_unfold[j][i] / x_prior[j]) * (delta - q_mat[i][k]);
                }
                T[j][k] = sum;
            }
        }
        J = M_unfold + T * J;

        x_prior = x_new;

        if (verbose)
            printf("  Iter %2d : total unfolded = %.2f\n", iter+1, x_new.Sum());
    }
    res.x_unfolded = x_prior;

    // ── Covariance via propagated Jacobian J = dx_N/db ───────────────────────
    // Cov(x_N) = J × diag(b_err²) × Jᵀ.  For n_iter=1 this is identical to the
    // previous approximation; for n_iter>1 it accounts for the nonlinear
    // dependence of the IBU prior on the data through the full chain-rule Jacobian.
    TMatrixD Cov(nt, nt);
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double cv = 0;
            for (int i = 0; i < nr; i++)
                cv += J[j1][i] * b_err2[i] * J[j2][i];
            Cov[j1][j2] = cv;
        }
    TVectorD x_err(nt);
    for (int j = 0; j < nt; j++)
        x_err[j] = (Cov[j][j] > 0) ? std::sqrt(Cov[j][j]) : 0.0;
    res.x_unfolded_err = x_err;
    res.cov_stat       = Cov;

    // ── Resolution matrix Ac = M_unfold × A (nt × nt) ────────────────────────
    // Shows how each true bin leaks into unfolded bins; (Ac - I) × x_prior = bias.
    TMatrixD Ac(nt, nt);
    Ac.Zero();
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double v = 0;
            for (int i = 0; i < nr; i++) v += M_unfold[j1][i] * A[i][j2];
            Ac[j1][j2] = v;
        }
    res.Ac = Ac;

    // ── Bias estimate: (Ac − I) × x_mc ───────────────────────────────────────
    TVectorD x_mc(nt);
    for (int j = 0; j < nt; j++) x_mc[j] = h_true_gen->GetBinContent(j+1);
    for (int j = 0; j < nt; j++) {
        double bj = -x_mc[j];
        for (int l = 0; l < nt; l++) bj += Ac[j][l] * x_mc[l];
        res.x_bias[j] = bj;
    }

    if (verbose) {
        std::cout << "\n  Unfolded spectrum (counts):\n";
        std::cout << "  " << std::left
                  << std::setw(8)  << "Bin"
                  << std::setw(12) << "Unfolded"
                  << std::setw(12) << "StatErr"
                  << std::setw(12) << "Efficiency"
                  << std::setw(12) << "Bias\n";
        for (int j = 0; j < nt; j++)
            printf("  [%2d]   %9.2f   %9.2f   %9.4f   %9.2f\n",
                   j, res.x_unfolded[j], res.x_unfolded_err[j], eff[j], res.x_bias[j]);
    }

    res.success = true;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// RunTikhonov — Tikhonov-regularized least-squares unfolding
// ─────────────────────────────────────────────────────────────────────────────
//
// Minimises the penalised residual:
//   χ²_τ(x) = ||A x − b||² + τ ||C x||²
//
// Closed-form solution:
//   x* = (AᵀA + τ CᵀC)⁻¹ Aᵀ b
//
// Resolution matrix:
//   Ac = (AᵀA + τ CᵀC)⁻¹ AᵀA
//
// Statistical covariance:
//   Cov = M V_b Mᵀ  where  M = (AᵀA + τ CᵀC)⁻¹ Aᵀ  (nt × nr)
//
// Parameters:
//   tau    Regularisation strength (≥ 0).  τ = 0 → ordinary least-squares.
//          Larger τ → smoother solution, larger bias, smaller variance.
//          Scan τ using the bias/variance tradeoff (L-curve or GCV criterion).
//   ctype  Regularisation matrix type: "identity"|"smooth1"|"smooth2".
//          "smooth2" (second-difference penalty) is the recommended default.
//
WienerSVDResult RunTikhonov(
    TH2D*              h_smear_sel,
    TH1D*              h_true_gen,
    TH1D*              h_data_reco,
    TH1D*              h_bkg_reco,
    double             tau   = 1.0,
    const std::string& ctype = "smooth2",
    bool               verbose = true)
{
    WienerSVDResult res;
    res.success = false;
    res.n_true  = 0;
    res.n_reco  = 0;

    const int nt = h_true_gen->GetNbinsX();
    const int nr = h_data_reco->GetNbinsX();
    res.n_true = nt;
    res.n_reco = nr;

    res.x_unfolded     .ResizeTo(nt);     res.x_unfolded     .Zero();
    res.x_unfolded_err .ResizeTo(nt);     res.x_unfolded_err .Zero();
    res.efficiency     .ResizeTo(nt);     res.efficiency     .Zero();
    res.x_bias         .ResizeTo(nt);     res.x_bias         .Zero();
    res.cov_stat       .ResizeTo(nt, nt); res.cov_stat       .Zero();
    res.response       .ResizeTo(nr, nt); res.response       .Zero();
    res.Ac             .ResizeTo(nt, nt); res.Ac             .Zero();
    res.data_minus_bkg .ResizeTo(nr);     res.data_minus_bkg .Zero();
    res.singular_values.ResizeTo(nt);     res.singular_values.Zero();
    res.wiener_filter  .ResizeTo(nt);     res.wiener_filter  .Zero();

    if (tau < 0) {
        std::cerr << "[Tikhonov] ERROR: tau must be >= 0 (got " << tau << ").\n";
        return res;
    }

    // ── Response matrix ───────────────────────────────────────────────────────
    TMatrixD A(nr, nt);
    A.Zero();
    TVectorD eff(nt);
    for (int j = 0; j < nt; j++) {
        const double n_gen_j = h_true_gen->GetBinContent(j+1);
        if (n_gen_j < 1e-12) { eff[j] = 0.0; continue; }
        for (int i = 0; i < nr; i++) {
            A[i][j] = h_smear_sel->GetBinContent(j+1, i+1) / n_gen_j;
            eff[j] += A[i][j];
        }
    }
    res.response   = A;
    res.efficiency = eff;

    // ── Background-subtracted data ────────────────────────────────────────────
    TVectorD b(nr), b_err2(nr);
    for (int i = 0; i < nr; i++) {
        b[i]      = h_data_reco->GetBinContent(i+1) - h_bkg_reco->GetBinContent(i+1);
        b_err2[i] = std::pow(h_data_reco->GetBinError(i+1), 2)
                  + std::pow(h_bkg_reco ->GetBinError(i+1), 2);
        if (b_err2[i] < 1e-12) b_err2[i] = std::max(1.0, std::abs(b[i]));
    }
    res.data_minus_bkg = b;

    // ── Regularisation matrix ─────────────────────────────────────────────────
    TMatrixD C   = BuildCMatrix(nt, ctype);
    TMatrixD Ct  = TMatrixD(TMatrixD::kTransposed, C);
    TMatrixD CtC = Ct * C;

    // ── Normal equations: H = AᵀA + τ CᵀC ───────────────────────────────────
    TMatrixD AtA(nt, nt);
    AtA.Zero();
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double v = 0;
            for (int i = 0; i < nr; i++) v += A[i][j1] * A[i][j2];
            AtA[j1][j2] = v;
        }

    TMatrixD H = AtA;
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++)
            H[j1][j2] += tau * CtC[j1][j2];

    TMatrixD H_inv = H;
    Double_t det   = 0.0;
    H_inv.Invert(&det);
    if (std::abs(det) < 1e-30) {
        std::cerr << "[Tikhonov] ERROR: regularised matrix H is singular "
                     "(det = " << det << ").\n"
                  << "  Try a larger tau or switch to ctype='identity'.\n";
        return res;
    }

    // ── x* = H⁻¹ Aᵀ b ───────────────────────────────────────────────────────
    TVectorD Atb(nt);
    Atb.Zero();
    for (int j = 0; j < nt; j++) {
        double v = 0;
        for (int i = 0; i < nr; i++) v += A[i][j] * b[i];
        Atb[j] = v;
    }

    TVectorD x_hat(nt);
    x_hat.Zero();
    for (int j = 0; j < nt; j++) {
        double v = 0;
        for (int k = 0; k < nt; k++) v += H_inv[j][k] * Atb[k];
        x_hat[j] = v;
    }
    res.x_unfolded = x_hat;

    // ── Unfolding matrix M = H⁻¹ Aᵀ (nt × nr) ───────────────────────────────
    TMatrixD M(nt, nr);
    for (int j = 0; j < nt; j++)
        for (int i = 0; i < nr; i++) {
            double v = 0;
            for (int k = 0; k < nt; k++) v += H_inv[j][k] * A[i][k];
            M[j][i] = v;
        }

    // ── Statistical covariance: Cov = M V_b Mᵀ ───────────────────────────────
    TMatrixD Cov(nt, nt);
    for (int j1 = 0; j1 < nt; j1++)
        for (int j2 = 0; j2 < nt; j2++) {
            double cv = 0;
            for (int i = 0; i < nr; i++)
                cv += M[j1][i] * b_err2[i] * M[j2][i];
            Cov[j1][j2] = cv;
        }
    TVectorD x_err(nt);
    for (int j = 0; j < nt; j++)
        x_err[j] = (Cov[j][j] > 0) ? std::sqrt(Cov[j][j]) : 0.0;
    res.x_unfolded_err = x_err;
    res.cov_stat       = Cov;

    // ── Resolution matrix: Ac = H⁻¹ AᵀA ─────────────────────────────────────
    // Ac[j][l] shows what fraction of true bin l ends up in unfolded bin j.
    // Perfectly unbiased: Ac = I.  Rows of (Ac - I) quantify bias per bin.
    TMatrixD Ac = H_inv * AtA;
    res.Ac = Ac;

    // ── Bias: (Ac − I) × x_mc ────────────────────────────────────────────────
    TVectorD x_mc(nt);
    for (int j = 0; j < nt; j++) x_mc[j] = h_true_gen->GetBinContent(j+1);
    for (int j = 0; j < nt; j++) {
        double bj = -x_mc[j];
        for (int l = 0; l < nt; l++) bj += Ac[j][l] * x_mc[l];
        res.x_bias[j] = bj;
    }

    // ── Singular values of H (effective regularisation spectrum) ──────────────
    // Large eigenvalues of H indicate well-constrained modes;
    // small ones show where τ CᵀC dominates (regularised away).
    {
        TDecompSVD svd_h(H);
        if (svd_h.Decompose()) {
            TVectorD sv = svd_h.GetSig();
            const int nk = std::min(sv.GetNrows(), nt);
            res.singular_values.ResizeTo(nk);
            for (int k = 0; k < nk; k++) res.singular_values[k] = sv[k];
        }
    }

    if (verbose) {
        std::cout << "\n╔════════════════════════════════════════╗\n"
                  << "║   Tikhonov Regularisation (CCπ⁺)       ║\n"
                  << "╚════════════════════════════════════════╝\n";
        printf("  τ (regularisation strength) : %.4e\n", tau);
        printf("  C-matrix type               : %s\n",  ctype.c_str());
        printf("  det(H = AᵀA + τCᵀC)        : %.6e\n", det);
        std::cout << "\n  Regularised solution (counts):\n";
        std::cout << "  " << std::left
                  << std::setw(8)  << "Bin"
                  << std::setw(12) << "Unfolded"
                  << std::setw(12) << "StatErr"
                  << std::setw(12) << "Efficiency"
                  << std::setw(12) << "Bias\n";
        for (int j = 0; j < nt; j++)
            printf("  [%2d]   %9.2f   %9.2f   %9.4f   %9.2f\n",
                   j, x_hat[j], x_err[j], eff[j], res.x_bias[j]);
    }

    res.success = true;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeXSec — convert unfolded counts to dσ/dp_μ
// ─────────────────────────────────────────────────────────────────────────────
//
// Returns TGraphErrors with:
//   x = bin centre of true p_μ [GeV/c]
//   y = dσ/dp_μ  [10⁻³⁸ cm²/(GeV/c)/Ar nucleus]
//   ex = half bin width
//   ey = statistical uncertainty
//
// Normalisation:
//   dσ/dp[j] = x_hat[j] / (Φ_total × N_T × Δp_j)
//
// where Φ_total = flux_per_pot × total_data_pot  [νμ/cm²]
//       N_T     = number of argon nuclei in FV
//
// Unit convention: both total_data_pot and flux_per_pot are in raw-POT units
// (no implicit factor of 10²⁰).
//
TGraphErrors* ComputeXSec(
    const WienerSVDResult& res,
    TH1D*   h_bin_template,    // any histogram with the same binning as true axis
    double  total_data_pot,    // total data exposure [raw POT]
    double  flux_per_pot       = XSecConst::FLUX_PER_POT,
    double  n_targets          = -1.0,   // <0 → compute from FV constants
    bool    verbose            = true)
{
    if (!res.success) {
        std::cerr << "[ComputeXSec] ERROR: WienerSVD result flagged as failed.\n";
        return nullptr;
    }

    if (n_targets < 0) n_targets = XSecConst::N_targets_Ar();

    if (!XSecConst::FLUX_CONFIRMED) {
        std::cerr << "[ComputeXSec] WARNING: FLUX_PER_POT has not been confirmed — "
                     "it is still the placeholder value.\n"
                  << "  The unfolded SHAPE is valid, but the ABSOLUTE cross-section "
                     "scale below is NOT trustworthy.\n"
                  << "  Integrate your PPFX νμ flux histogram over Eν to get νμ/cm²/POT,\n"
                  << "  update FLUX_PER_POT in WienerSVD.h, then set FLUX_CONFIRMED = true.\n"
                  << "  Proceeding anyway (output is for shape/closure studies only).\n";
    }

    const double Phi_total = flux_per_pot * total_data_pot;  // νμ/cm²
    const double norm_base = Phi_total * n_targets;          // νμ/cm² × targets

    // ── FIX-5: Guard against zero normalisation ───────────────────────────────
    // If TOTAL_DATA_POT or FLUX_PER_POT was left as zero in CCPiConfig.h,
    // norm_base = 0 and every cross-section point becomes Inf/NaN.
    if (norm_base < 1e-30) {
        std::cerr << "[ComputeXSec] ERROR: normalisation denominator is zero or "
                     "near-zero (norm_base = " << norm_base << ").\n"
                  << "  total_data_pot = " << total_data_pot << "  [raw POT]\n"
                  << "  flux_per_pot   = " << flux_per_pot   << "  [νμ/cm²/POT]\n"
                  << "  n_targets      = " << n_targets      << "\n"
                  << "  Set TOTAL_DATA_POT and FLUX_PER_POT in CCPiConfig.h.\n";
        return nullptr;
    }

    if (verbose) {
        printf("\n  ── Cross-section normalisation ──\n");
        printf("  Total data POT           : %.4e POT (= %.4f × 10²⁰)\n",
               total_data_pot, total_data_pot / 1e20);
        printf("  Flux per POT             : %.4e νμ/cm²/POT\n",   flux_per_pot);
        printf("  Total integrated flux    : %.4e νμ/cm²\n",       Phi_total);
        printf("  FV volume                : %.4e cm³\n",          XSecConst::FV_volume_cm3());
        printf("  N_targets (Ar nuclei)    : %.4e\n",              n_targets);
        printf("  Normalisation factor     : %.4e\n",              norm_base);
    }

    const int nt = res.n_true;
    std::vector<double> xc(nt), xe(nt), yc(nt), ye(nt);

    for (int j = 0; j < nt; j++) {
        const double lo = h_bin_template->GetBinLowEdge(j+1);
        const double hi = h_bin_template->GetBinLowEdge(j+2);
        const double dp = hi - lo;

        xc[j] = 0.5 * (lo + hi);
        xe[j] = 0.5 * dp;

        // dσ/dp [cm²/(GeV/c)/nucleon] → multiply by 10³⁸ to get [10⁻³⁸ cm²/(GeV/c)]
        const double factor = 1e38 / (norm_base * dp);
        yc[j] = res.x_unfolded[j]     * factor;
        ye[j] = res.x_unfolded_err[j] * factor;
    }

    TGraphErrors* gr = new TGraphErrors(nt,
        xc.data(), yc.data(), xe.data(), ye.data());
    gr->SetName("xsec_dsdpmu");
    gr->SetTitle(";p_{#mu} (GeV/c);"
                 "d#sigma/dp_{#mu} [10^{-38} cm^{2}/(GeV/c)/Ar]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetLineWidth(2);
    gr->SetLineColor(kBlack);
    gr->SetMarkerColor(kBlack);

    if (verbose) {
        printf("\n  ── dσ/dp_μ [10⁻³⁸ cm²/(GeV/c)/Ar] ──\n");
        printf("  %-20s  %10s  %10s  %10s\n",
               "p_μ range (GeV/c)", "dσ/dp_μ", "stat err", "bias");
        for (int j = 0; j < nt; j++) {
            const double bias_xs = res.x_bias[j]
                                 * 1e38 / (norm_base * (h_bin_template->GetBinLowEdge(j+2)
                                                       - h_bin_template->GetBinLowEdge(j+1)));
            printf("  [%.2f, %.2f]          %10.4f  %10.4f  %10.4f\n",
                   h_bin_template->GetBinLowEdge(j+1),
                   h_bin_template->GetBinLowEdge(j+2),
                   yc[j], ye[j], bias_xs);
        }
    }

    return gr;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeIntegratedXSec — total cross section by integrating the differential
// ─────────────────────────────────────────────────────────────────────────────
//
// σ = (Σ_j x_unfolded[j]) / (Φ_total × N_T), i.e. ∫(dσ/dx) dx.  This is the bin-
// width-independent integral of the differential result, so each observable's
// differential should yield a consistent σ — a built-in closure cross-check.
// The error uses the FULL covariance: σ_err² = (1e38/norm)² × Σ_ij cov[i][j].
//
struct IntegratedXSec { double value; double error; bool ok; };

inline IntegratedXSec ComputeIntegratedXSec(
    const WienerSVDResult& res,
    double total_data_pot,
    double flux_per_pot = XSecConst::FLUX_PER_POT,
    double n_targets    = -1.0,
    bool   verbose      = true)
{
    IntegratedXSec out{0.0, 0.0, false};
    if (!res.success) return out;
    if (n_targets < 0) n_targets = XSecConst::N_targets_Ar();

    const double norm_base = flux_per_pot * total_data_pot * n_targets;
    if (norm_base < 1e-30) {
        std::cerr << "[ComputeIntegratedXSec] ERROR: zero normalisation.\n";
        return out;
    }
    const double f = 1e38 / norm_base;   // → 10⁻³⁸ cm²/Ar

    const int nt = res.n_true;
    double tot = 0.0;
    for (int j = 0; j < nt; j++) tot += res.x_unfolded[j];
    double var = 0.0;
    for (int i = 0; i < nt; i++)
        for (int j = 0; j < nt; j++) var += res.cov_stat[i][j];

    out.value = tot * f;
    out.error = (var > 0 ? std::sqrt(var) : 0.0) * f;
    out.ok    = true;

    if (verbose)
        printf("\n  ── Integrated cross section ──\n"
               "  σ = %.4f ± %.4f  [10⁻³⁸ cm²/Ar]  (stat. only)\n",
               out.value, out.error);
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// PlotWienerSVDDiagnostics — QA canvases written to outdir
// ─────────────────────────────────────────────────────────────────────────────
//
// Produces four canvases:
//   1. response_matrix.png  — colour map of A[i][j]
//   2. svd_diagnostics.png  — singular values and Wiener filter coefficients
//   3. unfolding_check.png  — reco data vs prediction, unfolded vs MC truth
//   4. xsec_dsdpmu.png      — final differential cross section with stat errors
//
void PlotWienerSVDDiagnostics(
    const WienerSVDResult& res,
    TH1D*         h_data_reco,
    TH1D*         h_bkg_reco,
    TH1D*         h_mc_sig_reco,   // signal MC in reco space (for closure check)
    TH1D*         h_true_gen,      // MC truth in true space
    TH1D*         h_bin_template,
    TGraphErrors* gr_xsec,
    const TString& outdir,
    const TString& method_label = "Wiener-SVD")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    const int nt = res.n_true;
    const int nr = res.n_reco;

    // ── Canvas 1: response matrix ─────────────────────────────────────────────
    {
        TCanvas c("c_resp", "Response Matrix", 700, 600);
        c.SetRightMargin(0.15);

        TH2D h2("h_resp", "Response matrix A;True p_{#mu} (GeV/c);Reco p_{#mu} (GeV/c)",
                nt, h_bin_template->GetXaxis()->GetXmin(),
                    h_bin_template->GetXaxis()->GetXmax(),
                nr, h_bin_template->GetXaxis()->GetXmin(),
                    h_bin_template->GetXaxis()->GetXmax());
        for (int j = 0; j < nt; j++)
            for (int i = 0; i < nr; i++)
                h2.SetBinContent(j+1, i+1, res.response[i][j]);
        h2.SetMinimum(0);
        h2.Draw("COLZ");

        TLatex lat;
        lat.SetNDC(); lat.SetTextSize(0.04);
        lat.DrawLatex(0.15, 0.92, "MicroBooNE Simulation");
        c.Print(outdir + "response_matrix.png");
    }

    // ── Canvas 2: SVD / mode diagnostics ─────────────────────────────────────
    {
        TCanvas c("c_svd", "Mode Diagnostics", 1000, 450);
        c.Divide(2, 1);

        c.cd(1)->SetLogy();
        const int nk = res.singular_values.GetNrows();
        TString sv_title = "Singular values (" + method_label + ");Mode k;#sigma_{k}";
        TH1D h_sig("h_svd_sig", sv_title, nk, -0.5, nk-0.5);
        for (int k = 0; k < nk; k++) {
            const double v = res.singular_values[k];
            if (v > 0) h_sig.SetBinContent(k+1, v);
        }
        h_sig.SetFillColor(kAzure-9); h_sig.Draw("HIST");

        c.cd(2);
        const bool is_wiener = (res.wiener_filter.Max() > 1e-6);
        TString wf_title = is_wiener
            ? TString("Wiener filter;Mode k;w_{k}")
            : TString("Resolution matrix diagonal (" + method_label + ");True bin j;Ac[j][j]");
        const int nd = is_wiener ? nk : nt;
        TH1D h_wf("h_svd_wf", wf_title, nd, -0.5, nd-0.5);
        if (is_wiener) {
            for (int k = 0; k < nk; k++) h_wf.SetBinContent(k+1, res.wiener_filter[k]);
        } else {
            // For non-Wiener methods show the diagonal of the resolution matrix Ac
            for (int j = 0; j < nt; j++) h_wf.SetBinContent(j+1, res.Ac[j][j]);
        }
        h_wf.GetYaxis()->SetRangeUser(0, 1.3);
        h_wf.SetFillColor(kOrange-3); h_wf.Draw("HIST");

        c.Print(outdir + "svd_diagnostics.png");
    }

    // ── Canvas 3: reco closure + unfolded vs truth ────────────────────────────
    {
        TCanvas c("c_check", "Unfolding Check", 1100, 500);
        c.Divide(2, 1);

        // Left: reco data vs background vs prediction
        c.cd(1);
        TH1D* h_pred = (TH1D*)h_mc_sig_reco->Clone("h_pred");
        h_pred->Add(h_bkg_reco);
        h_pred->SetLineColor(kRed+1); h_pred->SetLineWidth(2);
        h_data_reco->SetMarkerStyle(20); h_data_reco->SetMarkerSize(0.9);
        h_pred->GetXaxis()->SetTitle("Reco p_{#mu} (GeV/c)");
        h_pred->GetYaxis()->SetTitle("Events (data POT equivalent)");
        h_pred->Draw("HIST");
        h_data_reco->Draw("EP SAME");
        h_bkg_reco->SetFillColor(kGray); h_bkg_reco->SetLineColor(kGray+1);
        h_bkg_reco->Draw("HIST SAME");

        TLegend lg1(0.55, 0.65, 0.88, 0.88);
        lg1.SetBorderSize(0); lg1.SetFillStyle(0);
        lg1.AddEntry(h_data_reco, "Data (pseudo)", "lep");
        lg1.AddEntry(h_pred,      "Signal + BG MC", "l");
        lg1.AddEntry(h_bkg_reco,  "Background",    "f");
        lg1.Draw();

        // Right: unfolded vs true MC
        c.cd(2);
        const double xlo = h_bin_template->GetXaxis()->GetXmin();
        const double xhi = h_bin_template->GetXaxis()->GetXmax();
        TH1D h_unf("h_unf_disp", ";True p_{#mu} (GeV/c);Events (data POT equiv.)",
                   nt, xlo, xhi);
        for (int j = 0; j < nt; j++) {
            h_unf.SetBinContent(j+1, res.x_unfolded[j]);
            h_unf.SetBinError  (j+1, res.x_unfolded_err[j]);
        }
        h_unf.SetMarkerStyle(20); h_unf.SetMarkerSize(0.9);
        h_true_gen->SetLineColor(kRed+1); h_true_gen->SetLineWidth(2);
        double ymax = std::max(h_unf.GetMaximum(), h_true_gen->GetMaximum()) * 1.3;
        h_true_gen->GetYaxis()->SetRangeUser(0, ymax);
        h_true_gen->Draw("HIST");
        h_unf.Draw("EP SAME");

        TLegend lg2(0.55, 0.75, 0.88, 0.88);
        lg2.SetBorderSize(0); lg2.SetFillStyle(0);
        lg2.AddEntry(&h_unf,    "Unfolded",       "lep");
        lg2.AddEntry(h_true_gen,"MC truth (gen.)", "l");
        lg2.Draw();

        c.Print(outdir + "unfolding_check.png");
        delete h_pred;
    }

    // ── Canvas 4: final cross section ─────────────────────────────────────────
    {
        TCanvas c("c_xsec", "Cross Section", 700, 600);
        c.SetLeftMargin(0.15);

        if (gr_xsec) {
            gr_xsec->GetXaxis()->SetLimits(
                h_bin_template->GetXaxis()->GetXmin(),
                h_bin_template->GetXaxis()->GetXmax());
            gr_xsec->Draw("AP");

            TLine zero(h_bin_template->GetXaxis()->GetXmin(), 0,
                       h_bin_template->GetXaxis()->GetXmax(), 0);
            zero.SetLineStyle(2); zero.SetLineColor(kGray+2);
            zero.Draw("SAME");
        }

        TLatex lat;
        lat.SetNDC(); lat.SetTextSize(0.038);
        lat.DrawLatex(0.18, 0.92, "MicroBooNE, NuMI FHC");
        lat.DrawLatex(0.18, 0.86,
            (TString("#nu_{#mu} CC#pi^{+}, ") + method_label).Data());

        c.Print(outdir + "xsec_dsdpmu.png");
    }
}

#endif // WIENER_SVD_H
