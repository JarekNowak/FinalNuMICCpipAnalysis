// CCPiConfig.h
// Single source of truth for constants shared between ccpi_selection.C
// and ccpi_xsec.C.  Include this header in both macros; edit only here.

#ifndef CCPI_CONFIG_H
#define CCPI_CONFIG_H

// ─── Muon-momentum binning for unfolding ─────────────────────────────────────
// Uniform 0.1 GeV/c bins up to 2.0 GeV/c, then three wider bins extending to
// 3.0 GeV/c.  The high-momentum extension matters: ~8.6% of true signal has
// p_μ > 2.0 GeV/c and the selected-signal MCS momentum reaches ~3.05 GeV/c.
// Without it those events fall in the true-axis overflow of h_true_gen /
// h_smear_sel (absent from the response matrix) yet still appear — smeared
// down — in the in-range reco/data spectrum, biasing the unfolded upper bins.
// Adjust edges here; both macros pick up the change automatically.
//
// NOTE (2026-06-09): the pμ unfolded spectrum has negative bins in the high-pμ
// tail (1.6–3.0 GeV) and its integrated σ sits ~20% above the other observables.
// This is NOT a binning bug — it is Run1-only statistics + a ~33% MC over-
// prediction (final-cut prediction 568 vs data 427), so in the sparse tail the
// MC background alone exceeds the observed data → negative subtracted signal.
// Merging bins does not fix it (data 4 vs bkg 13 over the whole 1.6–3.0 range).
// See [[pmu-lowbin-instability]] — needs the data/MC normalisation understood
// and more data (Run2/4/5), not a binning change.

static const int    N_MU_BINS         = 22;
// First edge raised 0.10 -> 0.15 to match the signal threshold (efficiency turns
// on at ~0.15 GeV/c; the [0.10,0.15] region is near-dead at ~2%).
static const double MU_EDGES[N_MU_BINS + 1] = {
    0.15, 0.20, 0.30, 0.40, 0.50,
    0.60, 0.70, 0.80, 0.90, 1.00,
    1.10, 1.20, 1.30, 1.40, 1.50,
    1.60, 1.70, 1.80, 1.90, 2.00,
    2.30, 2.60, 3.00
};   // GeV/c

// ─── Data exposure ────────────────────────────────────────────────────────────
// Sum of beam-on POT across all runs used in the analysis [units: raw POT].
// Cross-section results scale linearly with this number.
//
//   Run 1 FHC : 3.283 × 10²⁰ POT
//   Run 2 FHC : 1.268 × 10²⁰ POT
//   Run 4 FHC : 2.075 × 10²⁰ POT
//   Run 5 FHC : 2.231 × 10²⁰ POT
//
// NOTE (unit convention fix):  TOTAL_DATA_POT is now expressed in raw POT
// (not in units of 10²⁰ POT).  Pair this with XSecConst::FLUX_PER_POT, which
// gives the integrated νμ flux in νμ/cm² per single POT.  The previous
// convention (TOTAL_DATA_POT in units of 10²⁰ × flux/10²⁰POT) had a
// placeholder flux value that was off by ~10²⁰, producing nonsense
// cross-section magnitudes.

static const double TOTAL_DATA_POT =
    3.283 * 1e20;   // = 3.283 × 10²⁰ POT (Run 1 FHC only — data file is run1_beamon)

// ─── Differential cross-section observables ──────────────────────────────────
// The analysis extracts dσ/dx for five observables.  Each carries its own bin
// edges below; the selection, unfolding, and comparison code all iterate over
// this registry, so adding/retuning an observable means editing only this file.
//
// Angle convention (decided 2026-06-06): scattering angles are cosθ w.r.t. the
// detector z-axis (cosθ = p_z/|p|), for both truth and reco.  The μ–π opening
// angle is the angle between the two tracks, in radians (0..π).
//
// Bin edges were chosen from the MC truth signal distributions (5–95% spans:
// p_π 0.14–0.92, cosθ_μ −0.32–0.98, cosθ_π −0.72–0.98, θ_μπ 22°–146°) so that
// no true bin is empty (empty bins → zero-efficiency response columns).

enum Obs { OBS_PMU = 0, OBS_PPI, OBS_CTMU, OBS_CTPI, OBS_OA, N_OBS };

// p_μ reuses MU_EDGES / N_MU_BINS defined above.

// p_π binning starts at the 0.1 GeV/c signal threshold (no empty bin → no
// zero-efficiency column) and is coarse above 0.5 where the reco estimator
// saturates (pions stop/interact, so range-based momentum loses resolution).
// Reco estimator is trk_energy_proton_v (proton-hypothesis range KE), which
// tracks true p_π much better than the muon-hypothesis range momentum.
static const int    N_PPI_BINS = 5;
// First edge raised 0.10 -> 0.175 to match the signal threshold (efficiency turns
// on at ~0.175 GeV/c; the [0.10,0.15] region is near-dead at ~1%).  Kept at 5 bins
// so the imported framework covariance still aligns by bin index.
static const double PPI_EDGES[N_PPI_BINS + 1] = {
    0.175, 0.20, 0.30, 0.40, 0.50, 1.00
};   // GeV/c

static const int    N_CTM_BINS = 12;   // cosθ_μ — forward-peaked, finer near +1
static const double CTM_EDGES[N_CTM_BINS + 1] = {
    -1.0, -0.5, 0.0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.82, 0.88, 0.93, 0.97, 1.0
};

static const int    N_CTP_BINS = 12;   // cosθ_π — broader
static const double CTP_EDGES[N_CTP_BINS + 1] = {
    -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.35, 0.5, 0.65, 0.78, 0.88, 0.95, 1.0
};

static const int    N_OA_BINS = 9;     // θ_μπ — radians, measured range [0,2.6]
// Upper bound only: θμπ < 2.6 rad (rejects the large-angle cosmic tail).  No lower
// cut — the small-angle region is high-purity (~58%), so cutting it would only cost
// signal.  Matches the signal-level OA restriction and the reco OA cut (OA_CUT).
static const double OA_EDGES[N_OA_BINS + 1] = {
    0.0, 0.314, 0.628, 0.942, 1.257, 1.571, 1.885, 2.199, 2.513, 2.6
};

static const char* OBS_KEY[N_OBS]   = { "pmu", "ppi", "costhmu", "costhpi", "thmupi" };
static const char* OBS_TITLE[N_OBS] = {
    "p_{#mu} (GeV/c)", "p_{#pi} (GeV/c)",
    "cos#theta_{#mu}", "cos#theta_{#pi}", "#theta_{#mu#pi} (rad)"
};
static const char* OBS_XSECTITLE[N_OBS] = {
    "d#sigma/dp_{#mu} [10^{-38} cm^{2}/(GeV/c)/Ar]",
    "d#sigma/dp_{#pi} [10^{-38} cm^{2}/(GeV/c)/Ar]",
    "d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/Ar]",
    "d#sigma/dcos#theta_{#pi} [10^{-38} cm^{2}/Ar]",
    "d#sigma/d#theta_{#mu#pi} [10^{-38} cm^{2}/rad/Ar]"
};
static const int OBS_NBINS[N_OBS] = { N_MU_BINS, N_PPI_BINS, N_CTM_BINS, N_CTP_BINS, N_OA_BINS };

inline const double* OBS_EDGES(int o) {
    switch (o) {
        case OBS_PMU:  return MU_EDGES;
        case OBS_PPI:  return PPI_EDGES;
        case OBS_CTMU: return CTM_EDGES;
        case OBS_CTPI: return CTP_EDGES;
        case OBS_OA:   return OA_EDGES;
    }
    return MU_EDGES;
}

#endif // CCPI_CONFIG_H
