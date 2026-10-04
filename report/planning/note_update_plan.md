# Analysis-note update plan — new information from the 2026-08 session

Status legend: **[READY]** content is in hand · **[PENDING]** waits on the combined build finishing.

## A. Detector-variation systematic — Run1-FHC stand-in → Run4 native  **[READY]**
The whole analysis previously used ONE Run1-FHC detVar set applied globally. We now
have **native Run4 FHC and Run4 RHC** detVar (all 9 systcalc knobs each), used as the
detVar for all runs of each mode.
- **§Systematic breakdown (tab:systbreak, L2290)**: "Detector" row (9.6–14.6%) was
  Run1-FHC. Recompute from the new FHC build (numbers changed; pull from
  `xsec_FHC5_*` covariances / syst_breakdown plots).
- **§RHC (L2549–2551)**: delete "Native RHC … detector-variation samples are not yet
  available … FHC detector variations as a stand-in". Replace with: native Run4 RHC
  detVar (10 knobs, uniform ev/POT 3.80e-16, verified νμ overlay).
- **tab:rhc caption (L2584)**: remove "FHC detector variations used as a stand-in".
- **NEW methodological note**: Run5 FHC detVar Recomb2 & SCE were intrinsic-νe
  overlays (0% νμ, inflated POT) — unusable; quarantined. Run4 is therefore used as
  the single consistent detVar set per mode (detVar is detector-physics, ~run-indep).

## B. Combined flux systematic — RETRACT "broken"  **[READY]**
- **§RHC (L2552–2566)**: the paragraph "the joint flux systematic is not yet correct
  … requires a horn-mode-aware flux-universe combination" is **WRONG** and must be
  rewritten. Finding: all 14 numuMC files carry the SAME 600 PPFX universes;
  `UniverseMaker` chains them and accumulates universe *u* index-matched across both
  modes → per-universe summing IS the rigorously correct correlated treatment (the
  cross-horn hadron-production correlation is captured automatically). No fix needed.
- **§Known Issues (tab:issues, L2011)**: add a Resolved row for the combined flux
  systematic (was mis-flagged as broken).

## C. RHC results table (tab:rhc, L2586) — new numbers  **[READY]**
Rebuilt with native Run4 RHC detVar, full exposure 1.108e21:
| obs | flux(unf) | amp | closure χ²/ndf (p) |
|---|---|---|---|
| pμ | 27.3% | 1.78× | 1.37/7 (0.99) |
| pπ | 23.6% | 1.49× | 0.21/5 (0.999) |
| cosθμ | 25.9% | 1.76× | 0.43/5 (0.99) |
| cosθπ | 28.3% | 1.99× | 1.23/5 (0.94) |
| θμπ | 19.5% | 1.29× | 0.72/7 (0.998) |

## D. FHC full-exposure results (Run4 detVar) — new/updated table  **[READY]**
FHC {Run1,2,4,5}, exposure 8.857e20, native Run4 FHC detVar:
| obs | flux(unf) | closure χ²/ndf (p) |
|---|---|---|
| pμ | 24.0% | 1.73/7 (0.97) |
| pπ | 23.9% | 0.75/5 (0.98) |
| cosθμ | 24.0% | 0.54/5 (0.99) |
| cosθπ | 27.1% | 1.75/5 (0.88) |
| θμπ | 24.6% | 1.15/7 (0.99) |
(slots into §Multi-run overlay L2392 or a new FHC-full table.)

## E. Combined FHC+RHC measurement — NEW subsection  **[PENDING build]**
- **Combined generator predictions**: `combine_comb_fte.C`, fluence-weighted
  0.458 FHC / 0.542 RHC (Φ_tot×POT per mode); 4 generators. Reconcile with the
  existing "Φ_tot=6.663e-10 (42.7% ν̄μ)" line (L2549) — that was flux-only; the
  measurement weight is fluence×POT.
- **Per-mode normalization (KEY METHOD)**: combined covariance is absolute, so each
  mode is scaled to its OWN data exposure (D_FHC=8.857e20, D_RHC=1.108e21). Per-mode
  summed_pot (FHC 2.158e22, RHC 2.782e22) + per-mode Poisson throw (PS_FHC=0.0924,
  PS_RHC=0.0717 via TChain tree index); detVar native×(2.251 FHC / 1.799 RHC). A
  single summed_pot would mis-weight the modes (MC-POT 0.62:1 vs data 0.80:1).
  Verified: throw kept FHC=271750 (= standalone FHC) + RHC=289150.
- **Combined closure + flux table**: fill from the running combined build.

## F. EXT (cosmic) handling — mode-independent  **[READY]**
- **§Data exposure (L672) / §EXT subtraction (L1937)**: EXT is cosmic → NO FHC/RHC
  distinction; matched per-run by gate-count scaling. Per-run dedicated EXT processed
  (Run1,3,4,5); **no Run2 EXT** (Run1/Run3 stand-in scaled to Run2 gates). Note the
  old FP EXT trigger (3.82M) was an EVENT count, not gates; correct gate counts are
  ~10× larger (Run1 FHC 4.58M, Run3 RHC 32.88M, …).

## G. Known-Issues table (tab:issues, L2011) — status flips  **[READY]**
- Combined joint flux systematic → **Resolved** (per-universe summing is correct).
- Native RHC detVar "not available" → **Resolved** (Run4 RHC native).
- Add: Run5 detVar Recomb2/SCE νe-sample trap → **Resolved** (Run4 stand-in).
- "Runs 2,4,5 beam-on not available" → still open (blind); note Run3 RHC + full RHC
  exposure now processed on the MC/fake-data side.

## New figures to regenerate
- Systematic-breakdown plots (fw_systbreak_*) with Run4 detVar.
- Combined closure + differential plots (5 observables) once the build finishes.
- Optionally FHC/RHC/combined overlaid dσ/dx with the 4 generators.

# ============================================================
# CORRECTNESS AUDIT (statements that are wrong / inconsistent)
# ============================================================

## AUD-1. Fake-data sample described THREE inconsistent ways
The note never settles what "fake data" the main result uses:
- "NuWro fake data" x13 (reco spectra L1476, EXT subsec L1937-1959, sel-perf L1133, ...)
- "GENIE detVar-CV fake data" x3 (integrated xsec L1653, L1664; L1909)
- "Poisson-thrown fake data" x8 (multi-run/RHC/M_A^RES L1406, L2414+)
The CURRENT multi-run FHC/RHC/combined results use **Poisson-thrown CV fake data**
(throw_cv_*). Standardize the whole note to that; keep NuWro only inside the
explicitly-labelled NuWro cross-validation (Sec nuwroclosure).

## AUD-2. fig:reco_spectra caption (L1476-1479) — WRONG on two counts  [user-flagged]
"NuWro fake data (black points)... NuWro carries roughly twice the GENIE-predicted
signal, consistent with the generator spread in CCpi+ production."
  (a) We do NOT use NuWro fake data for the main spectra -> regenerate figs
      (fw_reco_*) with the current Poisson-CV fake data and relabel.
  (b) The ~2x is attributed to "generator spread" — CONTRADICTS Sec L2225 which
      calls the same ~2x a "pure POT/flux bookkeeping factor" (1.29x detVar-CV POT
      offset x ~1.3x A_C smearing). Memory data-2x-above-generators-explained agrees
      it is the POT offset, NOT physics. Fix the caption.

## AUD-3. The "2x excess" narrative is self-contradictory
- Sec nuwroclosure (L2123-2131): "flat ~2x offset is a GENUINE NuWro-vs-GENIE-tune
  model difference."
- Sec datamc/decomposition (L2225-2227): "~1.29x ... a PURE POT/flux bookkeeping
  factor ... NOT a physics/generator difference."
These are two different fake-data samples (NuWro closure vs GENIE detVar-CV main),
but the note reads as one story and contradicts itself. Reconcile: NuWro-closure 2x
= real generator diff of NuWro; main-result apparent 2x = POT bookkeeping. State
both distinctly, or drop the superseded framing.

## AUD-4. EXT-subtraction subsection (L1937-1959) framed around "NuWro fake data"
Mechanism (add EXT to pure-MC fake data, cancels in closure) is generic. Relabel to
the current Poisson-CV fake data; the ~3.7% EXT figure should be re-derived if it
was NuWro-specific.

## AUD-5. Selection-perf / phase-space-iteration table (L1133) — "NuWro fake data"
Legacy binning study. Verify whether to regenerate with current fake data or mark
explicitly as a historical binning scan (result — the 2.6 rad cut — still stands).

## AUD-6 (already in plan A/B). detVar "FHC stand-in" (L2549-2551,2584) and the
## "combined joint flux systematic not correct" paragraph (L2552-2566) are both
## now FALSE — native Run4 detVar exists; combined flux systematic is correct.

## Cross-check still-TRUE anchors (do NOT change)
- Flux Phi=6.81159e-10 cm^-2 POT^-1 (FHC) — still correct.
- Wiener-SVD, optimised binning, phase-space cuts — unchanged.
- NuWro cross-pipeline closure as a VALIDATION exercise — legitimate, keep labelled.

# ============================================================
# ALL TABLES & FIGURES FOR ALL THREE CONFIGS (FHC / RHC / combined)
# ============================================================

## DONE (results, 3 configs each)
- tab:fhc / tab:rhc / tab:comb (flux, amp, closure)
- tab:sigint_all (integrated sigma, 3 cols)
- tab:systbreak_fhc / _rhc / _comb (systematic breakdown)
- fig:dsigma_current / _rhc / _comb (differential xsec)
- fig:reco_spectra / _rhc / _comb (reco spectra)

## CONFIG-INDEPENDENT (one version is correct for all)
- tab:pot (updated, full FHC+RHC), tab:samples (updated)
- tab:mupid / tab:pipid (selection thresholds), fig:evdisp

## TO DO -- selection diagnostics, need selection-stored histos + reprocess
All currently FHC Run-1 only. Require the extended selection + a reprocess of
numuMC + EXT + dirt for each config, then POT-weighted aggregation.
  1. fig:cutflow (Fig 1) + tab:cutflow (Table 7)  -- cut-flow yields
       -> h_cutflow_tot / h_cutflow_sig counters ADDED (commit 561d5aa);
          macros/cutflow_yields.C written.
  2. fig:nm1_key + fig:nm1_tracks -- N-1 distributions (6 vars:
       TopologicalScore, MuonPID, PionPID, MuPiOpeningangle, MuonTrackLength,
       PionTrackLength), stacked by category. -> instrument selection.
  3. fig:final_mu / _pi / _topo -- final-cut data/MC distributions per category.
  4. tab:bkgcomp -- background decomposition at the final cut.
  5. tab:phasespace -- efficiency/purity in the measured phase space.
  6. fig:pid_proton, tab:pithresh(_eff) -- verify / regenerate per config.

## EXECUTION (one coordinated pass)
When the FHC binning re-run finishes (reprocessing unblocks):
  a. Finish instrumenting CC1mu1piXp with the N-1 / final-cut / bkg-category /
     efficiency histograms (per EventCategory) alongside the cut-flow counters.
  b. Rebuild; reprocess ALL file types (numuMC + EXT + dirt) with the new binning
     AND the full instrumentation; extend binning to RHC + combined.
  c. Aggregate per config (POT-weighted) -> all diagnostic tables/figures x3.
BLIND nuance: data == pred by construction (fake data), so the "Data" column /
points equal the MC prediction; present as the MC-prediction breakdown until
unblinding.
