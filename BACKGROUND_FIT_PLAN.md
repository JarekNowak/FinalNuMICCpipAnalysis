# Joint FHC+RHC background constraint from the control regions: plan

> **Status 2026-09-16 (implemented; see report/cr_data_note.tex §5–6).** The control regions had already been opened
> on beam-on data (2026-09-09), so everything below the governance line is post-inspection work. Code:
> `xsec_analyzer/macros/bkgfit_templates.C`, `xsec_analyzer/bkgfit/` (model, fit, gain, detvar, studies, figures).
> Findings that change this plan: (1) the sideband MC `processed/sb/` predated the beam-frame fix — the CR muon-angle
> shape failure was an artefact; templates now use `processed/sb_pi0/`; (2) the sideband detector variations were not
> POT-normalised (detector term ~10x inflated); corrected, with the statistical part removed by an event-matched
> bootstrap; the frozen global test then fails with FHC Run 1 (69.5/8) and passes without it; (3) nu/nubar is not
> separable with the two beams (CC0pi nubar share 11.8% FHC vs 15.4% RHC); (4) no normalisation correction is
> warranted (CC0pi 0.985 +- 0.191 without FHC Run 1), FHC Run-1 exposure factor 0.520 +- 0.020; (5) a conditional
> constraint would take the one-bin total from 35.8% to 11.3% (Asimov), entirely through the flux and detector
> correlations, while the frozen propagation rule gains nothing; (6) the remaining CR discrepancy is in p_mu.

*2026-09-16. Inclusive CC1μ1πXp only. Blind analysis: every step up to the go/no-go gate uses
Asimov and Poisson pseudo-data. No beam-on control-region data are opened by this plan.*

## 1. Starting point and hard constraints

**What exists already**
- Four control regions, disjoint from the signal region: CC0π, π⁰, multi-π and cosmic.
  They are carried in `processed/sb/` and the beam-on skim `XSEC_CR_SKIM=1`.
- `macros/sb_protocol.C` builds the full pre-fit covariance of the 8 (region, beam) totals:
  - flux, cross-section and re-interaction multisims, index-matched across FHC and RHC;
  - 8 detector knobs as correlated unisims;
  - MC, EXT, dirt and data statistics.
- Transfer factors (`sb_transfer.C`):

  | Class | Region | T | Target purity | T uncertainty |
  |---|---|---|---|---|
  | CC0π | CC0π | ≈ 0.03 | 72% | 31–34%, detector-dominated |
  | CC+NC π⁰ | π⁰ | ≈ 0.02 | 50–51% | 17–19%, at the MC-statistics floor |
  | multi-π | multi-π | ≈ 2.5 | — | — (not usable as a constraint) |

- The framework's `ConstrainedCalculator` applies sideband reco bins through a conditional
  covariance (DocDB 32672). It is unused so far.
- Fitting tools: ROOT has Minuit2 (OpenMP build) and RooFit; iminuit 2.32 is installable. No
  MCMC code exists yet.

**Protocol constraints.** These come from the frozen procedure, §cr_protocol. Anything beyond
them is an **amendment that must be approved before the regions are opened**.
- Permitted fit today: a one-parameter normalisation of a target class, fitted to FHC and RHC
  simultaneously, accepted by a likelihood-ratio test at p<0.01, and only if the region's shape
  test passes.
  - CC0π is pre-authorised.
  - A π⁰ fit may only be *proposed* after the pre-fit result.
- FHC and RHC share one background model. A class is never corrected in one beam alone.
- Signal contamination (8–24%) is held at the central-value tune in any fit.
- Prohibited: bin-by-bin reweighting, choosing a generator after inspection, any use of the
  signal-region beam-on sample.

**Why a constraint is worth having.** The flux term on the extracted cross section is
17% × (1 + B/S), about 30%, because the background is subtracted with the same flux as the
signal. Reducing the effective background uncertainty is the lever, not the unfolding.

## 2. What two beams actually buy (measure this first)

The intuition is that FHC and RHC separate ν from ν̄ backgrounds and probe different energies. A
first look (Run 1, unweighted) says the composition lever arm is small for this selection:

| Region | ν̄ fraction FHC | ν̄ fraction RHC |
|---|---|---|
| Signal region | 0.16 | 0.18 |
| CC0π | 0.12 | 0.16 |
| π⁰ | 0.13 | 0.13 |

The NuMI RHC flux at MicroBooNE carries a large ν_μ component (2.92 vs 3.52 ×10⁻¹⁰ ν̄_μ). On top
of that, ν̄ cross sections and efficiencies are lower. So ν and ν̄ normalisations fitted
separately would be almost fully anti-correlated.

What the joint fit *does* buy:
1. twice the statistics;
2. the correlated flux universes, which tie the two beams' normalisations together;
3. a consistency test between the beams (a parameter goodness-of-fit test);
4. whatever true-energy difference exists between the beams (to be quantified).

**Step 0 decides the parameterisation from Fisher information, not intuition.**

## 3. Options

| | A. Conditional-covariance constraint | B. Binned profile-likelihood fit | C. Bayesian MCMC |
|---|---|---|---|
| Idea | Gaussian update of the signal-region prediction given control-region data, using the joint covariance | Free class normalisations plus Gaussian-constrained systematic nuisances, fitted to control-region bins of both beams | The same likelihood, sampled with priors on every parameter |
| Tools | `ConstrainedCalculator` (in the framework) | Minuit2 MIGRAD/HESSE/MINOS in C++, iminuit for prototypes | NUTS/HMC (e.g. numpyro/Stan) or affine-invariant (emcee); MaCh3-style Metropolis–Hastings is overkill here |
| Correlations | All, automatically | All, through the covariance or the nuisances | All |
| Interpretable parameters | No: moves the whole prediction | Yes: θ_CC0π, θ_π⁰ with intervals | Yes: full posteriors |
| Fits the frozen protocol | No: it also constrains the common flux normalisation, which the protocol leaves to the signal region | **Yes**: generalises item 4(ii) with the same likelihood-ratio test | Only as a cross-check |
| Propagation into unfolding | Native | Reweight class templates in every universe, then add the fit covariance | Needs a Gaussian summary, then as B |
| Risks | Linear/Gaussian; completeness of the detector covariance; signal-model dependence through the shared flux | Choice of parameterisation; Gaussian nuisance approximation | Priors, convergence diagnostics, cost; hard to review |
| Cost | Low | Moderate | High |

**The methods agree in the limit that matters.** With linear templates, Gaussian systematics and
large counts, profiling B's nuisances analytically gives exactly A's constrained prediction. C's
posterior mode and credible intervals should equal B's best fit and MINOS intervals. The
choice is therefore about **interpretability, governance and robustness**, not the answer.

Minimiser choice within B:
- MIGRAD and HESSE for the fit and its covariance.
- MINOS for asymmetric intervals.
- A scan of the profile-likelihood ratio for the p<0.01 test.
- Global minimisers (simulated annealing, genetic) are unnecessary: the problem is low-dimensional
  and close to quadratic.
- Feldman–Cousins is unnecessary: the parameters sit far from boundaries and the counts are
  large (CC0π holds about 12k events per beam).

## 4. Recommendation

1. **Primary: B.** A joint FHC+RHC binned likelihood.
   - **Statistics term:** combined Neyman–Pearson Poisson (CNP), matching the collaboration's
     unfolding χ².
   - **Systematics:** a full covariance, scaled with the fitted templates, so the nuisances are
     profiled analytically.
   - **Parameters:** a *small* set of physically motivated class normalisations, fixed by step 0.
   - **Acceptance:** the protocol's likelihood-ratio test.
2. **Propagation.** Scale each fitted class in every systematic universe, so the (1 + B/S)
   coherence is kept. Replace the class prior with the fit covariance, as item 4(ii) already
   describes.
3. **Cross-check A.** Run `ConstrainedCalculator` on the same bins. Agreement confirms B's
   linear/Gaussian treatment. Any difference, in particular the flux-normalisation part A would
   also constrain, is reported and not used.
4. **Cross-check C.** Run an HMC fit on the Asimov dataset and one alternative-model fake-data
   set. Posterior intervals should match MINOS within 10%. C becomes primary only if step 0 shows
   the model needs explicit nuisance parameters, such as a principal-component flux
   decomposition plus the GENIE knobs (about 50–100 parameters), where sampling scales better than
   profiling.

## 5. Fit model (starting proposal; step 0 finalises it)

- **Data vector.** For each beam: the CC0π and π⁰ regions in p_μ at the released edges (7 bins),
  plus π⁰ leading-shower energy (10 bins). The cosmic-region total is included as an EXT
  normalisation check (not fitted). Multi-π stays out of the fit (T ≈ 2.5) but inside the
  goodness-of-fit.
- **Prediction.** μ_i = Σ_c θ_c B_ic + S_i + EXT_i + dirt_i, with the class templates B_ic taken
  from the same class definitions as `sb_transfer.C`.
- **Parameters, in order of priority:**
  - θ_CC0π, shared between the beams (pre-authorised form);
  - θ_π⁰ for CC+NC π⁰, only if a separate approval is granted;
  - candidate extensions, each admitted only if step 0 shows it identifiable (correlation below
    0.8 with its partner, expected uncertainty below 50% of its prior):
    - CC0π in 2 p_μ (energy-proxy) bins;
    - separate ν and ν̄ normalisation for CC0π (expected to fail given the table above).
- **Nuisances in the covariance:**
  - PPFX flux, GENIE multisim and unisims, re-interaction;
  - the 8 detector knobs, correlated between control and signal regions;
  - POT 2%, target count 1%;
  - MC statistics (Barlow–Beeston-lite), EXT and dirt statistics.
- **Fixed:** signal contamination at the central-value tune (with its systematics), EXT scaling,
  dirt.

## 6. Work plan

**Step 0: sensitivity and identifiability (MC only, cheap).**
1. Build class-split templates per (beam, region, bin) and per universe, extending
   `sb_protocol.C`/`sb_transfer.C` from totals to bins.
2. Compute the ν̄ fraction, mean true E_ν and class composition, weighted over all runs, for
   each beam and region.
3. Asimov fit for each candidate parameterisation: Fisher matrix, expected σ(θ), parameter
   correlations.
4. **Expected gain in the signal region.** Propagate the Asimov post-fit covariance into the
   one-bin total and one observable. Report the change in the flux term and the total.

*Gate:* if the total uncertainty on the signal-region cross section improves by less than about
3 percentage points, stop. Keep the frozen validation-only use and write that up.

**Step 1: implementation.**
1. A fit library in C++/Minuit2 next to `sb_protocol.C` (it reuses its covariance builder), with a
   thin iminuit prototype.
2. CNP statistics term, profile-likelihood-ratio test, MINOS intervals.
3. Split the background true bins by class in the bin configs (several `1 -1` bins with class
   selections), so the fitted scales can be applied universe by universe. This needs univmake
   rebuilds.

**Step 2: validation on pseudo-data.** All of these must pass before any request.
1. Asimov closure: θ̂ = 1, zero bias.
2. More than 1000 Poisson and systematic toys: pull mean |m| < 0.1, width 0.9–1.1, MINOS coverage
   68 ± 2%.
3. Injection tests:
   - each class shifted ±30%;
   - a shift in one beam only, which must be flagged by the parameter goodness-of-fit test, not
     absorbed;
   - an energy-dependent CC0π shift;
   - one detector knob at 1σ;
   - a flux shape distortion.
4. Signal injection of ±30% signal contamination: θ̂ must move by less than half its uncertainty.
5. Alternative-model fake data (the existing MaRES ± sets; NuWro/GiBUU-like reweights if
   available): measure the bias on the extracted signal-region cross section with and without the
   constraint.
6. Fit FHC-only, RHC-only and jointly, then run the parameter goodness-of-fit test between them.
7. Cross-check A on the Asimov dataset and two injections, and cross-check C on the Asimov
   dataset and one alternative model.

**Step 3: governance.**
1. Write the amendment to §cr_protocol item 4(ii): fit model, parameters, priors, test
   statistic, acceptance thresholds, propagation, and the validation results of step 2.
2. The π⁰ fit stays a separate proposal.
3. Obtain approval, then open the control regions under the amended frozen procedure.

**Step 4 (after approval): apply.**
1. Pre-fit result first, per the existing protocol.
2. Then the fit, the likelihood-ratio test and the propagation.
3. New release.

## 7. Out of scope / open

- **Proton-tagged sample.** No inverted proton-ID region exists, so the ~13% proton mis-tag
  background cannot be constrained (§cr_1p_scope).
- **Inverted pion-PID region.** This is the direct test of p→π misidentification, the
  mechanism behind the 31–34% CC0π transfer uncertainty. It is listed as a prerequisite for
  opening the signal region. Building it is the most effective way to make the CC0π constraint
  bite, and could be done in parallel with step 0.
- **Decisions for the user:**
  - whether an overall flux-normalisation constraint (what A does implicitly) is acceptable in
    principle;
  - whether to pursue the π⁰ fit amendment now or after the pre-fit result.
