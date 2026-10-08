# Phase 4 (first part): feasibility of the sideband constraint

2026-10-07. Phase 4 of `report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md`, items 1, 2 and 4, done as a study in Python
before the framework implementation (item 3). Step `constraint` of `scripts/mp_evt.py`; tables in
`phase4_constraint.md` (deployed particle classifier) and `phase4_constraint_llronly.md` (LLR-and-geometry classifier),
both at score > 0.92 with the PID event features.

## Method

The formalism of the framework's `ConstrainedCalculator`. Every systematic source changes the predicted signal plus
background in each region: the cross-section and reinteraction universes through the efficiency (smearceptance times the
CV true rate) and the background, the flux universes (`weight_ppfx_all`, 600) through the selected signal rate and the
background, the detector variations through the ratios of the Run-4 variation samples per horn mode (as in `detsyst`).
With POT and targets, the data statistics of every region and the MC statistics, the covariance of the regions is
conditioned on the sidebands; the one-bin uncertainty is sqrt(C(SR | sidebands)) / S, and the background estimate is
B_CV + K (N_sidebands - N_sidebands,CV). Test periods projected to the full exposure (COMB). Without a constraint the
same machinery reproduces the totals of Phase 2 (47.0% against 47.4%; the flux term from the PPFX universes, 32%,
against 31% from 17% x (1 + B/S)).

## Regions

Exclusive, inside the two-pion domain, in this order: SR (score > 0.92); three-pion-like (at least four tracks, three
pions more likely than two; the control sample of D8); proton-like (an assigned pion with P(p) > 0.5); score 0.5 to 0.92;
score below 0.5. The sidebands leave out the events of the blind single-pion signal region (6–10% of the proton-like and
score-band events; 0.4% of the three-pion-like ones), so they can be opened without it. Deployed classifier (COMB, full
exposure):

| Region | Events | Signal | Largest classes |
|---|---|---|---|
| SR | 896 | 499 (56%) | CC3π+ 13%, CC1π 11% |
| Three-pion-like | 220 | 29 (13%) | CC3π+ 38%, CCπ0 18% |
| Proton-like | 14,868 | 430 (3%) | CC0π 41%, CCπ0 16%, CC1π 16% |
| Score 0.5 to 0.92 | 3,878 | 582 (15%) | CCπ0 28%, CC1π 22% |
| Score below 0.5 | 10,889 | 140 (1%) | beam-off 46%, CC0π 11% |

Of the SR, 14% lies in the blind single-pion signal region (D7) and 12% in the opened multi-π control region (D3).

## Constraint (COMB, one-bin total)

| Constraint | Deployed: total | detector | flux | LLR and geometry: total | detector | flux |
|---|---|---|---|---|---|---|
| None | 47.0% | 30.6% | 32.3% | 45.5% | 20.6% | 35.4% |
| Three-pion-like | 25.0% | 17.2% | 12.6% | 24.2% | 11.7% | 13.6% |
| Proton-like | 36.0% | 29.3% | 11.5% | 31.9% | 22.9% | 13.0% |
| Score 0.5 to 0.92 | 32.8% | 29.5% | 6.6% | 22.5% | 17.2% | 6.2% |
| Score below 0.5 | 43.7% | 36.7% | 20.8% | 33.0% | 26.2% | 13.6% |
| All sidebands | 20.8% | 15.1% | 6.0% | 16.8% | 8.8% | 5.5% |

- The flux term falls from 32–35% to 6% with all sidebands: the flux scales the background and the sidebands alike.
- The detector term falls from 31% to 15% (deployed) and from 21% to 9% (LLR and geometry); what remains is close to
  the efficiency part found in the robustness study (14% and 8%), which no background sideband constrains.
- The three-pion-like sample alone (220 events) halves the total; its data statistics (6.7%) raise the data term from
  6% to 10%.
- With the constraint the LLR-and-geometry classifier is ahead (16.8% against 20.8%).

## Alternative-model closure (item 4, first part)

Pseudo-data from a reweighted model in every region, extracted with the CV model; extracted over true signal rate:

| Constraint | GENIE universe 545 | Δ→Nπ angular | GENIE multisim: rms bias | within the total |
|---|---|---|---|---|
| None (deployed) | 0.993 | 1.019 | 12.7% | 100% |
| All sidebands (deployed) | 0.983 (−0.08σ) | 1.003 (+0.01σ) | 4.0% | 100% |
| None (LLR and geometry) | 0.979 | 1.012 | 17.4% | 98% |
| All sidebands (LLR and geometry) | 0.952 (−0.29σ) | 0.994 (−0.03σ) | 4.5% | 100% |

Every single-sideband constraint closes as well (within 0.17σ). These models lie inside the variations that build the
covariance, so the test checks the procedure, not the model space; the independent-generator test (NuWro, D5) is open.
The score-band sideband holds more signal than the SR (582 against 499 events), so its constraint carries signal-model
information into the background estimate; it closes for these two models.

## Framework implementation (item 3, 2026-10-08)

- `CC1mu2piBDT` (`NewPIDSelections.hh`): the assignment, the 33 features, the event-classifier score (RBDT,
  `booster_decision_tree/mp_pid/evt_clf_n2_llronly_pid_rbdt.root`) and the region in the selection, on the
  LLR-and-geometry particle classifier; Selected is the signal region. On the 30k slice the features equal the Python
  ones exactly, the scores to 6e-8, the regions in every event; `CC1mu1piXp` stays bitwise identical to version 1.7.
- Production of the 28 release and 18 detector-variation samples (`slurm/slurm_mp2pi.sbatch`,
  `/data/uboone/processed/mp2pi`); file lists, links and per-run pseudo-data with dirt (`scripts/mp2pi_setup.py`,
  `macros/throw_mp2pi.C`; the framework adds the beam-off); bin configuration with the SR and the four sidebands
  (`configs/mp2pi_total_bin_config.txt`, sideband bins of type 1 in block 0 so that their signal is predicted).
- `Calculator Constrained` in the extraction configuration selects the `ConstrainedCalculator`. The extractor now
  propagates the ordinary block of each covariance, and with the constraint what remains of it after the conditioning
  on the sidebands, so the per-source terms add up to the constrained total; before, configurations with sideband
  bins returned zero uncertainties.

| One-bin total | Total | Detector | Flux | Cross section | Data statistics | Unfolded / pseudo-data truth |
|---|---|---|---|---|---|---|
| FHC, no constraint | 58.4% | 40.8% | 37.0% | 14.7% | 10.5% | 1.055 |
| FHC, constrained | 24.9% | 15.4% | 6.1% | 9.2% | 13.9% | 1.013 |
| RHC, no constraint | 46.6% | 28.4% | 31.5% | 15.5% | 8.7% | 1.035 |
| RHC, constrained | 18.5% | 11.1% | 5.0% | 6.1% | 11.2% | 1.015 |
| COMB, no constraint | 50.7% | 33.4% | 33.7% | 15.0% | 6.7% | 1.044 |
| COMB, constrained | 17.9% | 11.9% | 4.7% | 6.7% | 8.7% | 1.019 |

With the constraint the combined total agrees with the Python study (17.9% against 16.8%). Without it the framework
is higher (50.7% against 45.5%) through the detector term, which differs in RHC (28.4% against 17.1%; FHC 40.8%
against 38.0%).

## Detector variations: normalisation of the samples (2026-10-08)

The difference comes from the detector-variation samples. Per 1e20 POT, the Run-4 detector-variation CV overlays
predict 1.27 (FHC) and 1.51 (RHC) times the event rate of the nominal Run-4 overlays (RHC single-pion signal 781
against 517, all numu 28585 against 19906); their summed_pot equals the raw SubRun POT sum, so the samples, not the
processing, differ. The FHC and RHC variation samples have nearly the same rate and antineutrino fraction (0.20),
unlike the nominal modes (0.16 and 0.21). The framework's `DV` covariance takes the absolute difference between a
variation and the variation CV, scaled to the data POT, so it carries that rate ratio; the Python study takes the
ratio of the variation to the variation CV and applies it to the nominal prediction.

`DVfrac` in the systematics configuration (`SystematicsCalculator::make_cov_mat_frac`) does the latter: in every
covariance bin the difference is nominal x (variation / variation CV - 1). `configs/ccpi_systcalc_numi_dvfrac.conf`
is the release systematics configuration with the eight detector lines as `DVfrac`; the `*_dvfrac` two-pion
configurations use it.

| One-bin total, DVfrac | Detector | Total |
|---|---|---|
| FHC, no constraint | 29.8% | 51.3% |
| FHC, constrained | 11.9% | 22.4% |
| RHC, no constraint | 17.9% | 41.0% |
| RHC, constrained | 7.6% | 16.3% |
| COMB, no constraint | 23.0% | 44.5% |
| COMB, constrained | 8.5% | 15.7% |

The RHC detector term now agrees with the Python study (17.9% against 17.1%), and the combined constrained total
with it (15.7% against 16.8%).

### Impact on the released single-pion results (study only; version 1.7 unchanged)

`scripts/dvfrac_impact.sh` runs every released single-pion extraction (65: inclusive and proton-tagged, the one-bin
totals, the two two-dimensional results and W_πp) twice on copies of the release universe files, as released and with
`DVfrac`; `scripts/dvfrac_compare.py` writes `dvfrac_impact.md` and `dvfrac_impact.tsv` (this folder). The run as
released reproduces the released SYSTDUMP of all 65.

| One-bin total | Detector, released | DVfrac | Total, released | DVfrac |
|---|---|---|---|---|
| Inclusive FHC | 19.2% | 14.2% | 38.1% | 35.9% |
| Inclusive RHC | 22.1% | 14.3% | 38.1% | 34.1% |
| Inclusive COMB | 20.4% | 14.0% | 37.7% | 34.6% |
| Proton-tagged FHC | 11.7% | 8.6% | 41.4% | 40.6% |
| Proton-tagged RHC | 26.7% | 16.8% | 47.6% | 42.9% |
| Proton-tagged COMB | 16.2% | 10.9% | 42.3% | 40.5% |

- Per bin the detector uncertainty falls by a median factor 0.78 (FHC), 0.70 (RHC) and 0.71 (COMB), close to the
  inverse rate ratios of the variation samples (0.79 and 0.66). The total falls by a median factor 0.94; most where
  the detector term dominated (proton-tagged FHC δφ_T 112% to 47%, δp_T 83% to 50% on the integrated cross section;
  COMB p_n 97% to 67%).
- The central values of the differential results move, because the Wiener filter is built from the covariance:
  the integrated cross section of the 59 differential results changes by a median factor 1.026 (0.92 to 1.18; 16
  beyond 5%, the largest proton-tagged RHC θ_p), the one-bin totals not at all (no filter). The trace of A_C rises in
  most results (inclusive COMB cos θ_μ 3.81 to 4.23, θ_μπ 2.51 to 3.02). The largest bin change is the third
  proton-tagged FHC δφ_T bin, 0.024 ± 0.061 to 0.051 ± 0.031.
- The truth χ² of the pseudo-data rises slightly with the smaller covariance; the largest is 7.28/4 (proton-tagged
  RHC cos θ_π, 5.41/4 as released).

## Open

- Item 4, rest: CV pseudo-data and ensembles through the framework; NuWro if a NuMI sample can be validated.
- Item 5: the same constraint for the single-pion analysis, on pseudo-data.
- The choice of particle classifier (deployed or LLR and geometry), deferred to after this phase.
