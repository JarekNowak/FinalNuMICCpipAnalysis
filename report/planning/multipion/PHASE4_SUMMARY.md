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

## Open

- Item 3: the framework implementation. The assignment and the event classifier have to run in the selection (C++,
  RBDT for the event classifier) so that bin configurations with the SR and sideband reco bins can go through univmake
  and a `ConstrainedCalculator` switch in `UnfolderNuMI`.
- Item 4, rest: CV pseudo-data and ensembles through the framework; NuWro if a NuMI sample can be validated.
- Item 5: the same constraint for the single-pion analysis, on pseudo-data.
- The choice of particle classifier (deployed or LLR and geometry), deferred to after this phase.
