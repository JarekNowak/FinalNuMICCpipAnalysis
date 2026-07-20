# Muon Neutrino Charged Current Single Pion Cross Section on Argon using Run 1 NuMI Data

**Internal Note — Version 0.1 (DRAFT)**

---

> **DRAFT STATUS — READ FIRST**
>
> This note mirrors the structure of the BNB CC1π internal note
> (`BNBInternalNoteCC1pi_V0.91.pdf`, J. P. Detje, 4 Feb 2025) for the NuMI
> Run 1 FHC measurement.
>
> Sections 1–5 describe the analysis **as it is actually configured in the
> code today** and every number in them was read from the source or configs.
> Sections 6–7 (fake-data tests and results) are **structural placeholders**.
> No cross-section values, efficiencies, purities, χ² values or figures appear
> anywhere in this note, because the analysis has not been re-run since the
> corrections described in Section 2.4. Nothing here should be quoted as a
> result.
>
> Three things must be settled before this note can carry results at all —
> see Section 8.

---

## 1 Overview

This note describes the extraction of the νμ charged-current single charged
pion cross section on argon using MicroBooNE Run 1 FHC data from the NuMI
beam. It is the NuMI counterpart of the BNB analysis of the same final state,
and is intended to be read alongside that note, which contains fuller
descriptions of the selection variables and the boosted decision trees.

The measurement is performed with the common `xsec_analyzer` framework: event
selection via the `CC1mu1piXp` selection class, universe construction via
`univmake`, and unfolding via the Wiener-SVD implementation in
`CrossSectionExtractor`. Section 2 describes the software and samples,
Section 3 the signal definition and selection, Section 4 the uncertainty
treatment, Section 5 the cross-section extraction and unfolding, Section 6 the
fake-data tests, and Section 7 the results.

The principal differences from the BNB analysis are:

| | BNB analysis | This analysis |
|---|---|---|
| Beam | BNB | NuMI FHC |
| Runs | 1–5 | 1 only |
| Flux systematic | `weight_flux_all` | `weight_ppfx_all` (PPFX) |
| Beamline geometry systematic | n/a | NuMI-specific — **not currently included**, see §4.6 |
| Proton multiplicity | — | Xp (any number, no threshold) |
| Unfolding | Wiener-SVD | Wiener-SVD (D'Agostini available as a cross-check) |

## 2 Software and Samples

### 2.1 Framework

The analysis runs entirely within `xsec_analyzer`. The chain is:

1. `ProcessNTuples` — applies the `CC1mu1piXp` selection to raw PeLEE ntuples
   and writes `stv_tree` output (`run_process_ntuples.sh`).
2. `univmake` — builds systematic universe histograms from a bin config
   (`run_universe_maker.sh`).
3. `Unfolder` / `UnfolderNuMI` — extracts and unfolds the cross section
   (`run_unfolder.sh`).

### 2.2 Samples

Run 1 FHC, as configured in `configs/file_properties_numi.txt`:

| Sample type | Purpose |
|---|---|
| `onBNB` | beam-on data (see warning below) |
| `extBNB` | beam-off (EXT), 3 821 593 triggers |
| `numuMC` | NuMI overlay CV |
| `dirtMC` | dirt overlay |
| `detVarCV` + 8 variations | detector systematics |

> **The `onBNB` slot currently holds NuWro fake data**
> (`xsec-ana-numi_nuwro_overlay_pion_ntuples_run1_fhc.root`, 6.65 × 10²⁰ POT).
> The real beam-on sample
> (`xsec-ana-neutrinoselection_filt_run1_beamon_beamgood.root`, 3.283 × 10²⁰
> POT, 7 809 962 triggers) is **commented out** at
> `configs/file_properties_numi.txt:31`.
>
> The analysis is therefore configured as a **fake-data study**, not a
> real-data measurement. Section 7 cannot be filled until this is swapped
> back and the chain re-run.

### 2.3 Detector variation samples

Eight variations are used: `detVarLYdown`, `detVarLYrayl`, `detVarRecomb2`,
`detVarSCE`, `detVarWMAngleXZ`, `detVarWMAngleYZ`, `detVarWMX`, `detVarWMYZ`.

The BNB note additionally lists `detVarLYatten`; it is not present in the NuMI
sample set here. Like the BNB analysis, the WireMod dE/dx variation is not
used, in favour of the recombination variation.

### 2.4 Corrections applied since the previous iteration

The following were found and fixed during a code review of this analysis. All
of them affect the numbers a run would produce, which is why no results from
before them should be carried forward:

- **Muon momentum estimator.** The estimator combining range (contained) and
  MCS (uncontained) was computed but never written to the output tree, so the
  unfolding binned on MCS alone. It is now branched as
  `candidate_muon_mom_reco` and `ccpi_pmu_bin_config.txt` uses it.
- **Pion momentum estimator.** Reco pion momentum was a proton-hypothesis
  range *kinetic energy* binned against a true *momentum* on identical bin
  edges. It is now the muon-hypothesis range momentum.
- **Per-category systematics.** Six `MCFullCorrCategory` entries evaluated to
  identically zero and have been removed; see §4.6.
- **Background true bins.** `MakeConfig` had stopped emitting them, so any
  regenerated bin config would have had no background subtraction at all.
- **Diagnostic histograms** are now filled with the CV weight
  (`spline_weight_ × tuned_cv_weight_`).

## 3 Signal and Selection

### 3.1 Signal definition

An event is signal if, in truth:

- the neutrino vertex is inside the fiducial volume:
  10 < x < 246 cm, −101 < y < 101 cm, 10 < z < 986 cm;
- the interaction is charged-current;
- the neutrino is νμ (`|mc_nu_pdg| == 14`, see §8.3);
- there is exactly one muon above threshold;
- there is exactly one charged pion;
- there are no neutral pions, no kaons and no heavier mesons;
- any number of protons is allowed, with no momentum threshold (the "Xp" in
  `CC1mu1piXp`);

and the event falls in the measured phase space:

- p_μ > 0.15 GeV/c
- p_π > 0.175 GeV/c
- θ_μπ < 2.6 rad

The phase-space cuts are applied identically in `define_signal()` and
`categorize_event()`, so signal and category assignment cannot disagree.

### 3.2 Binning

Five differential observables are measured. Bin edges as configured:

| Observable | Config | Bins (true/reco) | Edges |
|---|---|---|---|
| p_μ | `ccpi_pmu_bin_config.txt` | 23/22 | 0.15, 0.20, 0.30 … 2.00, 2.30, 2.60, 3.00 GeV/c |
| p_π | `ccpi_ppi_bin_config.txt` | 6/5 | 0.175, 0.20, 0.30, 0.40, 0.50, 1.00 GeV/c |
| cos θ_μ | `ccpi_costhmu_bin_config.txt` | 13/12 | −1.0, −0.5, 0.0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.82, 0.88, 0.93, 0.97, 1.0 |
| cos θ_π | `ccpi_costhpi_bin_config.txt` | 13/12 | −1.0, −0.7, −0.4, −0.2, 0.0, 0.2, 0.35, 0.5, 0.65, 0.78, 0.88, 0.95, 1.0 |
| θ_μπ | `ccpi_thmupi_bin_config.txt` | 10/9 | 0.0 … 2.513, 2.600 rad |

Each true-bin count is one greater than the reco count: the extra true bin is
the background bin, defined by inverting the signal definition
(`!CC1mu1piXp_MC_Signal`). All configs are single-block with no sideband reco
bins.

### 3.3 Particle identification

Three BDTs are used, trained separately and applied via `TMVA::Reader`:

| Reader | Weights (`booster_decision_tree/`) | Role |
|---|---|---|
| `tmvaReader` | `dataset_MIP_BDT_no_len` | MIP identification |
| `tmvaReader_mu` | `dataset_muon_BDT` | muon — **booked but never evaluated**, see §8.4 |
| `tmvaReader_pi` | `dataset_pion_BDT_no_len` | pion |

Input variables: the three-plane Bragg likelihood ratios
(`trk_bragg_p_v`, `trk_bragg_mu_v`, `trk_bragg_mip_v`), the LLR PID score,
the track score, track length, and the SCE-corrected track end position. The
pion reader omits track length.

### 3.4 Event selection cuts

Applied in order; an event passes only if all are satisfied:

1. software trigger
2. reconstructed vertex in the fiducial volume
3. topological score
4. muon candidate is track-like
5. pion candidate is contained
6. muon candidate not in a detector gap
7. pion candidate not in a detector gap
8. shower cut
9. muon–pion opening angle cut
10. final track-multiplicity cut

### 3.5 Selection efficiency and purity

*To be filled from a complete `ProcessNTuples` pass.* `finalize()` prints
efficiency and purity; these are now CV-weighted (§2.4) and so differ from any
previously recorded values.

## 4 Uncertainties

Configured in `configs/ccpi_systcalc_numi.conf`. The total is

```
total = detVar_total + flux_total + reint + xsec_total
        + POT + numTargets + SimulationStats + DataStats
```

### 4.1 Flux

PPFX hadron-production multisims via `weight_ppfx_all`. A 2% fully-correlated
normalisation uncertainty is applied for POT (`POT MCFullCorr 0.02`).

### 4.2 Interaction

GENIE multisims (`weight_All_UBGenie`) plus nine unisim knobs
(`AxFFCCQEshape`, `DecayAngMEC`, `NormCCCOH`, `NormNCCOH`, `RPA_CCQE`,
`ThetaDelta2NRad`, `Theta_Delta2Npi`, `VecFFCCQEshape`, `XSecShape_CCMEC`) and
the two second-class-current form factors (`xsr_scc_Fa3_SCC`,
`xsr_scc_Fv3_SCC`), matching the BNB treatment.

### 4.3 Reinteraction

`weight_reint_all` multisims.

### 4.4 Detector

The eight detector variation samples of §2.3, evaluated as `DV` type.

### 4.5 Target

A 1% normalisation uncertainty on the number of argon nuclei in the fiducial
volume (`numTargets MCFullCorr 0.01`), as in the BNB analysis.

### 4.6 Known gaps in the uncertainty budget

Two components are **not** included, and the total above is correspondingly
under-estimated:

**NuMI beamline geometry.** The beamline variation weights
(`weight_Horn_2kA`, `weight_Horn1_x_3mm`, `weight_Beam_spot_*`, etc.) are
absent from the processed ntuples. `instructions_numi.txt` documents
`AddBeamlineGeometryWeights` as the tool that injects them, but no driver
script invokes it, and the required histogram file
(`NuMI_Geometry_Weights_Histograms.root`, canonical copy at
`/exp/uboone/data/users/pgreen/NuMIFlux/NewFluxFiles/`) is not available
locally. This systematic has no BNB counterpart and must be restored before
the measurement is complete.

**Background normalisation by category, and dirt normalisation.** Six
`MCFullCorrCategory` entries were removed because they evaluated to identically
zero — they matched true bins by exact string equality against `"category == N"`
while the bin configs declare a single `"!CC1mu1piXp_MC_Signal"` background
bin. Their category indices were also inherited from the NuMI CC1e analysis and
did not correspond to this selection's categories. Removing them changed no
number (verified: 152 histograms bit-identical). Full detail in
`configs/ccpi_systcalc_numi.README.md`.

Note the BNB uncertainty breakdown (Fig. 43–44 of that note) comprises
`EXTstats, MCstats, POT, detVar_total, flux, numTargets, reint, xsec_total` —
i.e. the same set retained here, with no per-category or dirt term.

## 5 Cross-section Extraction

The flux-integrated differential cross section in truth bin α is obtained by
subtracting background from the measured event rate, unfolding to truth space,
and dividing by the integrated flux Φ, the number of argon targets T, and the
bin width:

```
dσ/dx |α  =  (1 / (T · Φ · Δx_α)) · Σ_a U_αa (D_a − B_a)
```

where `D_a` is the measured rate in reco bin a, `B_a` the estimated background,
and `U_αa` the unfolding matrix, which absorbs both the reco→truth mapping and
the efficiency.

Targets are counted as **argon nuclei** in the fiducial volume, so the result
is per-nucleus, not per-nucleon — matching the BNB convention
(`T = 1.03 × 10³⁰` nuclei there). For NuMI the fiducial volume additionally
excludes the dead region 675.1 < z < 775.1 cm.

Results are quoted in units of **10⁻³⁸ cm² / Ar**, differential in the
observable (so 10⁻³⁸ cm² / GeV / Ar for momenta, 10⁻³⁸ cm² / rad / Ar for
angles).

> **Unresolved — see §8.1.** The code currently divides by `1e39` in NuMI
> mode and `1e38` otherwise, and `Unfolder.C` does not divide by bin width
> while `UnfolderNuMI.C` does. Both must be settled before any number in
> Section 7 is meaningful.

### 5.1 Unfolding

Wiener-SVD with the filter enabled and second-derivative regularisation
(`Unfold WienerSVD 1 second-deriv`), chosen over D'Agostini after a variant
study — it suppresses the bin-to-bin oscillation arising from sparse Run 1
statistics. D'Agostini with 4 iterations remains available as a cross-check
via `configs/ccpi_xsec_config_numi_range_dagost.txt`.

As in the BNB analysis, the unfolded result must be compared to predictions
smeared by the additional smearing matrix A_C.

## 6 Fake-Data Tests

*Structure only — no results.*

### 6.1 GENIE closure test

Treat the GENIE CV as data, add the proportionally scaled beam-off sample,
unfold, and confirm the GENIE CV truth distribution is recovered. Driven by
`configs/ccpi_xsec_config_numi.txt`.

### 6.2 NuWro fake-data study

Same procedure with the NuWro overlay
(`configs/ccpi_xsec_config_numi_nuwro.txt`, driven by
`run_nuwro_observables.sh` across all five observables). As in the BNB
analysis, only statistical and model uncertainties apply, since the samples
differ only in generator model.

*Per-observable results to be added: reco-space selection histograms, reco-
and truth-space bin correlations, confusion matrices, unfolded cross sections
and additional smearing matrices, for each of cos θ_μ, cos θ_π, p_μ, p_π,
θ_μπ.*

## 7 Results

*Empty pending a real-data run.* Requires, in order: the real beam-on sample
restored in `file_properties_numi.txt`; the full chain re-run from
`ProcessNTuples` (the corrected momentum branches do not exist in the current
processed files); and §8.1 resolved.

Should contain, per observable and for the total: the unfolded cross section
with the full uncertainty breakdown, and generator comparisons.

## 8 Open Items

**8.1 Cross-section units and bin-width normalisation.** `conversion_factor()`
divides by `1e39` in NuMI mode against `1e38` for BNB, while the BNB note
quotes 10⁻³⁸ cm²/Ar — a factor 10 discrepancy in convention. Separately,
`UnfolderNuMI.C` divides by bin width and `Unfolder.C` does not, so the two
binaries disagree on the same bin. The BNB note is unambiguous that the result
is differential and per-Ar.

**8.2 Beamline geometry systematic.** §4.6.

**8.3 Anti-neutrino contamination.** The signal uses
`|mc_nu_pdg| == 14`, admitting ν̄μ CC events into signal and the efficiency
denominator. All three sibling selections in the framework use an exact
comparison. If the flux normalisation is νμ-only, this biases the result.

**8.4 Unused muon BDT.** `tmvaReader_mu` is constructed and booked but never
evaluated, and its input variables are never filled — suggesting an intended
muon-PID cut was dropped.

---

*Prepared as a companion to `BNBInternalNoteCC1pi_V0.91.pdf`. Analysis code
state: see git history on branch `fix/systematics-warning`.*
