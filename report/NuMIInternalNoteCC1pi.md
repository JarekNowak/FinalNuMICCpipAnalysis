# Muon Neutrino Charged Current Single Pion Cross Section on Argon using Run 1 NuMI Data

**Internal Note — Version 0.2 (DRAFT)**

---

> **DRAFT STATUS — READ FIRST**
>
> This note mirrors the structure of the BNB CC1π internal note
> (`BNBInternalNoteCC1pi_V0.91.pdf`, J. P. Detje, 4 Feb 2025) for the NuMI
> Run 1 FHC measurement.
>
> Sections 1–5 describe the analysis **as it is actually configured in the
> code today** and every number in them was read from the source or configs.
> Section 6 (fake-data tests) now carries **real results** from the re-run
> chain: the GENIE closure χ² for all five observables and the four-generator
> comparison on the corrected flux. Section 7 (real-data results) is **still
> empty** — the real beam-on sample has not been restored and the chain has
> not been run on data (see §2.2 and §7).
>
> The single most important change since v0.1 is the **flux-normalisation
> correction** (§2.4): the NuMI flux constant was ~3.48× too large, which had
> made every extracted cross section ~3.5× too small. Every number in §6 uses
> the corrected flux. Nothing in §7 should be quoted as a data result.

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

> **The `onBNB` slot currently holds GENIE fake data.** It is the
> detector-variation central-value sample (a GENIE G18_10a production),
> symlinked as `xsec-ana-genie_detvarCV_fakedata_run1_fhc.root` and entered
> with its native POT (7.631 × 10²⁰) and beam-on-equivalent triggers
> (18 153 256). This sample is statistically **independent** of the numuMC
> overlay that builds the response matrix, so it is a genuine same-generator
> closure rather than a trivial self-closure. It replaced the earlier NuWro
> overlay fake data (`xsec-ana-numi_nuwro_overlay_pion_ntuples_run1_fhc.root`,
> 6.65 × 10²⁰ POT), which carried a real cos θ_μ *shape* difference versus any
> standalone NuWro (§6.4).
>
> The real beam-on sample
> (`xsec-ana-neutrinoselection_filt_run1_beamon_beamgood.root`, 3.283 × 10²⁰
> POT, 7 809 962 triggers) remains **commented out** at
> `configs/file_properties_numi.txt`. The analysis is therefore configured as
> a **fake-data study**, not a real-data measurement. Section 7 cannot be
> filled until this is swapped back and the chain re-run on data.

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

- **Flux normalisation (dominant correction).** The integrated NuMI flux
  constant in the cross-section conversion factor was
  2.372512 × 10⁻⁹ νμ/cm²/POT, taken from the superseded
  `uboone_numi_flux_histograms.root` (which included an unphysical ~25–30 MeV
  artifact spike). The authoritative value from the flux-file author, using the
  new GEANT4 flux integrated above 60 MeV, is
  4.43515 × 10⁻¹⁰ (νμ) + 2.37644 × 10⁻¹⁰ (ν̄μ) = 6.81159 × 10⁻¹⁰ νμ/cm²/POT —
  a factor **3.48 smaller**. Because the cross section divides by the flux, the
  old constant made every extracted cross section ~3.5× too *small*. The
  correction was confirmed four independent ways: standalone NuWro (2.7× low)
  and GENIE (3.2× low) run on the old flux; the flux author's numbers (3.48×);
  and a GENIE closure in which the framework tune scaled by 3.48× matches a
  standalone GENIE on the corrected flux to 3%. Fixed in
  `FiducialVolume.hh` (`integrated_numu_flux_in_FV`). The flux is the
  active-volume-averaged flux used as an approximation for the fiducial-volume
  flux — a documented ~few-% residual; no fiducial-volume correction is applied
  (the target count is done separately in the same conversion factor).
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

Both the NuMI unit factor (now `1e38`) and the bin-width division in
`Unfolder.C` were corrected to this convention; `Unfolder` and `UnfolderNuMI`
now agree bin-for-bin (§8.1).

### 5.1 Unfolding

Wiener-SVD with the filter enabled and second-derivative regularisation
(`Unfold WienerSVD 1 second-deriv`), chosen over D'Agostini after a variant
study — it suppresses the bin-to-bin oscillation arising from sparse Run 1
statistics. D'Agostini with 4 iterations remains available as a cross-check
via `configs/ccpi_xsec_config_numi_range_dagost.txt`.

As in the BNB analysis, the unfolded result must be compared to predictions
smeared by the additional smearing matrix A_C.

## 6 Fake-Data Tests

All results in this section use the corrected flux (§2.4) and are drawn from a
full re-run of `univmake` + `Unfolder`/`UnfolderNuMI` on all five observables.
Predictions are compared in the space of the unfolded result, i.e. after the
additional smearing matrix A_C is applied (§5.1).

### 6.1 GENIE closure test

The fake data is the detector-variation central-value GENIE sample (§2.2),
statistically independent of the numuMC overlay that builds the response
matrix. Treating it as data, subtracting the scaled beam-off sample and
unfolding should recover its own true signal distribution. The closure is
quantified by the χ² between the unfolded result and the A_C-smeared fake-data
truth, using the full covariance:

| Observable | χ² / ndf | p-value |
|---|---|---|
| p_μ | 7.23 / 22 | 0.999 |
| cos θ_μ | 4.40 / 12 | 0.975 |
| cos θ_π | 3.32 / 12 | 0.993 |
| p_π | 0.26 / 5 | 0.998 |
| θ_μπ | 2.26 / 9 | 0.987 |

All five observables close with p > 0.97: the unfolding machinery recovers the
fake-data truth. (The very high p-values reflect the large Run-1 covariance,
not overfitting.) This is a cleaner closure than the NuWro overlay could
provide, since data and response now share the same generator model.

### 6.2 Fake-data normalisation offset

The unfolded fake data sits ~1.29× above the framework GENIE tune (the numuMC
CV prediction), even though both are the same GENIE G18_10a tune. This was
traced to a **flat** difference in events-per-POT between the detVar-CV sample
and the numuMC overlay: the ratio is 1.293 integrated and 1.294 in the
dominant catch-all true bin, i.e. present across *all* event types, with the
signal-bin shape agreeing within Poisson statistics. It is therefore a POT /
flux bookkeeping difference from using the Run-3b detVar-CV sample as Run-1
fake data — consistent with the framework's own caveat that detVar samples
exist only for Run 3b and are applied globally (`SystematicsCalculator.cxx`) —
not a physics or generator difference and not a framework bug.

The familiar "data ≈ 2× above the generators" then decomposes as
**1.29× (fake-data normalisation) × ~1.3× (A_C smearing)**. It can be removed
by entering the fake data with an effective POT of 7.63 × 10²⁰ × 1.293 =
9.87 × 10²⁰ (and rebuilding the universes); this is documented but left
uncorrected in `configs/file_properties_numi.txt`, since it does not affect the
closure and vanishes entirely once the real beam-on sample replaces the
fake data.

### 6.3 Generator comparison

Standalone predictions for the same signal definition and phase space were
generated on the corrected "newg4" flux and overlaid on every unfolded plot.
Each generator was run in the SL7 container environment:

| Generator | Version / tune | CC1π total, cos θ_μ integral [10⁻³⁸ cm²/Ar] |
|---|---|---|
| GENIE | 3.04, tune G18_10a_02_11a | 0.82 |
| GiBUU | local release build, FSI on | 1.13 |
| NuWro | 21.09 | 1.25 |
| NEUT | local build, MPV-3 flux×σ sampling | 1.32 |

The four span 0.82–1.32 (GENIE lowest, NEUT highest) — a healthy generator
spread. Predictions are combined from νμ and ν̄μ runs with the authoritative
flux fractions (Φ_νμ = 4.43515 × 10⁻¹⁰, Φ_ν̄μ = 2.37644 × 10⁻¹⁰) and converted
to per-nucleus 10⁻³⁸ cm²/Ar (× A = 40). Two per-generator normalisation
conventions had to be handled explicitly: NEUT's `Totcrs` is per **nucleon**
(not per nucleus), and GiBUU's `perweight` is already in 10⁻³⁸ cm² units;
getting these wrong produced 40× and 10³⁸× errors respectively before the fix.

The five comparison figures (`unfold_output/plot_{costhetamu,costhetapi,pmu,
ppi,thetamupi}_0.pdf`) show, per observable: the unfolded fake data with its
uncertainty, the four generators and the MicroBooNE tune (all A_C-smeared), the
fake-data truth curve, and a data/MC ratio panel.

### 6.4 NuWro overlay vs standalone NuWro (shape study)

The NuWro overlay previously used as fake data was compared bin-by-bin against a
standalone NuWro 21.09 run on the same flux. Their cos θ_μ distributions differ
in **shape**, not just normalisation: the overlay/standalone ratio swings from
0.64 (backward) to 1.88 (mid-angle) and back to 0.65 (forward), while the
integral ratio is only ~1.1. The cause is *not* the flux — the two samples have
identical signal-weighted mean neutrino energies (1.93 vs 1.88 GeV, 3%) — but a
NuWro version/configuration difference in the muon angular distribution at fixed
energy. This is why the NuWro overlay was retired in favour of the
independent-GENIE fake data (§2.2), and it is a reminder that an overlay fake
sample is not a clean stand-in for an external generator prediction.

## 7 Results

*Empty pending a real-data run.* The fake-data machinery is now validated
end-to-end (§6): the flux normalisation is corrected, the GENIE closure passes
for all five observables, and four standalone generators are overlaid. What
remains before this section can carry a **data** measurement is, in order:

1. Restore the real beam-on sample in `file_properties_numi.txt` (currently the
   GENIE detVar-CV fake data occupies the `onBNB` slot, §2.2).
2. Re-run the full chain from `ProcessNTuples` — the corrected momentum branches
   (§2.4) do not exist in the currently processed files.
3. Restore the NuMI beamline-geometry systematic (§4.6, §8.2); the uncertainty
   budget is under-estimated without it.

The units/normalisation questions (§8.1) and the flux scale (§8.5) are resolved.

This section should then contain, per observable and for the total: the
unfolded cross section with the full uncertainty breakdown, and the
generator comparison of §6.3 against real data rather than fake data.

## 8 Open Items

**8.1 Cross-section units and bin-width normalisation. — RESOLVED.**
`conversion_factor()` divided by `1e39` in NuMI mode against `1e38` for BNB,
and `Unfolder.C` omitted the bin-width division that `UnfolderNuMI.C` performs,
so the two binaries disagreed on the same bin. Both are now fixed to the BNB
convention (10⁻³⁸ cm²/Ar, differential): the NuMI factor is `1e38`, and
`Unfolder.C` divides each bin by `conv_factor × bin_width × other_var_width`
via `SliceHistogram::transform()`, which also converts the covariance.
Verified on the costhmu NuWro config that `Unfolder` and `UnfolderNuMI` now
agree bin-for-bin (ratio 1.0000 across all 12 bins). Relative to earlier
output, every cross-section number is ~10× smaller and now differential.

**8.2 Beamline geometry systematic.** §4.6.

**8.3 Anti-neutrino contamination. — CONSISTENT, by construction.** The signal
uses `|mc_nu_pdg| == 14`, admitting ν̄μ CC events into signal and the efficiency
denominator (the three sibling selections use an exact comparison). With the
flux fix (§2.4) the normalising flux is now the **combined** νμ + ν̄μ flux
(4.43515 × 10⁻¹⁰ + 2.37644 × 10⁻¹⁰), and the standalone generators are combined
with the same fractions (§6.3). Signal, flux and predictions therefore all
treat νμ and ν̄μ together, so there is no longer a normalisation mismatch — the
measured quantity is an effective flux-averaged νμ+ν̄μ CC1π cross section. A
νμ-only measurement would instead require restricting the signal to
`mc_nu_pdg == 14` and using the νμ-only flux; that is a definition choice, not a
bug.

**8.4 Unused muon BDT.** `tmvaReader_mu` is constructed and booked but never
evaluated, and its input variables are never filled — suggesting an intended
muon-PID cut was dropped.

**8.5 Flux normalisation. — RESOLVED.** The NuMI flux constant was ~3.48× too
large (§2.4), making cross sections ~3.5× too small. Corrected to the flux
author's authoritative value (6.81159 × 10⁻¹⁰ νμ/cm²/POT, E > 60 MeV) in
`FiducialVolume.hh` and confirmed four independent ways. A residual ~few-%
active-volume-vs-fiducial-volume approximation is documented in that file.

**8.6 Fake-data normalisation offset. — DOCUMENTED, uncorrected.** The
detVar-CV fake data sits a flat ~1.29× above the numuMC GENIE tune from a
Run-3b-vs-Run-1 POT/flux bookkeeping difference (§6.2). It does not affect the
closure and disappears when the real beam-on sample replaces the fake data; the
POT correction to remove it (effective POT 9.87 × 10²⁰) is noted in
`file_properties_numi.txt`.

---

*Prepared as a companion to `BNBInternalNoteCC1pi_V0.91.pdf`. Analysis code
state: see git history on branch `fix/systematics-warning`. Fake-data figures:
`unfold_output/plot_*_0.pdf` and the per-observable closure dumps
`unfold_output/closure_hists_ccpi_Run1_*_xsec.root`. Standalone generator
predictions and their build scripts: `generator_predictions/newg4/`.*
