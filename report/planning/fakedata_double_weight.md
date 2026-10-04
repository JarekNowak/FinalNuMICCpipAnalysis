# The fake data is CV-weighted twice

Found 2026-08-30 while checking whether the unfolded results are compared against
the correct input.

## What the fake data is

Not a detVar-CV throw. `macros/throw_perrun_fhc.C` throws it from the **numuMC
overlay**, per run:

    cv = tuned_cv_weight * ppfx_cv_weight * normalisation_weight
    n  = Poisson( cv * D_run / MCPOT_run )     copies of each event

and writes those three weights back as **1.0** in the output file. The CV weight is
therefore already embedded in *how many events exist*.

## The double count

`SystematicsCalculator::build_universes()` treats a fake-data `onBNB` sample by
substituting the **weighted** CV histogram for the raw counts
(`dataContainsWeightedCV`). For a fake-data file that histogram is the thrown events
re-weighted by the CV weight a second time.

Measured in the per-file subdirectories of the univmake:

| config | thrown (unweighted_0_reco) | used as data (TunedCentralValue_0_reco) | factor |
|---|---|---|---|
| FHC5 | 1621.0 | 1721.9 | 1.062 |
| RHCFULL | 1726.0 | 1860.7 | 1.078 |
| COMB | 3347.0 | 3582.7 | 1.070 |

Cross-checks that this is real, not a bookkeeping artefact:

- 1621.0 equals a direct count of `CC1mu1piXp_Selected` in the four FHC fake-data
  files (558 + 242 + 409 + 412).
- 1721.9 equals `FakeDataMC_0_reco`, and equals `onBNB_reco - extBNB_reco`.
- In the run-1 subdirectory alone: `unweighted_0_reco` 558.0 with 558 entries
  (weight exactly 1), `weight_TunedCentralValue_UBGenie_0_reco` 599.9 with the same
  558 entries -- a mean weight of 1.075 applied to events that should carry 1.
- The throw's own expectation, `sum(tcv*pcv*nw) * D/MCPOT` over selected events
  summed across the four runs, is 1567.1 -- consistent with the 1621 thrown
  (1.3 sigma on Poisson) and inconsistent with the 1721.9 used.

## What it explains

The closure offset measured with the two least-distorting configurations is
1.127 (identity regularisation) and 1.125 (D'Agostini). The double weighting
accounts for 1.062-1.078 of that. It is the largest single contribution but not
all of it: roughly 5-6% remains after it is removed.

## Excluded along the way

- POT scaling: the throw's hardcoded MC POT matches each file's `summed_pot` to
  better than 0.01%.
- Stale inputs: the fake data (08-21) is newer than the MC it was thrown from (08-19).
- The CV weight definition: `apply_cv_correction_weights` gives the UBGenie universe
  `tune * ppfx * normalisation`, the same product the throw uses.
- The comparison target: `h_fakedata_truth` comes from the same files, and
  `trans_mat` is applied to unfolded and truth alike.

## What to do

Either the throw should not fold the CV weight into the event multiplicity, or
`build_universes` should use `unweighted_0_reco` for a fake-data `onBNB` sample
whose weights are already 1. The two conventions are individually defensible and
are being combined.

Until this is settled, the integrated cross sections carry a ~6-8% normalisation
error from this source alone, on top of the residual offset.

---

## Follow-up: is the proton-tagged residual caused by its own throw?

Asked 2026-08-30, after the fix left the inclusive family at unity but the
proton-tagged family at ~1.09.

**No. The two families are thrown from the same events.**

`throw_perrun_w.C` reads the same source ntuples with the same per-run POT values
as `throw_perrun_fhc.C` (3.283e20 / 2.3282e21 for run 1, and so on). The outputs
were checked directly: `xsec-ana-fakedata_fhc_run1.root` and
`w/xsec-ana-fakedata_fhc_run1.root` both hold 100274 entries and the same
`sum(reco_nu_vtx_sce_x)` to full floating-point precision. They are the same Poisson
realisation, processed with different selection branches (`CC1mu1piXp_Selected` in
one, `CC1mu1pi1p_Selected` in the other), which is why they cannot simply be
swapped.

Also excluded:

- **POT convention.** `processed/` and `processed/w/` carry identical `summed_pot`
  for the same run (2.3282e21, 2.4934e21), so the two-convention issue does not
  bite here.
- **Regularisation.** Running the proton-tagged family with identity regularisation
  gives 1.067 / 1.119 / 1.082 (mean 1.090 over 33 extractions), against 1.053 /
  1.118 / 1.100 with the second derivative. The residual is unchanged, so it is not
  a smoothing artefact -- unlike the observable-to-observable scatter, which is.
- **The double weighting itself.** It was present in the proton-tagged samples too,
  but smaller (factor 1.028 against 1.062 for the inclusive), and removing it moved
  the closure from ~1.15 to ~1.09.

So a ~9% residual remains, specific to the `CC1mu1pi1p` chain, and it is not the
fake data, the POT, or the regularisation. What is left is the selection itself:
its efficiency, its response matrix, or its background subtraction. That is where
to look next.

---

## Localising the proton-tagged residual: it is the input, fluctuating high

Decomposed 2026-08-30, FHC5, `costhmu`, proton-tagged.

Both ends of the chain check out independently:

- **The throw is correct.** Run 1 holds 252.0 selected events against a tree-level
  expectation of 251.3, computed as `sum(tcv*pcv*nw) * D/MCPOT` over
  `CC1mu1pi1p_Selected` events in the `w/` MC. Agreement 0.3%.
- **The framework's assembly is correct.** Summing the per-file numuMC CV
  histograms with their own POT scale factors gives 714.2 against an assembled
  `CV_reco` of 716.6. Agreement 0.3%.

What remains is that the thrown data sits above the MC prediction:

| run | thrown | MC | ratio | significance |
|---|---|---|---|---|
| run1 | 252.0 | 254.0 | 0.992 | −0.1σ |
| run2 | 131.0 | 108.1 | 1.212 | **+2.2σ** |
| run4 | 181.0 | 173.1 | 1.046 | +0.6σ |
| run5 | 191.0 | 178.9 | 1.068 | +0.9σ |
| total | 755.0 | 714.1 | 1.057 | +1.5σ |

And that ratio accounts for essentially all of the closure:

    input data/MC ratio      1.057
    closure (identity reg.)  1.067
    left unexplained         1.009

**So the unfolding chain is fine to better than 1%, and the proton-tagged residual
is the fake data itself having fluctuated high** -- by 1.5σ overall in FHC, driven
by run 2 at +2.2σ.

This is a statistical fluctuation in one Poisson throw, not a defect. The
independent checks are FHC (+1.5σ) and RHC (+1.8σ); COMB is their sum and is not a
third independent measurement, so the combined significance is around +2.3σ, not
the ~3σ suggested earlier in this session by treating the three as independent.

The test: re-throw with a different seed. If the closure moves toward unity, this is
confirmed; if it stays at 1.09, something systematic remains. That is cheap --
`throw_perrun_w.C(seed)` takes a seed argument -- and should be done before the
proton-tagged numbers are quoted.

---

## The re-throw test: confirmed statistical

Run 2026-08-30. `throw_perrun_w.C(7)` regenerated the proton-tagged FHC fake data
with a new seed; the `ccpi1p_FHC5_costhmu` universes were rebuilt from it into a
separate output, and the closure re-measured.

**The per-run excesses reshuffled.** A systematic defect would sit in the same run
each time; these did not:

| run | seed 1 | seed 7 | MC | seed1/MC | seed7/MC |
|---|---|---|---|---|---|
| run1 | 252 | 280 | 254.0 | 0.992 | 1.102 |
| run2 | 131 | 109 | 108.1 | **1.212** | 1.008 |
| run4 | 181 | 186 | 173.1 | 1.046 | 1.075 |
| run5 | 191 | 162 | 178.9 | 1.068 | 0.906 |
| total | 755 | 737 | 714.1 | 1.057 | 1.032 |

The run-2 excess that drove the seed-1 result (+2.2 sigma) is 1.008 under seed 7,
and the excess moved to run 1 while run 5 went the other way.

**The closure followed.** For `costhmu` FHC5:

| regularisation | seed 1 | seed 7 |
|---|---|---|
| identity | 1.0615 | 1.0258 |
| second derivative | 1.0615 | 0.9430 |

The input ratio moved 1.057 -> 1.032, predicting a closure near 1.036; identity
returned 1.026. The second derivative went to 0.943, i.e. below unity, which also
shows how much observable-level scatter the derivative regularisation adds on top.

**Conclusion.** The proton-tagged residual is a Poisson fluctuation in the fake-data
throw, not a defect. Combined with the decomposition above -- throw verified to
0.3%, framework assembly verified to 0.3%, under 1% left unexplained after the input
ratio -- the proton-tagged chain is sound.

The production fake data has been restored to the seed-1 throw (755 selected),
which is what the numbers in the note were derived from. The seed-1 backup remains
at `/data/uboone/processed/w/fakedata_seed1_backup/`, and the seed-7 test universes
at `ccpi1p_FHC5seed7_costhmu_univmake.root`.

Worth stating plainly for anyone quoting proton-tagged numbers: a single Poisson
throw carries a ~1.5 sigma normalisation wobble at this sample size, so a closure of
1.06 is not evidence of a bias. Averaging several throws, or quoting the closure
with its throw-to-throw spread, would make that visible rather than leaving it to be
rediscovered.
