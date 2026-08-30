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
