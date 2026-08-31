# Numerical data release — additional smearing matrices

## What these are

The Wiener-SVD unfolding in this analysis does not estimate the true distribution.  It
estimates an *additionally smeared* one,

    x_hat  ~  A_C . x_true

so a theory prediction cannot be compared with the published central values directly.
It must first be multiplied by the same `A_C`.  Comparing an unsmeared prediction to
these results is a category error and will produce a meaningless chi-square.

That makes `A_C` part of the measurement rather than an internal detail of it, which is
why it is released here.

## Files

- `A_C_<family>_<config>_<observable>.tsv` — one matrix per extraction.
  - `family` is `incl` (inclusive) or `1p` (proton-tagged).
  - `config` is `FHC5`, `RHCFULL` or `COMB`.
  - Row = smeared bin *i*, column = true bin *j*, value = `A_C[i][j]`.
  - Apply as `x_hat[i] = sum_j A_C[i][j] * p[j]`, summing over the **true** index.
- `index_A_C.tsv` — dimensions, row-sum range, and the source ROOT file for each matrix.

## Read the row sums before using these

`A_C` is **not** norm preserving, and the departure is large and observable-dependent.
The row-sum range in `index_A_C.tsv` quantifies it: it runs from `0.046`–`0.412` for
inclusive RHC `cos(theta_mu)` up to `0.334`–`1.269` for proton-tagged RHC `W_pi p`.

Two consequences:

1. **Do not integrate a differential result here and call it a total cross section.**
   Summing the Wiener-SVD differential result over bins gives an observable-dependent
   number, because each observable is suppressed differently by its own `A_C`.  The
   physical total is quoted separately in the analysis note from the cut-and-count
   extraction, which is not subject to this.
2. The spread between observables in the published integrated values is this effect, not
   a physical or acceptance difference.  The D'Agostini cross-check, which is flat across
   observables, is what establishes that.

## Coverage

The Wiener filter depends on the data covariance, so `A_C` depends on the central-value
weighting fix (commit `51af326`, 2026-08-30).  A matrix computed before that fix is not
the released measurement.

**45 of 48 extractions are released.**  The three withheld are the inclusive `p_pi`
matrices, whose sidecars predate the fix; they are pending re-extraction.
`export_matrices.C` refuses to write a matrix whose source predates the fix, so they will
appear automatically once re-run.

The proton-tagged matrices were regenerated on 2026-08-31 from universes rebuilt with the
detector variations included, and postdate the fix.

Two matrices have negative row sums (`1p_RHCFULL_Wpipr`, `1p_RHCFULL_costhmu`): the
smeared content of a bin is a net-negative combination of the truth.  That is permitted
for a regularisation operator but signals an ill-conditioned bin; smear predictions
through those rows with care.


## Covariance matrices

`cov/<family>_<config>_<observable>/` holds, for each extraction:

- `cov_<component>.txt` — covariance in cross-section units, one per systematic
  component (`flux`, `reint`, `xsec_total`, each `detVar*`, `POT`, `numTargets`, the
  statistical terms) plus the aggregates `PredTotal` and `total`, and the blockwise
  shape/norm/mixed decompositions.
- `unfolded_signal.txt` — the central values these covariances belong to, as
  `sigma` per bin (multiply-by-width already applied; summing the file gives the
  integrated cross section).
- `add_smear.txt`, `unfolding.txt`, `err_prop.txt` — the A_C, unfolding and
  error-propagation matrices for the same extraction.

Matrix files are `xbin ybin value`, with `numXbins`/`numYbins` on the first two lines.

**Verified**: for all 33 proton-tagged extractions, the mean per-bin fractional
uncertainty derived from `cov_PredTotal.txt` and `cov_detVar_total.txt` reproduces the
`PredTotal_pct` and `detVar_pct` columns of `report/current_results.tsv` to better than
0.1 percentage points.

These were produced by `dump_overall_results()` in `UnfolderNuMI.C`, which had been
written but never called — no covariance matrices were being produced at all. Any
`unfold_output/*.txt` file older than 2026-08-31 comes from a different build and does
not correspond to these results.

## Scope

`p_pi` matrices are the adopted **two-bin** scheme (`ppi2bin`).  The five-bin `p_pi`
binning is withdrawn and is deliberately not released.

Proton-tagged (`1p`) matrices correspond to results whose uncertainties do **not** yet
include detector variations; see the analysis note's status table before using them.

## Regenerating

    cd xsec_analyzer
    root -l -b -q 'macros/export_matrices.C("../report/data_release")'
