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

**All 51 extractions are released**, each with both its `A_C` matrix and its full
covariance decomposition.

`A_C` depends on the data covariance through the Wiener filter, so it depends on the
central-value weighting fix (commit `51af326`, 2026-08-30); `export_matrices.C` refuses
to write any matrix whose source sidecar predates that fix, so a stale one cannot reach
the release unnoticed.

Two matrices have negative row sums (`1p_RHCFULL_Wpipr`, `1p_RHCFULL_costhmu`): the
smeared content of a bin is a net-negative combination of the truth. That is permitted
for a regularisation operator but signals an ill-conditioned bin; smear predictions
through those rows with care.

## Scope

`p_pi` matrices are the adopted **two-bin** scheme (`ppi2bin`).  The five-bin `p_pi`
binning is withdrawn and is deliberately not released.

Proton-tagged (`1p`) matrices correspond to results whose uncertainties do **not** yet
include detector variations; see the analysis note's status table before using them.

## Regenerating

    cd xsec_analyzer
    root -l -b -q 'macros/export_matrices.C("../report/data_release")'
