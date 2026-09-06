# Numerical data release — additional smearing matrices

## What these are

The Wiener-SVD unfolding does not estimate the true distribution. It estimates an
*additionally smeared* one,

    x_hat  ~  A_C . x_true

so a prediction must be transformed by the same `A_C` before it can be compared with the
published values. Comparing an unsmeared prediction to these results is a category error.

## Start here: `curves_*.tsv`

**If you want to compare against the measurement, use these.** Each file holds, per bin,
the unfolded data with its uncertainty and every model curve shown in the note — all of
them *already smeared by `A_C`* and therefore directly comparable bin by bin. No
transformation is required on your side.

## `A_C_*.tsv` — for smearing your own model

For a prediction supplied as a per-bin cross section `p[j]` in true bin `j`:

    smeared_dsigma_dx[i] = ( sum_j A_C[i][j] * p[j] ) / width[i]

summing over the TRUE index `j`. This reproduces the published generator curves **exactly**
— worst per-bin deviation 1.2e-10 across the released matrices — so you can verify your
implementation against any model column of `curves_*.tsv` before trusting it on a new model.

Do not apply `A_C` to any column of `curves_*.tsv`: those are smeared already.

## `cov/` — covariance matrices

## Coverage

**All 51 differential extractions are released** (release 2026-09-06; `index_extractions.tsv`),
each with both its `A_C` matrix and its full covariance decomposition. The six one-bin total
cross sections (inclusive and proton-tagged x FHC/RHC/combined; one true bin over the full
cos(theta_mu) range, Wiener filter off, `A_C` = 1) are in `total_xsec.tsv` with their
uncertainty breakdown and realised fake-data truth. The differential W_pi-p shape is withdrawn
(status column of the index); its files are released for completeness only.

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

Proton-tagged (`1p`) matrices carry the full covariance including detector variations
(all 51 extractions, rebuilt 2026-08-31); the differential W_pi-p shape is withdrawn (status
column of `index_extractions.tsv`).

## Regenerating

    cd xsec_analyzer
    root -l -b -q 'macros/export_matrices.C("../report/data_release")'
