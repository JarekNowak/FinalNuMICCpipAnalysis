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

Use these when you have a prediction that is not one of the released curves.

    x_hat[i] = sum_j A_C[i][j] * p[j]      (sum over the TRUE index j)

**A caveat you must read before relying on this.** The released `A_C` is bit-identical to
the operator the analysis code applies internally (verified to 5e-11 against the
full-space matrix, and the internal identity `U.R = A_C` holds to 4.4e-6). What is *not*
currently demonstrated is the full path from a raw generator file to a published curve:
between the two sit a units conversion applied when a prediction is loaded and a
`1/(conversion_factor * width)` transform applied when it is plotted, and reproducing the
published curves from a raw generator file plus this matrix alone has not been shown to
work. Applying `A_C` to a raw generator histogram does **not** reproduce the plotted
curve for that generator.

So: `A_C` is correct as an operator, and the smeared curves are correct as results, but
the release does not yet let you rebuild the second from a raw model file. Until it does,
compare against `curves_*.tsv`, and treat `A_C` as the right operator for a prediction you
have already put in the analysis's own true-bin event units.

Do not apply `A_C` to `truth_smeared` or to any model column in `curves_*.tsv` — those are
smeared already, and smearing them again gives roughly 0.77x the correct answer.

## `cov/` — covariance matrices

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
