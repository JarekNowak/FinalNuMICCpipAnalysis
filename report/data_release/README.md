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

summing over the TRUE index `j`. This reproduces the published generator curves to the
precision of these text files (worst relative per-bin deviation 4e-8 across the released
matrices), so you can verify your implementation against any model column of `curves_*.tsv`
before trusting it on a new model.

Do not apply `A_C` to any column of `curves_*.tsv`: those are smeared already.

## `cov/` — covariance matrices

One directory per extraction (51, named as the `tag` column of `index_extractions.tsv`), each
with 42 files (44 for the inclusive extractions, which also carry the MCS momentum-scale term
`cov_MCSscale.txt` and `cov_total_plusMCS.txt`, identically zero except for `p_mu`): `cov_total.txt`,
`cov_PredTotal.txt`, the per-source terms (`cov_flux*`,
`cov_detVar*`, `cov_xsec_*`, `cov_reint`, `cov_POT`, `cov_numTargets`, `cov_MCstats`,
`cov_EXTstats`, `cov_DataStats`), the blockwise norm/shape/mixed decompositions, `err_prop.txt`,
`unfolding.txt`, `add_smear.txt` and `unfolded_signal.txt`. The covariance is the complete
configured covariance; the identified omissions (MCS momentum for the observables other than
`p_mu`, where it is identically zero; beam-off gate-ratio and Run-2 stand-in terms; beamline-geometry
flux) are listed in the note's Limitations.

Other files: `total_xsec.tsv` (the six one-bin totals with their own covariance breakdown),
`ensemble_*.tsv` (Poisson-ensemble pull widths and offsets, dated), `mcs_scale_*.tsv` (the MCS
momentum-scale term per extraction), `binreview_*.tsv` (the binning review of 2026-09-26),
`index_A_C.tsv` (row sums and conditioning of every released `A_C`), `index_curves.tsv`,
`ext_gates.tsv` (per-run beam-off gate counts and scale factors).

## Coverage

**All 51 differential extractions are released** (release 2026-09-26: fifteen inclusive, five
observables in three configurations, and thirty-six proton-tagged; `index_extractions.tsv`),
each with both its `A_C` matrix and its complete configured covariance decomposition. The six one-bin total
cross sections (inclusive and proton-tagged x FHC/RHC/combined; one true bin over the full
cos(theta_mu) range, Wiener filter off, `A_C` = 1) are in `total_xsec.tsv` with their
uncertainty breakdown and realised fake-data truth. The differential W_pi-p shape is withdrawn
(status column of the index); its files are released for completeness only.

`A_C` depends on the data covariance through the Wiener filter, so it depends on the
central-value weighting fix (commit `51af326`, 2026-08-30); `export_matrices.C` refuses
to write any matrix whose source sidecar predates that fix, so a stale one cannot reach
the release unnoticed.

One matrix has a negative row sum (`1p_RHCFULL_Wpipr`): the smeared content of a bin is a
net-negative combination of the truth. That is permitted for a regularisation operator but
signals an ill-conditioned bin; smear predictions through that row with care.

## Scope

The binnings are those adopted under the 0.50 migration criterion on 2026-09-26 (the note's
Appendix on observable binning): inclusive `p_pi` in **three regions** (`ppi3bin`, open top),
proton-tagged `p_pi` in two (`ppi2bin`); `p_mu` six bins, `cos(theta_mu)` eight; `W_had` two
regions; `delta phi_T` three bins; the proton momentum `pp` three bins. `theta_mu` is no longer
released (it is the same partition as `cos(theta_mu)` with less retained shape). The five-bin
`p_pi` binning is withdrawn and is deliberately not released. Files of the binnings retired on
2026-09-26 are kept in `retired_pre050_20260926/` for the record and are not to be quoted.

Proton-tagged (`1p`) matrices carry the complete configured covariance including detector variations
(first included on 2026-08-31); the differential W_pi-p shape is withdrawn (status
column of `index_extractions.tsv`).

## Regenerating

    cd xsec_analyzer
    root -l -b -q 'macros/export_matrices.C("../report/data_release")'
