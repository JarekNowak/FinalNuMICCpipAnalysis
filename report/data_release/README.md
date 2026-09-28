# Numerical data release

MicroBooNE NuMI charged-current single-charged-pion cross sections on argon, inclusive
(CC1mu1piXp) and proton-tagged (CC1mu1pi1p), in the FHC, RHC and combined configurations.
**All values are extracted from pseudo-data** (Poisson throws of the central-value simulation);
the signal region is blind.

## Conventions

- **Flux.** Cross sections are averaged over the nu_mu + nubar_mu flux of the new-G4 PPFX central
  value **integrated above 60 MeV**: FHC 6.81159, RHC 6.44646, combined 6.59691 x 1e-10
  cm^-2 POT^-1. The threshold includes the flux below the pion-production threshold, which lowers
  the flux-averaged value relative to a spectrum restricted to the pion-production region.
- **Exposure.** FHC 7.766e20, RHC 11.082e20, combined 18.848e20 POT.
- **Target.** 8.710e29 argon nuclei in the fiducial volume; cross sections are per argon nucleus,
  in 1e-38 cm^2/Ar (per unit of the variable for differential results).
- **Angles** are measured about the neutrino direction (the NuMI target-to-vertex direction at
  reconstruction level), not the detector z axis.

## Status of the results

Every extraction has one of four statuses (`status` column of `index_extractions.tsv`; product
level in `../status.tsv`, the single source used by the analysis documents):

- **Primary result**: inclusive p_mu, p_pi (three regions), cos(theta_mu), cos(theta_pi),
  theta_mu_pi and the inclusive one-bin total, FHC/RHC/combined.
- **Secondary result**: proton-tagged one-bin total; proton-tagged theta_p
  and theta_pi_p (FHC/RHC/combined); the two-dimensional theta_p x delta p_T (proton-tagged) and
  cos(theta_pi) x cos(theta_mu) (inclusive), combined configuration only.
- **Validation-only product**: every other proton-tagged extraction.
- **Not reported**: the differential W_pi_p (its files are kept for completeness; no differential
  W_pi_p result is reported, and the proton-tagged one-bin total is not a W_pi_p measurement).

## What the unfolded values are

The Wiener-SVD unfolding estimates an *additionally smeared* distribution,

    x_hat  ~  A_C . x_true

so a prediction must be transformed by the same `A_C` before it can be compared with the released
values. Comparing an unsmeared prediction with these results is a category error. The one-bin totals
have A_C = 1.

## Files

| File | Content |
|---|---|
| `curves_<tag>.tsv` | Per bin: `bin`, `low`, `high`, `width`, `unfolded_data`, `stat_err`, then `truth_smeared`, `tune_smeared` and the four generator columns, **already smeared by A_C**. Densities: multiply by width to integrate. Two-dimensional tags have `inner_low`, `inner_high`, `outer_low`, `outer_high`, `area` instead; analysis bins are ordered outer slice by outer slice, inner variable fastest. |
| `A_C_<tag>.tsv` | Additional-smearing matrix; row = smeared bin i, column = true bin j. |
| `cov/<tag>/cov_total_plusMCS.txt` | **The official covariance**: framework covariance plus the data-side MCS muon-momentum-scale term. Use this one. |
| `cov/<tag>/cov_total.txt` | Framework covariance alone. |
| `cov/<tag>/cov_MCSscale.txt` | The MCS term (identically zero for observables that do not use the muon momentum). |
| `cov/<tag>/cov_<source>.txt` | Per-source terms (`cov_flux*`, `cov_detVar*`, `cov_xsec_*`, `cov_reint`, `cov_POT`, `cov_numTargets`, `cov_MCstats`, `cov_EXTstats`, `cov_DataStats`), blockwise norm/shape/mixed decompositions. |
| `cov/<tag>/unfolded_signal.txt` | Unfolded result as bin integrals (the units of the covariances). |
| `total_xsec.tsv` | The six one-bin totals (inclusive and proton-tagged x FHC/RHC/combined) with their uncertainty breakdown and realised pseudo-data truth. |
| `index_extractions.tsv` | Every released differential extraction with its status, bins, integral, uncertainty and closure. |
| `mcs_scale_*.tsv`, `ensemble_*.tsv`, `binreview_*.tsv`, `ext_gates.tsv`, `index_A_C.tsv`, `index_curves.tsv` | Supporting summaries. |

Tags are `<family>_<CFG>_<observable>` with family `incl` or `1p` and CFG in `FHC5`, `RHCFULL`,
`COMB`; the two-dimensional tags are `1p_COMB_thetap_dpt2d` and `incl_COMB_costhpi_costhmu2d`.

## Comparing a prediction

**Start with `curves_*.tsv`**: every model shown in the notes is there, already smeared, and directly
comparable with `unfolded_data` bin by bin using `cov_total_plusMCS.txt` (bin-integral units: multiply
the density difference by the width or area).

For a new prediction given as a cross section `p[j]` integrated over true bin `j`:

    smeared_dsigma_dx[i] = ( sum_j A_C[i][j] * p[j] ) / width[i]

summing over the TRUE index `j`. This reproduces the published generator curves to the precision of
these text files (worst relative per-bin deviation 4e-8), so an implementation can be checked against
any model column of `curves_*.tsv` before it is applied to a new model. Do not apply `A_C` to any column
of `curves_*.tsv`: those are smeared already.

## Terms outside the official covariance

The beamline-geometry flux term (1.6-3.3%), the beam-off gate-ratio term and the Run-2 beam-off
stand-in term are bounded in the analysis note but not yet in the covariance (prerequisite U7).

## Regenerating

    cd xsec_analyzer
    root -l -b -q 'macros/export_curves.C("../report/data_release")'
    root -l -b -q 'macros/export_matrices.C("../report/data_release")'
    python3 ../report/tools/export_2d.py
    python3 ../report/tools/mcs_eval.py all incl; python3 ../report/tools/mcs_eval.py all 1p
    python3 ../report/tools/index_extractions.py

`export_matrices.C` refuses any matrix whose source predates the central-value weighting fix
(commit `51af326`, 2026-08-30).
