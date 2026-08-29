# Results diff: note as printed vs current extraction

Generated 2026-08-30 from `report/current_results.tsv` (51 rows, all re-extracted
2026-08-29/30 after the COMB detector-systematic fix, 80d613d).

## 1. `tab:sigint_all` (analysis_note.tex:3663)

Integrated cross section, $10^{-38}$ cm$^2$/Ar.

| observable | FHC note | FHC now | %  | RHC note | RHC now | %  | COMB note | COMB now | %  |
|---|---|---|---|---|---|---|---|---|---|
| `pmu` | $0.956$ | **0.918** | -4% | $0.717$ | **0.628** | -12% | $0.736$ | **0.780** | +6% |
| `ppi` | $0.965$ | **0.923** | -4% | $0.882$ | **0.850** | -4% | $0.880$ | **0.875** | -1% |
| `costhmu` | $0.966$ | **0.795** | -18% | $0.757$ | **0.711** | -6% | $0.636$ | **0.812** | +28% |
| `costhpi` | $0.829$ | **0.824** | -1% | $0.766$ | **0.779** | +2% | $0.729$ | **0.868** | +19% |
| `thmupi` | $0.987$ | **0.930** | -6% | $0.889$ | **0.878** | -1% | $0.931$ | **0.938** | +1% |
| `thetamu` | $0.861$ | **0.934** | +9% | $0.781$ | **0.917** | +17% | $0.650$ | **0.953** | +47% |
| `thetapi` | $0.878$ | — | — | $0.739$ | — | — | $0.748$ | — | — |

**Every one of the 18 values has moved.** Largest shifts: `costhmu` COMB
$0.636\to0.812$ (+28%), `thetamu` COMB $0.650\to0.953$ (+47%), `costhmu` FHC
$0.966\to0.795$ (−18%).

`theta_pi` is still a row in this table but **is no longer built**. It was dropped
from the driver's observable list (e6d270c): it existed only as the symmetric partner
of `theta_mu`, is not unfolded by stage 3, and cost ~4.7 h per pass. The row must go,
or the observable must be reinstated.

## 2. Claims in the abstract that the current numbers refute

These matter more than the table values, because they are conclusions rather than
inputs. All three are in the abstract (analysis_note.tex:119-140).

### 2.1 "RHC consistent with (or slightly above) FHC" — **refuted**

RHC is below FHC in **all six** inclusive observables:

| observable | FHC | RHC | RHC below FHC by |
|---|---|---|---|
| `pmu` | 0.9179 | 0.6279 | 31.6% |
| `costhmu` | 0.7946 | 0.7113 | 10.5% |
| `ppi` | 0.9231 | 0.8496 | 8.0% |
| `thmupi` | 0.9304 | 0.8782 | 5.6% |
| `costhpi` | 0.8240 | 0.7788 | 5.5% |
| `thetamu` | 0.9344 | 0.9173 | 1.8% |

### 2.2 "an earlier apparent RHC deficit was a per-run POT normalization bug ... now corrected" — **refuted**

The deficit is present in the fake-data **truth** (RHC/FHC $=0.785$) and in the GENIE
tune ($0.779$) — two independent model inputs that agree. It is a property of the
sample, not a normalisation artefact. POT sums exactly
($8.857\times10^{20}+1.1082\times10^{21}=1.9939\times10^{21}$) and the per-configuration
flux values reproduce the POT-weighted combination to 0.01%.

The sentence has the sign of the correction backwards: the measured ratio (0.900)
*understates* the deficit, because the extraction distorts it toward unity.

### 2.3 "self-closure ... recovers the truth to within $\approx6\%$ for every observable" — **stale**

Measured unfolded/truth is **1.13 in all three configurations** — 13% high, not 6%.
This is the single most important open item in the analysis: a uniform Wiener-SVD
normalisation offset on fake data, the same order as the physics differences being
measured, documented but unexplained.

Note that 1.13 being *identical* across configurations is new. Before the detVar fix
the ratios scattered (1.14 / 1.03 / 1.08); the note's own text already asserted the
offset was "identical across configurations", and only now is that true.

## 3. `tab:dagostini` (analysis_note.tex:3462)

Quotes Wiener-SVD means FHC $0.923$ / RHC $0.836$ / COMB $0.852$. Current inclusive
means are **0.887 / 0.794 / 0.871**. The COMB value moved most (the detVar fix), so
the DAg/WSVD ratios and the "20-40% above Wiener-SVD" claim in the abstract both need
restating. **Requires a D'Agostini re-run** — those numbers cannot be recovered from
existing outputs.

## 4. What is confirmed rather than changed

- The observable-to-observable spread within a configuration is real and is the $A_C$
  smearing, as the note says. Current inclusive spreads: FHC 15.8%, RHC 36.5%,
  COMB 17.5%.
- The blind-analysis statement is intact: every "data" point is Poisson-thrown from
  the central-value MC; no beam-on data is unblinded.
- Per-configuration flux normalisation is correct and unchanged.

## 5. Source data

`report/current_results.tsv` — 51 rows: 18 inclusive (6 observables x 3 configs) and
33 proton-tagged (11 x 3). Columns: sigma_int, PredTotal, detVar, flux, xsec
systematics, chi2 vs truth and its p-value. Regenerated 2026-08-29/30 with zero
failures. Sources are `logs/systdump/*.dump`, `logs/systdump/1p_*.raw` and
`/data/uboone/processed/closure_hists_*.root`.
