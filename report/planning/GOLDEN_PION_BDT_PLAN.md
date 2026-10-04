# Plan: a golden-pion BDT for the $p_\pi$ measurement

Written 2026-08-17. Every number below is measured on Run-1 FHC, not assumed — the
feasibility checks are already done and recorded in the session; what remains is execution.

## Why

`candidate_pion_mom_reco` is a range momentum that **saturates near 0.34 GeV/c** because
charged pions interact hadronically before stopping. The consequence, measured on the five
analysis bins, is a nearly lower-triangular response:

| true bin | diagonal | → reco bin 1 |
|---|---|---|
| 0.175–0.250 | **94.3%** | 94.3% |
| 0.250–0.320 | 36.5% | 59.0% |
| 0.320–0.420 | 23.1% | 49.0% |
| 0.420–0.550 | 10.4% | 43.3% |
| > 0.550 | **6.2%** | 32.6% |

The BNB CC1π note requires a smearing diagonal **> 0.68** for a bin to be used. One of our
five passes. Their fix is a *golden pion* — one that ranges out rather than interacting —
selected by a dedicated BDT applied to the $p_\pi$ cross section only.

## What is already established

1. **Wiggliness exists.** `trk_avg_deflection_stdev_v` (plus `_mean_v`,
   `_separation_mean_v`) is in the PeLEE ntuples and fully populated. Now passed through
   by `ProcessNTuples` (commit `a82fe82`), verified populated on 100% of selected events
   with the selection bit-identical.
2. **It discriminates, with the opposite sign to the naive reading.** Golden pions are
   ~2× *wigglier* (a stopping pion scatters hard through its final slow centimetres; an
   interacting one is truncated while fast and stays straight). Survives binning in track
   length, so it is not a momentum confound. Single cut: golden fraction 39.9% → 53.0% at
   86.8% efficiency.
3. **A clean truth label exists.** `mc_end_p < 0.01 GeV/c` = the pion ranged out.
4. **The truth join is exact.** `backtracked_px/py/pz` → `mc_px/py/pz` matches **100%** of
   candidate pion tracks with **0%** ambiguity.
5. **The label targets the right thing** — it is not circular, being defined on truth alone:

   | | ⟨p_corr/p_true⟩ | within 20% |
   |---|---|---|
   | golden | 0.770 | **58.3%** |
   | non-golden | 0.489 | **7.9%** |

6. **Statistics are ample**: ~43k labelled pion tracks in Run-1 FHC, ~173k across all FHC.

## Steps

### 1. Build the training sample
New macro `macros/golden_pion_train.C`, modelled on the existing
`booster_decision_tree/mp_pion_bdt/train_bdt.C` (same TMVA scaffolding, so this is a
variation on working code rather than something new).

- Loop the raw PeLEE ntuple; for each track-like PFP with `backtracked_pdg == ±211`,
  purity > 0.5, `trk_score > 0.5`, `trk_len > 5` cm.
- Join to truth via the momentum match; label `golden = mc_end_p < 0.01`.
- Restrict to the analysis phase space (`p_true > 0.175 GeV/c`) so the classifier is
  trained where it will be applied.
- Split 50/50 train/test, and hold out **entire runs** rather than random events, so the
  quoted performance is not inflated by run-correlated detector conditions.

### 2. Inputs
Start from the BNB's four and add what we have:

| variable | status |
|---|---|
| `trk_avg_deflection_stdev_v` | verified discriminating |
| `trk_avg_deflection_mean_v`, `_separation_mean_v` | available |
| `trk_bragg_pion / _mip / _mu / _p` | already used by the multi-pion BDT |
| `trk_score_v` | available |
| `pfp_trk_daughters_v`, `pfp_shr_daughters_v` | available |
| `trk_len_v` | **include with care** — see below |

**Track length is the one to watch.** It correlates directly with momentum, so a BDT that
leans on it produces a strongly momentum-dependent efficiency. That is not fatal (the
response matrix carries the efficiency) but it inflates model dependence in exactly the
observable we are trying to protect. Train twice, with and without length, and compare the
efficiency-vs-$p_\pi^{\rm true}$ curves; prefer the flatter one unless the gain is large.

### 3. Choose the working point on the right metric
Not ROC area. The metric that matters is **the migration-matrix diagonal**, because that is
what the BNB criterion is stated in and what governs whether a bin is measurable.

For each candidate cut, rebuild the true→reco migration on the corrected estimator and read
off the five diagonals. Target: as many bins as possible above 0.68, at the least efficiency
cost. The BNB accepted 20.3% → 9.5% efficiency for 36.3% → 67.7% golden; that is the
precedent for how much efficiency is worth trading.

### 4. Rebin $p_\pi$
Once the diagonals are known under the cut, re-derive the binning against both BNB criteria:
≥100 predicted MC events per bin **and** diagonal > 0.68. Expect something closer to their
shape — resolved bins only where the estimator tracks truth, then a wide bin and an
overflow — rather than our current five edges spread through the saturated region.

### 5. Wire it in as a subset, not a global cut
Follow the BNB structure exactly: the golden cut applies to the **$p_\pi$ differential cross
section only**. The total, angular and W/TKI results stay on the generic selection. Mechanically
that means a new bin config plus its own univmake and unfold, in the same way the analysis
already carries separate configs per observable.

### 6. Validate
- **Closure**: unfold the fake data with the golden subset and confirm it recovers truth.
- **Blind-safe data/MC**: compare the BDT response distribution in data and MC. This is a
  reco-level distribution, so it can be checked before unblinding — and it is the check the
  BNB note shows (their Figure 41).
- **Efficiency shape**: efficiency vs true $p_\pi$, flat-ish preferred.
- **Detector systematics**: rerun the detVar samples through the cut. A BDT trained on CV
  simulation can respond differently under detector variations, and this cut sits upstream
  of the measurement, so its detector systematic must be evaluated, not assumed small.

## A parallel fix worth doing regardless

Even *golden* pions reconstruct at ⟨p_corr/p_true⟩ = **0.770** — still 23% low. The
kinetic-energy-preserving μ→π mass swap in the bin config is an approximation: the correct
transformation holds the **range** fixed, not the kinetic energy, and inverting the CSDA
range–energy relation properly is a small, self-contained change. It would improve every
pion momentum, golden or not, and is independent of the BDT. Worth doing first, since it may
shift the diagonals on its own and change what binning the BDT then has to achieve.

## Effort and risk

| step | effort | risk |
|---|---|---|
| training sample + labels | low — join and label verified | low |
| BDT training | low — existing TMVA scaffolding | low |
| working point + rebinning | medium | low |
| univmake + unfold for the subset | ~2 h compute per config | low |
| detector systematics for the cut | medium | **the real one** |
| exact range→momentum conversion | low | low, and helps regardless |

Nothing here changes the total, angular or W/TKI results. The honest summary is that this
makes $p_\pi$ a measurement rather than something largely inferred from the response model
above ~0.5 GeV/c — and that the parallel mass-conversion fix should be tried first because
it is cheaper and may move the diagonals on its own.
