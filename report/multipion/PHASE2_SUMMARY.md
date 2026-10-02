# Phase 2 (first part): the particle classifier in the selection

2026-10-01. Phase 2 of `report/MULTIPION_BNB_ADAPTATION_PLAN.md`, item 1, and the question whether the
classifier would also serve the single-pion analysis.

## Integration

- `AnalysisEvent` binds the 21 classifier inputs that nothing else read (per-plane Bragg likelihoods,
  PIDA, chi2 under four hypotheses, calorimetric energy), each only when the branch exists
  (`Branches.hh`).
- `CC1mu1piXp::mp_pid_features()` builds the 43 inputs of a candidate with the conventions of the
  training dumper; `define_constants()` loads the model (`booster_decision_tree/mp_pid/
  mp_pid_rbdt_unmatched_common.root`, or `$CC1MU1PIXP_MP_PID_MODEL`; `TMVA::Experimental::RBDT`,
  linked through `libTMVAUtils`).
- Switches in `CC1mu1piXp.hh`: `eval_new_pid()` computes the four class probabilities per candidate
  (on for N > 1, stored as `pic_pid_mu/pi/p/other`); `use_new_pid()` makes P(π) > `new_pid_cut()` the
  pion identification, with tracks no longer than `new_pid_min_length()` rejected.
- Study selections (`NewPIDSelections.hh`, registered in the factory): `CC1mu1piXpNewPID` (P(π) >
  0.3087 and length > 20 cm, the pion efficiency of the single-pion identification on its candidate
  pool), `CC1mu2piNewPID` and `CC1mu3piNewPID` (P(π) > 0.2276, the efficiency of the split
  `mp_pion_bdt`). Every other cut is the parent's; `<name>_tNN` builds the variant with P(π) > NN/100.
- Fix found on the way: `RBDT` opens and closes the model file, which left ROOT's current directory
  away from the output file; the processed tree then lost its last baskets (29,652 of 30,000 events)
  with 569 write errors. The loading now restores the directory (`TDirectory::TContext`).

## Validation on the 30k-event Run-1 slice

- The released `CC1mu1piXp` and the `CC1mu2pi` outputs are unchanged: all 234 shared branches are
  identical event by event to the Phase 0 reference.
- The class probabilities computed in the selection equal those of the dumper and XGBoost for all
  2,101 matched candidates (largest difference 0), for the first model and for the 43-input one.

## Would the classifier serve the single-pion analysis? Track level

The single-pion pion identification (contained, vertex distance <= 4 cm, length > 20 cm, LLR > 0.1, MIP
and pion BDTs > -0.1, Bragg-pion >= 0.08) against P(π) with the same length cut, on its candidate pool
in events with a muon candidate, held-out runs 3 and 5 (107,426 tracks, 24,919 true pions):

| | Pion efficiency | Purity of the counted tracks | Protons passing | Muons passing | Unmatched passing |
|---|---|---|---|---|---|
| Released single-pion identification | 93.1% | 52.6% | 10.5% | 92.3% | 44.5% |
| P(π), same efficiency | 93.1% | 67.8% | 3.1% | 54.3% | 29.6% |
| P(π), same purity | 99.7% | 55.2% | 8.5% | 91.1% | 47.4% |

Composition of the counted tracks, released against P(π) at the same efficiency: pions 52.6% / 67.8%,
protons 13.5% / 5.2%, muons 14.5% / 11.0%, unmatched 17.7% / 15.1%. The muons in this pool are true
muons of events whose muon candidate is another track. The pion efficiency above 0.175 GeV/c is 93.3%
and 93.1%.

Muon candidate. Choosing, among the tracks the selection already admits as muon candidates, the one with
the highest P(μ) instead of the highest muon-BDT score picks a true muon in 78.0% of the events instead
of 76.6% (a true pion in 8.5% instead of 9.6%; 252,134 events of the held-out runs). The gain is small;
the classifier matters for the pion identification.

## Event-level comparison on the held-out runs

`scripts/mp_pid_eval.py` (full tables in `phase2_eval.md`), FHC Run 5 and RHC Run 3 at their data
exposure, beam-off and dirt included, classifier on the 43 common inputs:

| Selection | Signal | Background (ν + beam-off + dirt) | Efficiency | Purity | Stat. unc. (one bin) |
|---|---|---|---|---|---|
| 1π released | 697.9 | 504.5 | 17.6% | 58.0% | 4.97% |
| 1π with the classifier | 660.4 | 360.0 | 16.6% | 64.7% | 4.84% |
| 2π current | 203.4 | 565.7 | 18.3% | 26.4% | 13.6% |
| 2π with the classifier | 177.7 | 322.6 | 16.0% | 35.5% | 12.6% |
| 3π current | 21.2 | 175.9 | 10.1% | 10.8% | 66% |
| 3π with the classifier | 14.7 | 46.1 | 7.0% | 24.2% | 53% |

For 1π the background falls by 29% for 5% of the signal, the same in both periods (purity 59.8% to
67.0% in FHC Run 5, 57.2% to 63.6% in RHC Run 3); the proton-as-pion background falls from 199 to 78
events, the true-pion backgrounds are unchanged (210 and 214). The released purity on these runs,
58.0%, reproduces the published 57.9%. The statistical gain is small; the larger gain is systematic:
the flux term scales as 17% x (1 + B/S), and B/S falls from 0.72 to 0.55 (2π: 2.78 to 1.82).

The thresholds above were set per track. Scan at event level (SLURM 3450945, `mp_pid_eval.py --scan`,
`phase2_eval_scan.md`), both periods:

| Selection | Pion identification | Efficiency | Purity | B/S | Stat. unc. (one bin) |
|---|---|---|---|---|---|
| 1π | released | 17.57% | 58.0% | 0.72 | 4.97% |
| 1π | P(π) > 0.20 | 17.58% | 63.3% | 0.58 | 4.76% |
| 1π | P(π) > 0.3087 | 16.63% | 64.7% | 0.55 | 4.84% |
| 1π | P(π) > 0.40 | 15.75% | 65.8% | 0.52 | 4.93% |
| 2π | current | 18.25% | 26.4% | 2.78 | 13.6% |
| 2π | P(π) > 0.10 | 19.32% | 28.2% | 2.55 | 12.8% |
| 2π | P(π) > 0.15 | 17.81% | 31.5% | 2.18 | 12.7% |
| 2π | P(π) > 0.2276 | 15.95% | 35.5% | 1.82 | 12.6% |
| 2π | P(π) > 0.30 | 14.00% | 37.9% | 1.64 | 13.0% |
| 3π | current | 10.12% | 10.8% | 8.28 | 66% |
| 3π | P(π) > 0.10 | 8.95% | 16.5% | 5.05 | 57% |
| 3π | P(π) > 0.2276 | 7.02% | 24.2% | 3.13 | 53% |
| 3π | P(π) > 0.30 | 5.75% | 26.7% | 2.75 | 56% |

At P(π) > 0.20 the 1π selection keeps the released efficiency with purity 63.3% and the smallest
statistical uncertainty. The 2π and 3π statistical minima lie at the per-track thresholds, and the
minimum is shallow. Raising the threshold lowers B/S further in all three; the working point is to
be chosen on the total uncertainty (plan, Phase 2 item 4).

## Inputs missing from the older ntuple productions (found 2026-10-01)

The beam-on, beam-off and dirt ntuples of Runs 1–3 (and the Run-4b/c/d and Run-5 beam-off and dirt in
part) lack up to 20 track branches that every overlay carries: the Bragg-pion likelihoods, the
forward-direction flags, truncated dE/dx, deflection, end spacepoints, hit counts and descendants.
The first classifier used them, so its dirt and beam-off yields were artificially zero. It was
retrained on the 43 inputs every sample carries (`--unmatched --common`): on the 2π pool, π-vs-rest
ROC area 0.938 (63 inputs 0.946, current 0.868); at the released 1π efficiency the counted-pion purity
rises from 52.6% to 64.8%. Thresholds 0.3087 (1π) and 0.2276 (2π, 3π). The event-level comparison
above uses this model (SLURM 3450935; first-model outputs kept in `mp_pid_eval_v1_full63`). The
classifier output is cached per event and track for all selections of a job, so threshold variants
cost no extra model evaluations (outputs identical with and without the cache, 379 branches).

The released 1π selection cuts on Bragg-pion >= 0.08 and passes every track where the branch is
absent, so in Runs 1–3 data, beam-off and dirt are selected without the cut while the simulation has
it. The cut removes 3.9% of pions, 8.1% of protons and 4.2% of muons; without it about 5% more events
have exactly one counted pion. The released results (pseudo-data from the overlays) are unaffected;
at unblinding it would bias Runs 1–3. The current `mp_pion_bdt` also takes Bragg-pion as an input
(1.0 where absent). Decision of the user: investigate before fixing.

Investigation, step 1 (`scripts/cr_pid_production_check.py`, `cr_pid_production_check.md`): the per-track
variables both productions carry (LLR score and per plane, proton chi2, track score, length,
MCS/range), beam-on data against the prediction in the opened control regions, per period. They are
mis-modelled in every period (non-muon LLR score shape chi2 160–250/19, mean shifts -0.007 to -0.031),
with no pattern separating the older-production periods (FHC R1, RHC R1, RHC R3) from the newer ones;
sentinel fractions agree between data and prediction everywhere. Not yet checked: the Bragg
likelihoods behind the MIP and pion BDTs, which the processed control-region files do not keep; they
need the raw beam-on files restricted to control-region events (the skim carries no event IDs, so the
restriction has to come from re-running the selection in skim mode with the Bragg branches passed through).

Investigation, step 2 (2026-10-02). `trk_bragg_pion_v` is absent from exactly the older-production
samples: beam-on FHC Run 1, RHC Run 1 and Run 3b, beam-off Runs 1 and 3b, dirt Runs 1 and 3b. Every
overlay, every detector-variation sample and all Run-4 and Run-5 beam-on, beam-off and dirt carry it.
`CC1mu1piXpPIDDiag`, a twin of the released selection with identical decisions (all shared branches
identical on the 30k slice), stores for every pion-pool candidate and for the muon candidate the Bragg
likelihoods, the MIP, pion and muon BDT scores, the LLR and the classifier probabilities. The 8
beam-on files and the matching overlays, beam-off and dirt (34 inputs, SLURM 3450955,
`slurm/slurm_pid_diag_skim.sbatch`) were processed as the control-region skim; every output passed
the skim guards. Comparison: `scripts/cr_pid_diag_check.py`, `cr_pid_diag_check.md`.

- The variables both productions carry are mis-modelled alike in all six periods, with no pattern
  separating the older-production ones: the proton Bragg likelihood of the pion candidates is lower in
  data (mean -0.04 to -0.06, shape chi2 101–710 for 12–16 degrees of freedom) and with it the classifier's P(p); the MIP and
  pion BDT scores are higher (+0.001 to +0.019); the muon candidate's muon Bragg likelihood and P(μ)
  are lower. P(π) is modelled to chi2 23–36/19 on the candidates the identification acts on.
- Pass fractions of those candidates (contained, vertex distance <= 4 cm, length > 20 cm), data over
  prediction:

| | FHC R1 (old) | FHC R4 | FHC R5 | RHC R1 (old) | RHC R3 (old) | RHC R4 |
|---|---|---|---|---|---|---|
| Released identification | 1.23 ± 0.05 | 1.04 ± 0.04 | 1.06 ± 0.04 | 1.04 ± 0.08 | 1.15 ± 0.03 | 0.99 ± 0.03 |
| The same without Bragg-pion | 1.11 ± 0.04 | 1.04 ± 0.04 | 1.10 ± 0.04 | 0.97 ± 0.07 | 1.06 ± 0.03 | 0.98 ± 0.03 |
| P(π) > 0.3087 | 1.14 ± 0.05 | 1.04 ± 0.05 | 1.08 ± 0.05 | 0.91 ± 0.08 | 1.06 ± 0.03 | 1.05 ± 0.04 |

  In the older-production periods the released identification passes 7–11% more data candidates than
  predicted on top of the ratio without the cut: the data skip the cut, the overlay applies it and
  removes 7–11% of the candidates that pass the other three cuts. Where the data carry the branch, the
  cut after the other three passes 90.1, 86.4 and 91.6% of the data candidates against 89.8, 89.6 and
  91.0% predicted (FHC R4, R5, RHC R4). Without the cut, or with P(π), the older and newer periods agree.

The older production is therefore not otherwise different, and the fail-open cut shows up in the
opened control regions at the size the simulation predicts.

Event level (`CC1mu1piXpNoBragg`, `CC1mu1pi1pNoBragg`: the released selections without the cut in every
sample; release samples at data exposure, SLURM 3450989, `scripts/bragg_cut_eval.py`,
`bragg_cut_eval.md`; on the 30k slice the extra counted pions are exactly the candidates that fail only
the cut):

| COMB, released / without the cut | Efficiency | Purity | B/S | Stat. unc. (one bin) | Bias at unblinding, cut as released |
|---|---|---|---|---|---|
| Inclusive | 17.55 / 18.39% | 57.9 / 56.9% | 0.73 / 0.76 | 3.06 / 3.01% | +4.7% |
| Proton-tagged | 11.45 / 12.06% | 50.7 / 50.4% | 0.97 / 0.98 | 5.23 / 5.10% | +4.6% |

The last column is the overlay change in the older-production periods over the released signal: the
excess the one-bin cross section would carry if the older-production data are selected without the cut
and the overlay with it (FHC +3.4%, RHC +5.7%; +5.5% and +8.4% if Run-2 data come from the older
production). Per period the data would exceed the prediction by 4.3–7.3%. The released results
(pseudo-data from the overlays) are unaffected. Dropping the cut everywhere costs one point of purity
in the inclusive selection and 0.3 in the proton-tagged one, raises B/S by 0.03 and 0.01, and lowers
the statistical uncertainty slightly. Decision of the user (2026-10-02): drop the cut in every sample of
the released single-pion selections, with a new release before unblinding.
