# Phase 2 (first part): the particle classifier in the selection

2026-10-01. Phase 2 of `report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md`, item 1, and the question whether the
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

## Track assignment and event classifier (2026-10-06)

Phase 2, items 1–4, on the branch merged with `main` at version 1.7 (the merge sets the Bragg-pion
default to off; on the 30k-event slice the outputs of `CC1mu1piXp` and `CC1mu1pi1p` are bitwise
identical to version 1.7, 229 branches). Script `scripts/mp_evt.py` (steps `build`, `build_detvar`,
`train`, `syst`, `detsyst`, `wp`); full tables in `phase2_event.md`, `phase2_wp.md` and their `_3pi`
versions.

### Event-level track dump

`CC1mu2piEvt` (`NewPIDSelections.hh`) is `CC1mu2pi` plus, in events with one neutrino slice and the
vertex in the fiducial volume, every track of the pion pool with track score >= 0.5, the muon candidate
included (`etrk_*`: the four class probabilities, length, vertex distance, containment, range and MCS
momentum, direction, backtracked truth). Tracks shorter than the 5 cm of the classifier's training
sample are kept (1,262 of 12,932 on the slice). On the slice its 86 shared branches equal those of
`CC1mu2pi`, and the probabilities equal those of the counted candidates (1,965 tracks, largest
difference 0). Processed: the 28 release samples (SLURM 3520138, `slurm/slurm_mp_evt.sbatch`,
`/data/uboone/processed/mp_evt`) and the 18 Run-4 detector-variation samples (SLURM 3532497,
`slurm/mp_evt_detvar_manifest.list`). The normalisation reproduces the comparison above: `CC1mu2pi`
203.4 signal and 565.8 background events on the held-out periods, `CC1mu2piNewPID` 177.7 and 322.6.

### Assignment (item 1)

For "one muon, N pions, the rest protons or other" the log-likelihood is
L_N = Σ_rest log(1 − P_μ − P_π) + log P_μ(muon) + Σ_pions log P_π. Every choice of the muon and the N
pions is enumerated (at most 8 tracks) and the best and the runner-up are kept for N = 1, 2, 3; a
brute-force check on 200 events agrees exactly. On the slice, 41 of the 136 preselected two-pion signal
events (30%) have a true muon and both true pions among the tracks, and in 67 only one pion is a track,
so in most signal events no assignment can be fully right. Keeping the events whose most likely pion
count is two gives 21.9% efficiency at 32.5% purity (COMB).

### Event classifier (items 2 and 3)

XGBoost on 33 features (PID set: probabilities of the assigned tracks, L_2 − L_1, L_2 − L_3, the margin
over the runner-up, class sums and counts over all tracks, track and shower counts, topological score,
CosmicIP) or 48 (full set: adding lengths, vertex distances, containment, track scores and opening angles
of the assigned tracks). Domain: one slice, vertex in the fiducial volume, software trigger, at least three
tracks, and θ(μ, longest assigned pion) < 2.6 rad (D6). Trained on run periods 1, 2 and 4 of both horn
modes (206,152 events, 13,558 signal; the Run-3b dirt, shared with a test period, left out), tested on
FHC Run 5 and RHC Run 3, which the particle classifier did not see either (104,814 events, 7,253 signal).
ROC area on the test periods 0.887 (PID set) and 0.898 (full set). The largest gains: the sum of P(π) over
the tracks, L_2, the number of tracks with P(π) > 0.5, the topological score, the number with P(μ) > 0.5.
At the efficiency of `CC1mu2pi` the purity doubles (score > 0.90: 19.9% and 51%, against 18.2% and 26%).

### Working point (item 4)

Expected total uncertainty of the one-bin total, with the test-period yields projected to the full
exposure: data statistics; flux 17% × (1 + B/S); POT and targets (2% and 1%) × (1 + B/S); cross section
as `configs/ccpi_systcalc_numi.conf` (GENIE multisim, RPA and SCC multisims, eight GENIE unisims);
reinteraction; detector (the eight Run-4 variations, FHC and RHC correlated); MC statistics. The current
selections carry the same terms. COMB:

| Selection | Efficiency | Purity | Stat. | Flux | POT, targets | Cross section | Reint. | Detector | MC stat. | Total |
|---|---|---|---|---|---|---|---|---|---|---|
| `CC1mu2pi` | 18.2% | 26% | 8% | 65% | 9% | 47% | 17% | 34% | 5% | 89% |
| `CC1mu2piNewPID` | 15.8% | 35% | 8% | 48% | 6% | 27% | 13% | 21% | 4% | 61% |
| Classifier, PID set, score > 0.92 | 17.1% | 56% | 6% | 31% | 4% | 12% | 5% | 33% | 3% | 47% |
| Classifier, full set, score > 0.94 | 14.4% | 62% | 6% | 27% | 4% | 10% | 5% | 33% | 3% | 45% |
| Classifier, PID set, score > 0.9725 | 3.9% | 76% | 11% | 22% | 3% | 6% | 4% | 78% | 6% | 83% |

- Statistics and flux alone put the minimum at score > 0.9725 with 3.9% efficiency. The detector term
  rises steeply above 0.95 (43% at 0.95, 78% at 0.9725), so the minimum with all terms is at 0.92 (PID
  set), and the total is flat between 0.90 and 0.94 to within the statistics of the detector samples
  (±3%).
- At 0.92 the detector term comes from Recomb2 (22%), WMAngleYZ (18%), WMX (12%) and WMYZ (8%), the
  variations that change the calorimetry behind the PID inputs; light yield and SCE stay below 5%.
  `CC1mu2pi` shows the same pattern (Recomb2 23%, WMAngleYZ −19%).
- The cross-section term is the background model (GENIE multisim 47% for `CC1mu2pi`, 12% at 0.92) and
  falls with the purity; the model dependence of the efficiency stays below 6% even at 4% efficiency.
- The full set gains 2% in total, within the noise of the detector term; the PID set does not use the
  kinematics of the observables.

At score > 0.92 (PID set): 499 signal and 398 background events at full exposure (FHC 200 and 157, RHC 298
and 241). Background: CC with three or more charged pions 116, CC1π 96, CCπ0 73, CC2π outside the signal
26, NC 25, CC other 22, outside the fiducial volume 22, dirt 7, CC0π 6, beam-off 4. In selected signal the
assigned muon and pions are correct in 91–92% of events; in the neutrino background the assigned pions are
protons in 12–14%. The efficiency against the true pion momenta follows that of `CC1mu2pi` except in the
lowest leading-pion bin (0.10–0.175 GeV/c: 6.9% against 10.1%). 55% of the selected signal is shared with
`CC1mu2pi`. Of all selected events, 14% lie in the blind single-pion signal region (D7) and 12% in the
multi-π control region of the single-pion analysis, which has been opened (D3); for `CC1mu2pi` the
fractions are 16% and 12%.

### Three or more pions (D8)

Events of the domain with at least four tracks in which three pions are more likely than two: 234 at full
exposure (COMB), of which CC with three or more pions 88, CCπ0 41, two-pion signal 35, CC other 13,
beam-off 13; 6% of them pass the working point. A three-pion classifier built the same way reaches a ROC
area of 0.855 (PID set), but with all terms every three-pion selection has a total above 100% (`CC1mu3pi`
245%, `CC1mu3piNewPID` 117%, the classifier at best 105%).

### Open

- Working point and feature set for Phase 3 (proposed: PID set, score > 0.92).
- The detector term (33%) is the largest. The levers are a classifier less sensitive to the calorimetric
  inputs, or a constraint from the control samples (Phase 4).
- The 12% overlap with the opened multi-π region (D3).

## Detector robustness (2026-10-07)

The detector term leads the two-pion working point. Inputs: the 27 Run-4 detector-variation files dumped track by track
(SLURM 3533973, `/data/uboone/processed/mp_pid_detvar`), tracks matched across the samples by event and backtracked particle
(`scripts/mp_pid_detsens.py`); particle-classifier variants (`scripts/mp_pid_train.py --tag --drop --aug`,
`slurm/slurm_mp_pid_train.sbatch`); two of them deployed and the release and detector-variation samples reprocessed with each
(SLURM 3538383–3538386, `/data/uboone/processed/mp_evt_{nochi,llronly}`, `CC1MU1PIXP_MP_PID_MODEL`), then the event chain of
`scripts/mp_evt.py` on each (`MP_EVT_DIR`, `--tag`).

### Statistics of the detector-variation samples

Paired bootstrap (`mp_evt.py detboot`, 200 replicas, one Poisson weight per physical event shared by the CV and the eight
variations; `phase2_detboot.md`): at score > 0.92 the statistics contribute 7.7% of the 32.7% (31.8% without them). Recomb2
+22.1 ± 3.6%, WMAngleYZ +17.9 ± 3.2%, WMX +12.2 ± 2.6%, WMYZ +8.4 ± 2.7%; WMAngleXZ, SCE and the light-yield variations are
consistent with zero. `CC1mu2pi`: 33.7% (11.2% from statistics), `CC1mu2piNewPID`: 20.6% (9.9%).

### Where the sensitivity sits

- Inputs (`detsens_inputs.md`): χ², PIDA, the per-plane Bragg likelihoods and the calorimetric energy move by up to 0.5 of
  their spread (protons, Recomb2 and WMAngleYZ; at most 0.16 for pions); the LLR score moves by at most 0.04 and the geometric
  inputs not at all.
- Deployed classifier at 80% pion efficiency: the pion efficiency moves by at most 2.4%; muons passing as pions by up to
  +24% and protons by up to −27% (Recomb2).
- Event level, score > 0.92: Recomb2 and WMAngleYZ lower the selected signal by about 8% and the background by 10–18%, and
  both raise the extracted total. The detector term is an efficiency part (14.3%) and a background part (17.8%) that add in
  each variation.

### What was tried (COMB)

| Particle classifier | Event features | Score > | Efficiency | Purity | Detector (efficiency part, background part) | Total |
|---|---|---|---|---|---|---|
| Deployed | PID set | 0.92 | 17.1% | 56% | 32.7% (14.3%, 17.8%) | 47% |
| Deployed | full set | 0.94 | 14.4% | 62% | 33.0% (18.2%, 14.8%) | 45% |
| Without χ², PIDA, calorimetric energy | PID set | 0.93 | 14.4% | 57% | 35.1% (14.9%, 19.0%) | 48% |
| Without χ², PIDA, calorimetric energy | full set | 0.94 | 13.4% | 61% | 29.1% (14.3%, 15.2%) | 43% |
| LLR and geometry only | PID set | 0.92 | 13.0% | 51% | 21.4% (7.8%, 15.2%) | 44% |
| LLR and geometry only | full set | 0.93 | 12.0% | 55% | 20.3% (8.4%, 12.9%) | 42% |

At track level the LLR-and-geometry classifier has a pion ROC area of 0.892 against 0.926, and lets 9.0% of protons pass at
80% pion efficiency against 3.9%. Also tried, without gain: the event classifier without the topological score, CosmicIP and
the track and shower counts (detector 54%, total 70%), and the event classifier trained with the detector-variation events
(even event numbers; on the odd ones its best total is 44% against 47% for the nominal classifier, the PID set 51–54% against
53%). The particle classifier trained with the detector-variation tracks changes nothing at track level.

- No variant removes the detector term. The totals of all variants lie between 42% and 48% at their best thresholds, within
  the 3% scatter of the scans.
- The LLR-and-geometry classifier halves the efficiency part (14.3% to 7.8%), the part a constraint of the background from
  the control samples (Phase 4) cannot remove; its background part is similar (15.2% against 17.8%).
- Proposed for Phases 3 and 4: the LLR-and-geometry particle classifier with the PID event features at score > 0.92, and a
  constraint of the background from the control samples.

### Comparison with the BNB note

The BNB CC2π±Np note (`report/internalDocs/Internal_Note_v2.pdf`, section 13, Table 16) quotes a detector term of 11.0%
(Recomb2 7.1%, WMX 5.5%, SCE 4.1%, LY attenuation 4.1%, LY down 2.0%, WMYZ 1.5%, LY Rayleigh 1.0%). It is the
fractional uncertainty of the predicted selected event count in the signal region, from seven Run-4d variations without
the wire-angle ones, in a preliminary study with the Run-4b open data. In the same convention (change of the predicted
selected events, beam-off unchanged, COMB) the two-pion selections here have:

| Selection | All eight variations | The note's variations |
|---|---|---|
| Deployed classifier, PID set, score > 0.92 | 17.0% | 14.0% |
| LLR-and-geometry classifier, PID set, score > 0.92 | 10.6% | 8.4% |
| `CC1mu2pi` | 8.8% | 6.6% |
| `CC1mu2piNewPID` | 7.2% | 6.8% |

The detector terms quoted above (21–34%) are those of the extracted total, in which the background subtraction and the
efficiency correction amplify the change of the predicted events by a factor of two to four (most for `CC1mu2pi`, B/S
2.8). Recomb2 leads in both analyses.
