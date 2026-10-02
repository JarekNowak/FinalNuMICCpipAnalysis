# Phase 1: particle classifier for the multi-pion selections

2026-10-01. Phase 1 of `report/MULTIPION_BNB_ADAPTATION_PLAN.md`. Full tables: `phase1_classifier.md`
(training on tracks with backtracked purity > 0.5) and `phase1_classifier_unmatched.md` (tracks with
purity <= 0.5 kept as "other"); all numbers in the matching `.json` files.

## What was done

- Training dumper `xsec_analyzer/macros/mp_pid/dump_pid_tracks.C`. It reads the raw overlays directly,
  because the processed files do not carry the chi2, PIDA, truncated dE/dx, calorimetric energy or
  end-spacepoint vectors. It applies the event preselection of the selection (one slice, SCE-corrected
  vertex in the reco fiducial volume outside the dead z region) and the pion-candidate pool of
  `CC1mu1piXp` (generation 2, track score >= 0.3, LLR and Bragg values in range). Per track it writes
  63 PID features, the backtracked truth and the scores of the six BDTs the selections use today.
- Check of the dumper on the 30k-event Run-1 FHC slice of Phase 0: all 2,101 pion candidates have the
  LLR, vertex distance, length and BDT score of the selection exactly, and all 3,747 muon candidates
  are the same tracks.
- The 14 release overlays dumped (SLURM 3450911, `slurm/slurm_mp_pid_dump.sbatch`): 3.13 million
  tracks in `/data/uboone/processed/mp_pid`.
- `xsec_analyzer/scripts/mp_pid_train.py`: XGBoost, four classes (muon, pion, proton, other), balanced
  by weight. Trained on run periods 1, 2 and 4 of both horn modes, tested on the held-out periods 3
  (RHC) and 5 (FHC). Missing per-plane values are the sentinel -9999, not NaN.
- `xsec_analyzer/scripts/mp_pid_export_rbdt.py`: the model in the format of ROOT's
  `TMVA::Experimental::RBDT`. RBDT reproduces XGBoost exactly (largest difference 0 on 20,000 Run-5
  tracks), so the selection can evaluate the model in C++; the TMVA retraining of the plan is not needed.

## Results on the held-out runs

Population: the selection's pion-candidate pool, as CC1mu2pi counts pions (not the muon candidate,
vertex distance <= 9.5 cm), with track score >= 0.5 and length >= 5 cm.

ROC area, tracks with purity > 0.5:

| Pion identification | π vs rest | π vs p | π vs μ |
|---|---|---|---|
| New classifier (purity > 0.5 training) | 0.954 | 0.985 | 0.895 |
| `mp_pion_bdt` split, as CC1mu2pi applies it | 0.888 | 0.945 | 0.804 |
| `mp_pion_bdt` single | 0.886 | 0.942 | 0.804 |
| Inclusive pion BDT | 0.818 | 0.961 | 0.558 |
| MIP BDT | 0.718 | 0.957 | 0.271 |
| LLR score | 0.670 | 0.893 | 0.251 |

The gain is the same in FHC and RHC (0.953 and 0.954) and for π+ and π− (0.953 and 0.955); for
uncontained tracks it is 0.880 against 0.817.

At the pion efficiency of the CC1mu2pi cuts (91.8%), all tracks of the pool counted:

| | Pion purity, matched tracks | Pion purity, all tracks | Protons passing | Muons passing | Unmatched passing |
|---|---|---|---|---|---|
| `mp_pion_bdt` split (current) | 37.3% | 27.6% | 13.4% | 52.9% | 44.9% |
| New, purity > 0.5 training | 53.4% | 40.1% | 4.4% | 34.9% | 29.4% |
| New, unmatched as other | 54.8% | 43.3% | 4.1% | 31.8% | 22.8% |

Pion efficiency by true momentum at that working point (unmatched-as-other model, current BDT in
parentheses): 91.8% (88.8%) at 0.10–0.175 GeV/c, 95.6% (96.6%), 93.9% (96.0%), 91.8% (93.5%),
89.8% (88.8%) and 86.0% (84.1%) in the bins 0.175, 0.25, 0.35, 0.5, 0.75 and above. The new classifier
does not sculpt the momentum distribution more than the current one and is better at the threshold.

## Findings

1. Tracks without a truth match (purity <= 0.5, mostly cosmic rays of the overlay) are 23.5% of the
   candidate pool. Training on them as "other" lowers their pass rate from 29.4% to 22.8% and is the
   candidate model for Phase 2 (`/data/uboone/processed/mp_pid/mp_pid_rbdt_unmatched.root`). It stopped
   at the 2,000-tree limit with the validation loss still falling, so a longer training has some room.
   Phase 2 replaced it by the same training on the 43 inputs every sample carries
   (`mp_pid_rbdt_unmatched_common.root`; `PHASE2_SUMMARY.md`).
2. Muon against pion is the hard separation: 0.895, and 0.73 below 20 cm. About a third of the true
   muons in the pool (the muon candidate itself is excluded) still pass at 91.8% pion efficiency.
3. The largest inputs by gain are the proton chi2, the LLR, the track score, the MIP Bragg likelihood,
   containment, the proton-hypothesis energy, the number of descendants and the length. Length, range
   momentum and proton energy are kinematic quantities; their effect on the response is a Phase 3 item.

## Next (Phase 2)

1. Bind the 63 inputs in `AnalysisEvent` (guarded, as `swtrig_pre` is) and evaluate the RBDT model on
   each candidate when a multi-pion selection asks for it; store the four class scores with the
   candidate bookkeeping.
2. Compare the leading inputs (proton chi2, LLR, track score, MIP Bragg) between beam-on data and the
   prediction in the opened control regions before relying on them.
3. Assignment of tracks to μ, N × π and p, the event classifier, and the working point chosen on the
   expected uncertainty of the one-bin total.
