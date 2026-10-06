# Phase 0: baseline of the NuMI multi-pion selections

2026-10-01. Phase 0 of `report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md`. The full tables are in
`phase0_tables.md` (made by `xsec_analyzer/scripts/mp_phase0_tables.py` from `phase0_baseline.json`);
the cut-flow alone is in `phase0_cutflow.tsv`. These numbers replace the raw Run-1 counts of
`report/other_notes/multipion_note.tex`.

## What was done

- Selection code (`CC1mu1piXp.hh/.cxx`). Selections with N > 1 now write:
  - every counted pion candidate, longest first: index, containment, length, range momentum (muon and pion hypothesis), MCS momentum, cos θ about the reco ν direction, LLR, BDT score, vertex distance, and the backtracked PDG, momentum and cos θ;
  - the true primary charged pions, hardest first;
  - true pion counts above 0.10 and 0.175 GeV/c, protons above 0.3 GeV/c, and kaons;
  - the ≥ 1p signal flag `MC_Signal_1p` (decision D1);
  - a truth topology class `mc_topology`;
  - per-stage pass bits `cutflow_bits`.

  Nothing is written for N = 1. On a 30k-event Run-1 FHC slice, all 262 existing branches are identical event by event before and after the change. The 52 new branches appear only under `CC1mu2pi_` and `CC1mu3pi_`.
- Reco opening-angle cut (decision of 2026-10-01). CC1mu2pi and CC1mu3pi now require θ(μ, longest pion candidate) < 2.6 rad, the reco counterpart of the truth cut in the signal definition (`apply_opening_angle_cut()` in `CC1mu2pi.hh` and `CC1mu3pi.hh`; the angle is stored as `mu_leadpi_opening_angle`). The 1π selection is unchanged.
- Reprocessing. 28 inputs (14 overlays, 7 beam-off, 7 dirt) were processed in one pass with `CC1mu1piXp,CC1mu2pi,CC1mu3pi`, into `/data/uboone/processed/mp_phase0` (SLURM array 3439300, `slurm/slurm_mp_phase0.sbatch`). All 28 passed the branch and entry-count checks. No beam-on file was read. The first pass without the reco cut (array 3439270) is in `/data/uboone/processed/mp_phase0_v1`; its tables are kept as `phase0_*_v1.*`.
- Normalisation. `scripts/mp_phase0_baseline.py` applies the framework rules: per-run data POT over the distinct `summed_pot` of the overlays and dirt, and 0.98 × beam-on triggers over beam-off gates. Exposures, triggers and gates come from `configs/file_properties_numi_comb_w.txt`. The `summed_pot` of every reprocessed file equals its release copy; the Run-1 dirt copies differ by 2 ppm from float accumulation.
- Check. The inclusive selection in the same pass reproduces the released final cut-flow row:

  | Config | Signal | ν background | Beam-off |
  |---|---|---|---|
  | FHC | 832.9 (released 833) | 560.3 (560) | 23.2 (23) |
  | RHC | 1015.4 (1015) | 719.8 (720) | 23.2 (23) |

## Baseline at the full exposure (selections of August plus the reco opening-angle cut)

| | Generated signal | Selected signal | Selected total | Efficiency | Purity |
|---|---|---|---|---|---|
| 2π FHC | 1210 | 221 | 836 | 18.2% | 26.4% |
| 2π RHC | 1672 | 303 | 1168 | 18.1% | 26.0% |
| 2π combined | 2882 | 524 | 2004 | 18.2% | 26.1% |
| 3π FHC | 204 | 21.7 | 203 | 10.6% | 10.7% |
| 3π RHC | 347 | 32.8 | 314 | 9.5% | 10.4% |
| 3π combined | 551 | 54.5 | 517 | 9.9% | 10.5% |

Purities include beam-off and dirt. The reco cut removes 1.2% of the selected 2π signal, 40% of its beam-off and 48% of its dirt; without it the combined 2π sample was 530 signal events in 2132 (18.4%, 24.9%) and the 3π sample 54.8 in 543 (10.0%, 10.1%).

## Findings

1. 2π background. In the combined configuration the largest backgrounds are:
   - CC1π± with a second counted pion, 22% of the sample;
   - events with a π0, 21% (no shower veto is applied);
   - CC≥3π±, 7%;
   - NC, 6%;
   - beam-off, 2.5%.

   Adding the 1π shower veto would raise the purity to 32% and lower the efficiency to 12.6%.
2. 3π. The 3π sample is 10.5% pure: 54.5 signal events among 517. The largest background is a 2π event with a third counted candidate, 30%. The current selection cannot support a 3π measurement.
3. Pion assignment in selected 2π signal events:
   - Only 60% have both counted candidates backtracked to true charged pions.
   - The longest candidate is the true leading pion in 52%.
   - The stored highest-LLR candidate is the true leading pion in 51%.
   - For 3π: 45% have all three candidates true pions, and the longest candidate is the true leading pion in 45%.

   Observables ordered by pion cannot be used with the present candidate choice. This is the case for the Phase 1 and Phase 2 classifiers.
4. Uncontained pions. 23% of selected 2π signal events count an uncontained pion, which has no range momentum.
5. Signal-definition details (combined configuration):
   - The truth cut θ(μ, leading π) < 2.6 rad removes 2.2% of the generated 2π signal (63.5 events). Without a reco counterpart 9.1 of them were selected as background with a signal topology; the reco cut on the longest candidate leaves 4.7.
   - Pions all above 0.175 GeV/c: 65% of the generated 2π signal and 53% of the 3π.
   - ≥ 1p (> 0.3 GeV/c): 81% of the 2π signal and 75% of the 3π. The current efficiency for the 2π ≥ 1p subsample is 17.7%.
6. Exclusivity:
   - The true 1π, 2π and 3π signals never overlap.
   - The reco 2π and 3π selections overlap: 121 events in the combined configuration, 23% of the 3π sample. This happens because 2π uses the split pion BDT and 3π the single one.
   - The 2π selection overlaps the blind 1π signal region: 311 events (15.5% of the 2π sample), containing 19.7% of the selected 2π signal. Decision of 2026-10-01: the overlap is kept, and the overlap events are opened only after the 1π signal region.
7. Opened control region (D3). 79% of the predicted 1π multi-π control-region events (224 of 283, combined) fall in the 2π selection.
   - They are 11% of the 2π sample and contain 15% of its selected signal.
   - The beam-on data of that region have been compared with the prediction (`cr_data_note.tex`).
   - The rest of the 2π region remains blind.
8. Side finding for the 1π documents, now corrected on branch `fix/dirt-scale`. The count-level and control-region macros scaled the dirt overlay by the data POT over the overlay POT (0.081020 FHC, 0.071666 RHC) instead of over the dirt sample's own POT, mostly with the 0.65 weight applied twice.
   - The framework uses the dirt POT, so the extraction is unaffected.
   - Signal region (commit 3177067): dirt 6.6/9.5/16.0 instead of 1.1/1.0/2.1; FHC purity 58.5% rather than 58.8%.
   - Control regions (2026-10-01, `logs/crdirtfix/`): dirt seven to twenty times larger; the 1π multi-π region total is 282.7 events (combined).
