# TODO — state at 2026-09-24 (updated end of second session)


Last commit: see git log. Untracked NuWro generator dirs remain (NuWro closure dropped; leave or delete).

## Done 2026-09-24 (second session)

1. **Ensemble refresh: DONE** (1e64f8f). All six ensembles are complete at 100 members. The first
   refresh exposed a race: concurrent ensembles re-threw the same `ens/fakedata_*_t<T>.root`, and
   10 members read truncated trees. They were quarantined in `ens/corrupt_race_2026-09-20/` and re-run.
   The throw macros now write to a temporary file and rename it. Final: widths 0.92–1.07, offsets
   0.4–1.1% at ≤2.0σ (the RHC total is −0.81% at 2.0σ).
2. **Rebuilt binaries verified**: make is a true no-op (sources older than binaries); an FHC p_mu
   re-unfold reproduces all 38 released covariances and the closure sidecar bit-for-bit;
   check_tables and check_staleness are clean.

## Done 2026-09-24/25 (third session)

- Inclusive analysis note rewritten in analysis order, without em dashes, self-description or release
  history (3aeaf8a); stale pre-Run-1-fix values corrected (see change log, "Values corrected while
  rewriting the analysis note"). Combined generator weights fixed (8.857e20 -> 7.766e20) and the result
  set re-unfolded. MCS term shown to be exactly zero for the non-p_mu inclusive observables.

## Follow-ups from that session

- Proton-tagged theta_p / theta_pip COMB study extractions: their combined generator files
  (gen2d/*_1p_ext_comb_fte.root) were regenerated with the corrected weights (-0.8%); re-unfold them
  (UnfolderNuMI on the existing univmakes) and refresh the proton-tagged note tables.
- Proton-tagged ensembles (SLURM 3427889-91, 3x100: 1p_fhc5 Whad, 1p_fhc5 dpt2bin, 1p_comb pn2bin):
  when done, `bash slurm/ens_statcov.sh 1p_fhc5 Whad` etc., then ensemble_stat_pulls.evaluate('1p_fhc5_Whad').
- The technical supplement and proton-tagged note have not had the same style pass as the note.

## Open analysis items (need a decision or work)

3. **FHC Run-1 residual (gating for unblinding).** On the corrected exposure the frozen global test
   passes (8.1/8) but period by period FHC Run 1 sits at data/pred ≈ 1.0 while every other period
   sits 12–25 % high: a Run-1 factor 0.869 ± 0.022 is preferred (Δχ² = 37; a free factor on any
   other period gains ≤ 15). Cause not identified. Checks: (a) the 2.192e20 POT against the run
   list actually in `beam_on/neutrinoselection_filt_run1_beamon_beamgood.root` (963 runs, 82 200
   subruns; the ntuple carries no per-subrun POT — needs the beam database / run list);
   (b) completeness of that file against the good-run list; (c) the Run-1 overlay normalisation
   (MC POT 2.3282e21) and whether the Run-1 MC models the Run-1 detector conditions; (d) EXT is
   not implicated (cosmic region 0.95 ± 0.17). Documented in cr_data_note (tab:perperiod_corrected),
   analysis_note §cr_outcome + status box, change_log open items.

4. **MCS momentum-scale term: propagate.** Now released for p_mu in all three configurations
   (`cov_MCSscale.txt`; `report/tools/mcs_eval.py CFG`). FHC and combined are ≤0.35σ. RHC is ≤0.32σ
   except the open top bin, which is +21 %/−44 % (0.85σ): a candidate for the p_mu binning or the
   prior discussion. Other inclusive observables: the term is exactly ZERO (2026-09-24): the selection has no
   reco-p_mu cut and no other inclusive observable reads the muon momentum, so the scaled fake data give
   bit-identical results (FHC5 cosθμ/cosθπ/p_π/total checked; SLURM 3427873 cancelled after 9/36).
   Still to do: the proton-tagged TKI observables (δp_T, δα_T, p_n use the muon momentum vector; needs
   mcs_scale_fakedata on the w/ fake data), and inclusion inside the framework covariance rather than as a released
   add-on. Decide whether ±5 % is the right prior (the note argues from the −10 % estimator bias).

5. **Proton-tagged model sensitivity.** Low-imbalance δp_T and p_n bins are −1.1/−1.2σ under the
   coherent GENIE variation (truth moved 1 %; the shift enters through response/background). Before
   those bins are quoted: carry it in the covariance or exclude it with a validated independent
   generator. Remaining proton-tagged prerequisites: inverted-proton-PID control region; W_had edge
   scan; proton-tagged ensembles (statistical coverage); candidate-mistag bound for θ_p; 2D
   effective-rank reassessment.

6. **Independent-generator closure** still pending: the NuWro overlay is unvalidated (rate 1.76×
   GENIE vs ~1.6 implied) and was dropped; needs a validated NuWro or GiBUU overlay.

7. **Absolute flux normalisation** to be confirmed with the flux contact (unblinding prerequisite).

8. **Beamline-geometry flux term** (1.6–3.3 %) bounded but not in the covariance; beam-off
   gate-ratio and Run-2 stand-in terms (≲0.4 %) likewise.

## Reviewer "strongly recommended" items not done

9. W_had → W_cal rename (your call); larger figure legibility pass (the 2D maps and montages are
   dense at page width); reviewer's request for a run-by-run beam-on/beam-off table in the CR note.

## Housekeeping

10. Study trees on the old exposure that were *not* rebuilt: the 0.50-criterion candidates other
    than θ_p (rebuild_c50; conclusions are conditioning-based and do not depend on the throw).
11. Old-exposure artefacts kept for the record: `xsec-ana-fakedata_alt*_fhc_run*_pot3283.root`,
    `univmake_backup_job3351878`, `_job3352400`, `univmake_backup_study_pre_run1fix`, and the
    pre-fix background-fit templates `bkgfit/templates_run1old.root` — delete once the release is
    frozen.
12. Memory files to trust: `release-2026-09-19-state`, `hardcoded-run-scales`, `far-sidebands`
    (corrected), `altmodel-closure`, `nuwro-closure-fails` (MCS section superseded by the data-side
    result). After any exposure/flux change grep macros for 0.14101 / 3.283e20 / 9846635 / 6.60865.

## Later (before unblinding; beam-on files not used yet)

13. **Beam-on run lists vs exposure (Run 4/5).** Our beam-on files for Run 4b, 4d and 5 are event-level
    SUPERSETS of the good-run files in the FNAL copy `/data/uboone/temp/custom_pelee_ntuples`
    (+164 runs / +5.5 %, +52 runs / +1.9 %, +166 runs / +5.5 %; Run 4d equals the folder's `_all` file;
    good-run list in the folder's `selectGoodRuns.cc`). Check whether their POT/trigger values
    (Run 4 2.075e20 / 4,131,149; Run 5 2.231e20 / 5,154,196) were counted on the good-run lists; if so,
    apply the good-run filter to the beam-on (and beam-off) files. The per-directory POT/trigger text
    files did not arrive with the rsync -- re-sync them first. Run 1 is settled (user, 2026-09-24):
    keep our 829-run beam-good file, ignore the folder's 923-run file.
