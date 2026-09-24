# TODO — state at 2026-09-24 (session handover)

Last commit: e3e59fa. Working tree clean apart from `slurm/ensemble_jobs_2026-09-19.txt`
(job bookkeeping, commit with the next item) and the untracked NuWro generator dirs
(NuWro closure dropped; leave or delete).

## Immediate (mechanical, ~1 h)

1. **Final ensemble refresh.** All 600 members are in (the 5 transient failures were re-run:
   SLURM 3356048/3356049). Run, from `xsec_analyzer/`:
   `bash slurm/ens_statcov_all.sh` (incremental, ~40 min) then
   `python3 ../report/tools/ensemble_tables.py` (patches the coverage tables of the note and
   supplement, writes `data_release/ensemble_2026-09-20.tsv`). Then in `analysis_note.tex` remove
   the sentence "Two ensembles are still filling ($51$ and $46$ members); the others have
   $96$--$100$." and in the supplement/note prose replace "$46$--$100$ members" by "$100$ members"
   and "(−0.24 for the RHC total on 46 members)" by the refreshed value; check the change-log
   open item "Ensemble offset" numbers; recompile (note, supplement, change_log); commit.

2. **Verify the rebuilt binaries** (the reason for the restart): after `make`, check the mtimes of
   `bin/univmake` and `bin/UnfolderNuMI` (bare `make` was once a no-op), then re-run one release
   unfold (e.g. `slurm/unfold_incl_local.sh` FHC p_mu) and confirm `closure_summary.tsv` /
   `current_results.tsv` are unchanged (`python3 report/check_tables.py`, `check_staleness.py`).

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

4. **MCS momentum-scale term: propagate.** Released only as `cov_MCSscale.txt` for FHC p_mu (up to
   0.34σ per bin; integral −6.8 % / +2.1 %). Still to do: RHC and combined p_mu (same data-side
   recipe: `macros/mcs_scale_fakedata.C` on the RHC fake data), the other observables (enters via
   the p_mu acceptance), and inclusion inside the framework covariance rather than as a released
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
