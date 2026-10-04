# Plan: can the CC0pi flipped-track method remove backward-going muons here?

Written 2026-08-18. Sources: `report/internalDocs/BNB_CC0pi_InternalNote_v6.pdf` (Sec. 3,
Fig. 4) and `rutgersnu/xsec_analyzer@cczeropi_sel_mods`, `cc0pi/cc0pi_analyzer.C:1741`.

## What the CC0pi analysis does

Pandora sometimes reconstructs a track with the wrong direction relative to the vertex. The
CC0pi note is explicit that the usual culprit is **vertex** mis-reconstruction rather than the
track itself, and that it is overwhelmingly a **single-track** problem: a lone muon with no
second prong has no vertex "vee", so its direction carries a 180 degree ambiguity. Their
Fig. 4a shows the flipped fraction reaching ~14-16% for one-track events and falling steeply
with track multiplicity.

Their fix uses calorimetry: `trk_bragg_mu_fwd_preferred_v` flags whether the Bragg-peak
likelihood prefers the reconstructed direction or the reverse. Cutting on that alone removes
too much signal, so they require it *in coincidence* with a poor forward-muon fit:

```cpp
// cc0pi_analyzer.C:1741
if ( track_length_->size() == 1
     && trk_bragg_mu_fwd_preferred_v_->at(muon_candidate_idx_) == 0
     && 6. < track_chi2_muon_->at(muon_candidate_idx_) ) {
  muon_mom = LOW_FLOAT;          // reject the event
}
```

Separately, for **uncontained** muons they reject `cos(theta)_reco < -0.9` outright,
reporting negligible signal loss.

## Why it barely applies to us, measured rather than assumed

Our signal definition requires exactly one muon **and** at least one charged pion, so a
selected event always has >=2 tracks. **The single-track condition never fires.** The method
as written is a no-op on this analysis.

The residual question is whether the underlying problem survives at >=2 tracks. Measured on
CC1mu1piXp selected signal, FHC Runs 1/2/4 (8712 events), via
`macros/flipped_muon_check.C`:

| | rate |
|---|---|
| hard flip (reco cos < -0.5, true > +0.5) | **0.22%** |
| sign flip (reco < 0, true > 0) | 0.83% |
| reco cos < -0.9 | 0.24% |
| **true** cos < -0.5 (genuine backward muons) | 1.81% |

So the flip rate is ~0.2%, roughly two orders of magnitude below the CC0pi single-track case,
exactly as their multiplicity trend predicts. Note also that genuine backward-going muons
(1.81%) outnumber flipped ones by ~8:1, so a blanket `cos_reco < -0.9` rejection would remove
more real signal than mis-reconstruction, the opposite of the CC0pi situation.

## Does it help with EXT rejection?

Part of the surviving beam-off sample is thought to be broken cosmic muon tracks, where the
break is reconstructed as a vertex. That leaves a directional signature, and it is present:

| reco cos(theta_mu) | EXT | MC signal | enrichment |
|---|---|---|---|
| < 0 | 31.8% | 4.1% | **7.8x** |
| < -0.5 | 9.1% | 1.4% | 6.4x |

The enrichment is real, but the arithmetic does not support a cut. At the final cut the FHC
sample is 974 signal, 128 EXT, 1746 predicted. A `cos < 0` requirement removes ~41 EXT and
~40 signal, close to one for one; and since EXT is subtracted data-driven and unbiased,
removing part of it buys only a smaller *statistical* uncertainty on the subtraction,
sqrt(128)~11 -> sqrt(87)~9 events. Spending 40 signal events to gain ~2 is a net loss.

**The ceiling on any cosmic cut**: EXT is 7.3% of the selected sample and already subtracted
without bias, so removing ALL of it for ZERO signal loss would gain ~11 events of statistical
uncertainty on 1746. Any cut costing more than ~1% of signal is a net loss however
efficiently it tags cosmics.

The `fwd_preferred && chi2>6` tagger is more surgical than an angular cut, since it fires on
tracks whose calorimetry disagrees with the assigned direction rather than on backward tracks
as such, and so would spare the genuine backward muons that outnumber flipped ones ~8:1. It
is worth measuring under the same 11-event ceiling, but it cannot become a significant
improvement to the measurement.

## What is still worth doing

Not the cut. Two smaller things:

1. **Bound the effect and record it.** 0.22% is small enough to state as a negligible-bias
   argument in the note, which is worth having when a reviewer asks why we do not apply the
   CC0pi correction. Cost: the measurement above, already done.

2. **Check the theta_mu tail specifically.** The theta_mu resolution panel (Sec. 8.7, panel g)
   shows the mean residual drifting to -0.15 rad above ~2.5 rad. That is the backward region,
   and it is the same observable whose FHC A_C is rank-1. If flipped tracks populate that tail
   disproportionately, cleaning them could sharpen the top theta_mu bin. Worth an hour, not more.

## Steps, if we do step 2

1. **Pass the branches through.** `trk_bragg_mu_fwd_preferred_v` and `trk_pid_chimu_v` are
   present in our PeLEE ntuples (verified) but not written by `ProcessNTuples`. Add them as
   optional binds, exactly as `trk_avg_deflection_stdev_v` was added in `a82fe82`. The
   selection must stay bit-identical: check the cut-flow counts before and after.
2. **Measure the tagger on our sample.** For selected signal, cross-tabulate
   (`fwd_preferred == 0` && `chi2_mu > 6`) against the truth-level flip. Report the efficiency
   for tagging true flips and the signal loss. Against 0.22% contamination the tagger has to
   be very pure to be worth applying.
3. **Decide on the evidence.** Apply only if it removes a clear majority of the 19 hard-flip
   events while costing well under ~0.5% of signal.
4. **If applied, re-run the response** and compare theta_mu and cos(theta_mu) migration
   diagonals and the A_C singular-value ratio before and after. That is the only outcome that
   would change a physics conclusion.
5. **Systematics.** A new reco-level cut needs the detVar samples re-run through it. That is
   the real cost, and it is why steps 1-3 should gate step 4 rather than the reverse.

## Recommendation

**Do step 1 (record the bound), and step 2 only if the theta_mu tail matters for the A_C work.**
Do not port the CC0pi cut: its single-track condition cannot fire here, and its
`cos_reco < -0.9` companion would cut ~8 genuine backward muons for every flipped one. The
CC0pi analysis needed this because a 0p final state gives a lone track; requiring a pion
removes the ambiguity for free.
