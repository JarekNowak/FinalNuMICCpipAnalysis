# Plan: a Sigma0 selection in the xsec_analyzer framework

Ported from `/home/t2k/nowak/MicroBooNE/Sigma/pelee` (Niam Patel's PeLEE-based
Sigma0 selection). That code already reads the same PeLEE ntuple structure as the
framework, so this is a port of physics into the framework's class interface, not a
rewrite of the reconstruction.

## What the existing selection does

Signal: `nu_pdg == -14 && ccnc == 0` with a Sigma0 (PDG 3212) in `mc_pdg`, true
vertex in the fiducial volume. The decay chain is Sigma0 -> Lambda gamma,
Lambda -> p pi-, so the reconstructed topology is muon + proton + pion + photon.

Cut chain, nine cumulative stages:

    AllEvents -> FiducialV -> CosmicVeto -> NShowers -> MuonCandidate
              -> ProtonCandidate -> PionCandidate -> NTracksCut -> BDTG

with `topological_score > 0.20`, `0 < nshowers < 4`, `ntracks < 5`.

Candidate assignment:
  * muon   - longest track with trk_score > 0.8, dist-to-vertex < 20 cm,
             length > 10 cm, LLR PID > 0.7
  * p / pi - joint discriminant D = LLR/0.591 + (chi2_p - chi2_pi)/124.5,
             proton = argmin D, pion = argmax D. AUC 0.873 for true p vs pi,
             against 0.838 for chi2_p alone.
  * shower - the two PFPs closest to the vertex with trk_score < 0.5

Derived observables: Lambda mass from (p, pi), Sigma0 mass from (Lambda, gamma)
using the per-PFP shower energy `shr_energy_y_v`, Lambda and photon
vertex-pointing angles and impact parameters, and a pi0 di-photon mass for
background rejection.

Final discriminant: a BDTG over 49 inputs, evaluated with two-fold cross-validation
(folds trained on even/odd event number, each applied to the events it did not
train on). Per-fold ROC ~0.81.

Reported performance at BDTG > -0.90, scaled to RHC data POT 1.1526e21:
efficiency 12.1%, purity 0.108%, S/sqrt(S+B) 0.0537 - a ~20% gain over the
cut-based selection at equal signal.

## Framework gaps to close first

### 1. Nine unbound branches
All nine are present in the ntuples - verified directly on the signal file - but
none is bound by `Branches.hh`, so `AnalysisEvent` cannot see them:

    trk_pid_chipi_v     shr_energy_y_v    shr_px_v   shr_py_v   shr_pz_v
    pi0_energy1_Y       pi0_energy2_Y     pi0_gammadot          shr_moliere_avg_v

Add members to `AnalysisEvent.hh` and guarded bindings in `Branches.hh`, following
the pattern already used for `trk_avg_deflection_stdev_v` - `if (etree.GetBranch(...))`
so that samples lacking them still load. `trk_pid_chipr_v` is already bound.

### 2. A different fiducial volume
Sigma0 uses x [3.00, 253.35], y [-112.53, 114.47], z [3.1, 1033.0] with the same
675.1-775.1 dead region, which is much looser than the CC1mu1piXp box
(10-246, +/-101, 10-986). `SelectionBase` already supports this per selection via
`define_reco_FV()` / `define_true_FV()`, so no framework change is needed - just do
not inherit the CCpi numbers by accident.

### 3. Registration
Add the header include and an `else if (selection_name == "Sigma0")` branch to
`SelectionFactory.cxx`.

## Scope question that shapes everything else

This is not a differential cross-section measurement. At 0.1% purity and 12%
efficiency the existing analysis produces a **cross-section upper limit**
(`CrossSection.C`, `Sigma0_CrossSection_UpperLimit.pdf`), whereas the framework is
built end-to-end for differential extraction: `univmake` builds systematic universes
per analysis bin, `UnfolderNuMI` regularises and unfolds them.

Nothing in the framework computes a limit - no Feldman-Cousins, no CLs. So the
honest framing is:

  * **In scope, and worth doing.** The selection itself, the cut flow, the truth
    categorisation, the reco/true observables, and above all the systematic
    universes. The framework's flux/xsec/reinteraction/detector machinery is exactly
    what a limit needs for its background prediction, and it is a lot of
    infrastructure to get for free.
  * **Out of scope.** The limit-setting step. That stays in `CrossSection.C`,
    reading the framework's covariance rather than recomputing it.

An unfolded Sigma0 spectrum would be meaningless at this purity, so I would not
build bin configs aiming at one. A single-bin scheme is the right target: it makes
`univmake` produce a total-rate prediction with a full covariance, which is precisely
the input to a limit.

## Proposed work

**Step 1 - branches.** Add the nine members and guarded bindings. Verify against
both the hyperon signal file and an RHC overlay, since the overlay is what the
background prediction is built from.

**Step 2 - `Sigma0.hh` / `Sigma0.cxx`.** Implement the eight pure virtuals:

    define_constants()        FV boxes, PID thresholds, the D-discriminant scales
    define_signal()           nu_pdg==-14, ccnc==0, 3212 in mc_pdg, true vtx in FV
    selection()               the nine-stage chain, storing the candidate indices
    compute_reco_observables()  Lambda mass, Sigma0 mass, pointing angles, impacts,
                                pi0 mass, inter-candidate distances and angles
    compute_true_observables()  true Sigma0 / Lambda kinematics where available
    categorize_event()        signal / other hyperon / RES pi0 / DIS / other / cosmic
    define_output_branches()  everything above, plus the BDT inputs
    define_category_map()     labels and colours for the stacked plots

Structure it on `CC1mu1piXp`, which already loads TMVA readers, so the BDTG fits the
existing pattern rather than needing new machinery.

**Step 3 - BDT.** Copy the trained folds from `pelee/SigmaBDT_A` and `SigmaBDT_B`
into `booster_decision_tree/`, and evaluate by event-number parity exactly as the
original does. The cross-validation must be preserved: applying a fold to the events
it trained on would inflate the ROC.

**Step 4 - validate the port.** Run `ProcessNTuples` on the signal file and compare
the cut flow stage by stage against the original's table (650 raw signal at
NTracksCut, 502 after BDTG). Any disagreement is a porting bug, and this is the step
that makes the port trustworthy - it should be done before any config work.

**Step 5 - configs.** `file_properties` for the hyperon signal, the five RHC overlay
runs, dirt and fake data; a single-bin scheme; an xsec config. Then `univmake` for
the covariance.

## Things to check with Niam before starting

1. **The fake-data slot is the signal file itself.** `SigmaSelection.C` uses the same
   NuWro hyperon sample for both sample 0 and sample 4, with a comment that a real
   closure test needs a different-generator hyperon ntuple that does not exist
   locally. So the current "data" is not an independent test. Worth confirming that
   is still the intent.
2. **Dirt is Run1-RHC only**, scaled up to full RHC POT, as the config comments note.
3. ~~The signal file POT.~~ Checked: the SubRun tree has 10352 entries summing to
   2.17640e23, matching the config exactly. No action needed.
4. **Whether a framework port is wanted at all**, given the limit-setting stays
   outside it. The gain is the systematics machinery; the cost is maintaining the
   selection in two places unless the original is retired.
