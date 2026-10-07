# Plan: adapting the BNB CC2π±Np note to the NuMI 2π and 3π measurements

Written 2026-10-01. Source: `report/internalDocs/Internal_Note_v2.pdf`, N. Majeed (KSU), "Study of production
mechanisms of CC events with two charged pions and at least one proton in the final state",
May 2026, 118 pp. (78 pp. analysis, the rest event displays). Every number below is quoted from
that note, from our own documents, or was checked in the code and ntuples today.

## 1. What the BNB note does

| Element | BNB note |
|---|---|
| Samples | MCC9.10 "surprise" BNB Run 4a–5 (Pandora + NuGraph); dedicated `cc1mu2piNp` signal-enhanced overlays (Runs 1–5) for training; open data Run 4b ($4.03\times10^{19}$ POT); centrally produced NuWro overlays Run 4a/4c/5; detVar Run 4d |
| Signal | $\nu_\mu$ CC only; 1 μ (KE > 10 MeV), exactly 2 π± (KE > 10 MeV), ≥ 1 p (KE > 30 MeV); no π0/kaon condition stated |
| Pre-selection | nslice = 1, CosmicIP > 10, topological score > 0.1, vertex in FV, generation 2, vertex distance < 4 cm, track score > 0.5 |
| Particle ID | three XGBoost (DART) classifiers, μ / π / p each against all other tracks, trained on truth-matched tracks of the enhanced samples, 70/30 random split |
| Event ID | per-track thresholds (μ > 0.8, π > 0.6, p > 0.8); the track assignment with the largest product of scores is the event hypothesis; a second "meta" BDT on the assigned scores, their product, per-type sums and counts above threshold; AUC 0.88 |
| Final cuts | event BDT ≥ 0.6, opening angle < 160°, μ start contained, all π and p fully contained |
| Performance | event BDT takes purity 1.1–1.2% → 41–46% and efficiency 86–90% → 9–11%; containment takes efficiency 6.8–8.1% → 3.9–5.6%; final purity 47–58% |
| Observables | invariant masses $M_{p\pi_1}$, $M_{p\pi_2}$, $M_{\pi\pi}$, $M_{p\pi\pi}$, $M_{\mu p\pi\pi}$ and the $p_T$ of each system; μ, p, π1, π2 momentum and angles |
| Sidebands | (1) event BDT < 0.4, per-run data/MC with χ² for single, two, three and four-particle quantities and the vertex; (2) CC1π±Np (N ≥ 2) selected with the particle BDTs |
| Systematics | one bin only: total 26.6%, of which open-data statistics 20.6%, detVar 11.0% (Recomb2 7.1%, WMX 5.5%), flux 8.1%, cross section 6.9%, MC statistics 5.9% |
| Fake data | GENIE closure; NuWro fake data with statistical and cross-section terms only: unfolded $0.45\pm0.12$ against NuWro truth 0.37 and tune 0.62 ($10^{-38}$ cm²/Ar) |
| Background constraint | conditional covariance (`ConstrainedCalculator`) with 1μ1π2p (0.25 < BDT < 0.6, purity 35%), 1μ3π (21 events, purity 62%) and 1μ3p (1043 events, purity 74.5%) sidebands; each alone lowers the uncertainty (27.3% → 21.7 / 25.9 / 25.2%); all combined give 21.5% but move the NuWro result to $0.51\pm0.11$, away from its truth; the note does not apply the constraint |

## 2. Where our multi-pion work stands

Checked in the code on 2026-10-01.

- `CC1mu2pi` and `CC1mu3pi` are thin subclasses of `CC1mu1piXp`, registered and compiled. Nothing
  downstream exists: no bin, univmake, xsec or systcalc configuration and no results. Every number
  predates the beam-frame fix (08-19), the CV-weight and EXT-gate fixes (09-03), the Run-1 exposure
  correction (09-19) and the 0.50 criterion (09-26).
- Truth signal: $|\nu_\mathrm{pdg}|=14$, $p_\mu>0.15$ GeV/c, exactly N charged pions counted at any
  momentum with the softest above 0.10 GeV/c, π0, kaon and other-meson veto, any nucleons,
  θ(μ, leading π) < 2.6 rad at truth only (reco opening-angle cut, shower veto and wire-gap cuts off).
- Reco: the inclusive chain with a looser pion ID (vertex distance ≤ 9.5 cm, `mp_pion_bdt`), N pion
  candidates with ≥ 1 contained and up to 2 uncontained.
- Only one reco pion is stored, the contained candidate with the highest LLR, while the truth
  candidate is the hardest pion. For N > 1 the reco and true $p_\pi$, $\cos\theta_\pi$ and
  $\theta_{\mu\pi}$ refer to different particles. No second or third pion kinematics, no masses, and
  every $\nu_\mu$CC non-signal event is in one category (`kUnknown`).
- `report/other_notes/multipion_note.tex` (Draft 0.1, 08-13) uses the Run-1 FHC overlay only, without EXT or
  dirt: 2π ε = 15.5%, P = 23.2% (579 signal); 3π ε = 6.1%, P = 11.6% (40 signal). The counts are raw
  overlay counts (weight spline × tune, MC exposure $2.33\times10^{21}$ POT), so the "670 (66) events
  at full exposure" of its abstract corresponds to about 190–225 (13–22) selected signal events at
  the FHC data exposure, before EXT. The threshold study keeps 90% (2π) and 82% (3π) of the signal at
  a softest-pion threshold of 0.10 GeV/c, and 70% and 49% at 0.175.
- `mp_pion_bdt` (TMVA BDTG, 8 inputs, Run-1 FHC, random split, evaluated on its own sample; the
  training dumper is not committed) has ROC 0.81 and already runs on every pion candidate.
- The June standalone `NPiSec` study (NuMI Run 1, exclusive 2π, 4-class XGBoost with π-vs-p AUC
  0.979) reached ε = 3.1%, P = 45.2%; the "exactly two pions" step alone took the efficiency from
  61.7% to 9.1%. Its POT, EXT scale and angle convention are stale.
- The multi-π control region of the 1π analysis (strict pion ID, effectively exactly 2π) has been
  compared with beam-on data (`cr_data_note.tex`): FHC 109 events against 109.6 predicted, RHC 136
  against 125.0. How many of these events a new 2π selection would keep is not known.
- Ntuple inputs: the raw NuMI PeLEE ntuples carry `trk_bragg_{p,mu,pion,mip}` per plane,
  `trk_pida_v`, `trk_pid_chipr_v`, `trk_pid_chipi_v`, `trk_trunk_dEdx_*_v`, `trk_calo_energy_*_v`,
  `trk_end_spacepoints_v`, `trk_llr_pid_*` and the deflection variables. None of the NuGraph fractions
  (`pfng2mipfrac`, `pfng2hipfrac`) is present. `AnalysisEvent` binds the deflection and Bragg-pion
  branches but not χ², PIDA, trunk dE/dx, calorimetry or end spacepoints.
- The framework already has `ConstrainedCalculator` and the `kSidebandRecoBin` reco-bin type, the
  same tool the BNB note used. Nothing calls it: `UnfolderNuMI` uses `MCC9SystematicsCalculator`.
- The Bragg-pion requirement (≥ 0.08), dropped from the released single-pion selections on
  `fix/bragg-pion` (2026-10-02), does not act in `CC1mu2pi` and `CC1mu3pi`: with `use_pion_bdt()` their
  pion identification is the BDT alone (`require_bragg_with_bdt()` false). On this branch
  `apply_bragg_pion_cut()` defaulted to true for the single-pion selection until `main` (version 1.7)
  was merged into this branch on 2026-10-06; the default is now false here as well.

## 3. Assessment

### Adopt

| Item | Why it helps us |
|---|---|
| Two-stage identification: per-track μ/π/p classifiers, best-assignment hypothesis, event classifier | The event classifier is where the BNB purity comes from (1% → 44%). Our selection counts pion candidates and keeps one; NPiSec lost 85% of its signal at "exactly two pions". For 3π, assigning 4–6 tracks is the whole problem |
| Assignment by the largest product of particle scores | Gives every track a role in one pass: the second and third pion and the proton candidate become available for observables, categories and the response |
| Per-track features: Bragg log-ratio p/MIP and MIP likelihood, LLR, trunk dE/dx, track score, length, deflection spread, end spacepoints, reco containment | All but NuGraph exist in our raw ntuples. Their pion classifier ranks the Bragg p/MIP ratio first by a factor of five |
| Background categories by topology (CC1π with extra protons, CC0π with ≥ 3 p, CC(N+1)π±, CCπ0, ...) | Today every νμCC background is `kUnknown`; no sideband can be designed or checked without these |
| Sidebands each aimed at one migration path, defined from the classifier outputs, with purity and efficiency of the targeted class quoted | Same structure as our transfer-factor tables; the 1μ3p sideband (74.5% purity, no signal-region overlap) is the cleanest example |
| Conditional-covariance background constraint, adopted only after alternative-model closure | Available in our framework. It is also the one lever on the (1 + B/S) amplification of the flux term identified in the 1π analysis |
| Centrally produced NuWro overlays as fake data | Our independent-generator closure is the open validation item of the 1π analysis. MCC9.10 NuWro overlays exist for BNB; whether a NuMI equivalent exists has to be checked |
| Run-by-run stability plots and an event-display scan of the EXT events that pass the selection | Cheap, and at four or more tracks beam-off events are the background the classifier has seen least |
| One-bin total first, differential second | The order the BNB analysis and our 1π analysis both followed |

### Adapt with changes

| Item | Change |
|---|---|
| Signal definition | Keep our conventions ($|\nu_\mathrm{pdg}|=14$, $p_\mu>0.15$ GeV/c, per-pion momentum threshold, π0/kaon veto, θ cut), not KE > 10 MeV and $\nu_\mu$ only; NuMI is 35–55% $\bar\nu_\mu$ |
| Training sample | No enhanced samples exist for NuMI. Use all 14 overlays, both horn modes, and hold out whole runs for testing. The event classifier must see the nominal background mixture, EXT and dirt included, with CV weights; the BNB meta-BDT was trained on enhanced samples only |
| Containment | Contained tracks are needed for any range momentum, not for angles. Decide per observable rather than requiring every π and p contained, which costs a third of the BNB efficiency |
| Pion ordering | The BNB note orders reco pions by length and truth pions by energy; we order truth by momentum and keep the highest-LLR reco pion. Order truth and reco by the same quantity, or use observables symmetric in the pions |
| Working point | Choose the event-classifier threshold on the expected total uncertainty of the one-bin total, not on efficiency × purity: the flux term scales as 17% × (1 + B/S) |
| Invariant masses and RES/DIS separation | Study products, not results, until a binning passes the 0.50 criterion and the trace rule. The BNB pion-energy response saturates (reco $E_\pi \lesssim 0.4$ GeV for true values up to 2 GeV; leading-pion residual peaked near −0.5), so the reco masses are squeezed towards threshold: $M_{p\pi\pi}$ is $1.34\pm0.07$ GeV at reco against $1.55\pm0.16$ (RES) and $1.65\pm0.27$ (DIS) at truth. The reco LLR test ($p\sim10^{-9}$) compares two MC samples through the same response; it does not show that the separation survives unfolding. This is the $W_{\pi p}$ limit of our supplement |
| Detector systematics | Check the knob-by-knob detVar of the multi-track selection (Recomb2 leads in the BNB note). Our Run-5 FHC Recomb2 and SCE samples are $\nu_e$ samples and stay excluded |

### Do not adopt

- KE > 10 MeV thresholds: most of that phase space is not reconstructed (BNB efficiency 4–6%), and
  the extrapolation is model dependent.
- The muon classifier as built: `is_contained` dominates it (importance 17015 against 1708 for the
  next feature), and the muon training tree has a different containment preselection from the
  others, so the classifier can learn the preselection. Any containment feature must be reco-level,
  with the same preselection for signal and background tracks.
- The combined-sideband constraint without closure (the note's own conclusion).
- The sequential-resonance table (1–5 events per chain).
- NuGraph features, unless NuGraph-processed NuMI samples turn up (it leads their proton classifier).

## 4. Plan

### Phase 0: definitions and an honest baseline

Status 2026-10-01: items 2–4 done, re-run with the reco opening-angle cut (D6); results and open choices in `report/planning/multipion/PHASE0_SUMMARY.md`.

1. Fix the signal family with D1 and D2 (section 5): inclusive Xp with a ≥ 1 p subsample, per-pion
   threshold 0.10 GeV/c; then the leading-pion convention and the exclusivity between 1π, 2π and 3π.
2. Store every pion candidate as vectors: range and MCS momentum, direction in the beam frame,
   length, containment, PID scores, backtracked PDG and true momentum. Store the true pions
   ordered consistently with the reco ordering.
3. Add topology categories to the multi-pion selections: CC0π, CC1π±, CC2π±, CC≥3π±, CCπ0, NC,
   $\nu_e$CC, out-of-FV, EXT, dirt.
4. Re-run `CC1mu2pi` and `CC1mu3pi` unchanged on all 14 overlays with EXT and dirt and the current
   normalisation. This gives the baseline cut-flow for FHC, RHC and combined that replaces the
   raw Run-1 counts of the multi-pion note.

### Phase 1: particle classifier

Status 2026-10-01: items 2–3 done with XGBoost, deployed through ROOT's RBDT instead of a TMVA retraining; item 1 moves to Phase 2, where the selection evaluates the model. Results in `report/planning/multipion/PHASE1_SUMMARY.md`.

1. Pass through `ProcessNTuples` and bind in `AnalysisEvent`: `trk_bragg_{p,mu,pion,mip}_v` per
   plane, `trk_pida_v`, `trk_pid_chipr_v`, `trk_pid_chipi_v`, `trk_trunk_dEdx_{u,v,y}_v`,
   `trk_end_spacepoints_v`, `trk_calo_energy_{u,v,y}_v`, guarded as the `swtrig_pre` binding is.
2. Commit a training dumper (the `mp_pion_bdt` one never was): truth-matched primary tracks
   (purity > 0.5, track score ≥ 0.5, length ≥ 5 cm), labels μ / π / p / other, with run, horn mode,
   charge and containment recorded.
3. Train a multiclass TMVA BDTG, which the selection can already evaluate in C++; keep XGBoost as a
   benchmark. Test on held-out runs. Compare with `mp_pion_bdt` and the inclusive MIP and pion BDTs
   on the same tracks, separately for FHC and RHC and for π+ and π−.

### Phase 2: assignment and event classifier

Status 2026-10-01: the classifier is evaluated in the selection (Phase 1 item 1) on the 43 inputs every ntuple production carries, and replaces the pion identification in study selections of the single-, two- and three-pion selections; items 1–4 done on 2026-10-06 (track assignment, event classifier, working point on the expected total uncertainty with all terms: score > 0.92 proposed, total 47% against 89% for `CC1mu2pi`). Detector robustness studied on 2026-10-07: the detector term (33%) is real (7.7% of it from the statistics of the variation samples); a particle classifier on the LLR score and geometry only halves its efficiency part (14% to 8%) at 13% efficiency. Results in `report/planning/multipion/PHASE2_SUMMARY.md`.

1. Enumerate the assignments of primary tracks to {μ, N × π, optional p, other}; keep the one with
   the largest summed log score and the runner-up.
2. Event features: the BNB set (assigned scores, product, per-type sums and counts above threshold)
   plus the shower count, topological and CosmicIP scores, and opening angles.
3. Train separate 2π and 3π event classifiers on the nominal mixture (overlay, dirt, EXT, CV
   weights), holding out whole runs.
4. Choose the working point on the expected total uncertainty of the one-bin total.

### Phase 3: selection performance and response

1. Cut-flow, categories and run-by-run stability for FHC, RHC and combined; scan the selected EXT
   events by eye.
2. Response of the candidate observables: one bin; $p_\mu$, $\cos\theta_\mu$; pion angles with the
   convention of Phase 0; the opening angle between pions; proton multiplicity. Apply the 0.50
   criterion and the trace rule with `report/tools/maximin_binning.py`.
3. Keep $M_{\pi\pi}$, $M_{p\pi}$ and $M_{p\pi\pi}$ as study products unless they pass.

### Phase 4: control regions and the constraint

1. Define the sidebands from the classifier outputs and quote the purity and efficiency of the
   targeted class and the overlaps:
   - proton taken as a pion: event score just below the working point;
   - pion candidates that are proton-like (the 1μ3p analogue);
   - one pion candidate too many (the ≥ 3π sample, D8);
   - inverted event score.
2. Draw one reco-level region map for the 1π, 2π and 3π analyses together:
   - A 1μ1π + protons sideband for the 2π analysis lies inside the still-blind 1π signal region.
   - Measure in simulation the overlap of the 2π signal region with the already-opened 1π multi-π
     region, and report it with the frozen procedure.
3. Add a calculator switch to `UnfolderNuMI` for `ConstrainedCalculator`, with the sideband reco
   bins after the ordinary ones. Check that it handles the per-run normalisation, the pooled detVar
   samples and the data-side MCS term.
4. Validate each sideband constraint, alone and combined, on CV pseudo-data, on the two reweighted
   models (GENIE u545, Δ→Nπ angular) and on NuWro if available. Adopt a constraint only if the
   alternative truth stays within 1σ and the result is not pulled towards the tune.
5. As a separate study on pseudo-data only, run the same constraint for the 1π analysis.

### Phase 5: extraction

1. Configurations for the 2π one-bin total in FHC, RHC and combined first; the 3π one-bin total in
   the combined configuration only as an exploratory study (D8). Generator predictions use the same signal definition.
2. Validation as for 1π: identity closure, CV closure, 100-member ensembles, reweighted models, and
   an independent generator if a sample exists.
3. Then the 2π differential observables that pass the binning rule.

Each univmake pass takes about 1.5–2 h per configuration. Ensembles multiply that by the number of
members.

## 5. Decisions

Taken 2026-10-01:

| # | Decision |
|---|---|
| D1 | Signal: inclusive CC2π±Xp (and CC3π±Xp), with a ≥ 1 p subsample for comparison with BNB |
| D2 | Per-pion threshold 0.10 GeV/c for the signal, totals and angles; pion-momentum bins only above 0.175 GeV/c, where the range estimator is unbiased |
| D3 | The 2π signal region is treated as blind. Its overlap with the opened 1π multi-π control region is not known and is measured in simulation (Phase 4, item 2) |
| D6 | The truth cut θ(μ, leading π) < 2.6 rad gets a reco counterpart: θ(μ, longest pion candidate) < 2.6 rad in CC1mu2pi and CC1mu3pi |
| D7 | The overlap with the blind 1π signal region is kept (no exclusivity cut); those events are opened only after the 1π signal region |
| D8 | (2026-10-02) The measurements stop at 2π; 3π stays exploratory. Events with three or more pion candidates are a control sample that constrains the 2π background, as the BNB note's 1μ3π sideband does. Four and five pions are not measurable at the full exposure (`report/planning/multipion/highn_feasibility.md`, `scripts/mp_highn_feasibility.py`): combined, 3.1 selected four-pion signal events at 4.8% purity and 0.1 five-pion events; P(reco bin \| true bin) is 0.84, 0.48 and 0.33 for 1, 2 and ≥ 3 pions |

Open:

| # | Item | Proposal |
|---|---|---|
| D4 | NuGraph inputs | proceed without them; ask whether NuGraph-processed NuMI samples exist |
| D5 | NuWro for NuMI | ask whether a centrally produced NuMI NuWro overlay exists (the BNB ones are under `/pnfs/uboone/persistent/users/uboonepro/surprise/`, not mounted here); it would also close the open validation item of the 1π analysis |
