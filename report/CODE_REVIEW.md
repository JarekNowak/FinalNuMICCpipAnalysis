# Code review — xsec_analyzer framework, macros and drivers

Review date 2026-08-16, branch `fix/systematics-warning`. Findings are ordered by how much
they could change a published number. Everything below was checked against the code and,
where a number is quoted, against a live run — not inferred from naming.

---

## 1. ~~The supporting tables use a different EXT model~~ — CORRECTED: a 2% EXT occupancy factor

**Status: FIXED (commit 9892a24). The original finding below was wrong; see the correction.**

> **Correction.** The eight per-run beam-off files are all *symlinks to*
> `beamoff_run1Andrun3.root` with identical trigger counts, so the framework's per-run sum
> reduces exactly to the macros' single scalar (`22667109/3821593 = 5.9313`, the hardcoded
> FHC value). The two are mathematically equivalent — I asserted a difference from the
> filenames without checking they resolve to one inode. The **real** discrepancy is that the
> framework applies the 2% NuMI beam-occupancy factor (`* 0.98`) and the macros did not, so
> EXT was 2% over-counted: 130.5 → 127.9 (FHC), 135.5 → 132.8 (RHC), 266.0 → 260.7 (comb).
> The cut-and-count σ does **not** move (under blinding it is background-independent);
> only the background breakdown shifts.

The framework reads **eight per-run EXT files** and scales each by that run's trigger ratio
(`SystematicsCalculator.cxx:613-627`, `bnb_trigs/ext_trigs * 0.98`):

```
run 1  xsec-ana-beamoff_fhc_run1.root      run 11  xsec-ana-beamoff_rhc_run1.root
run 2  xsec-ana-beamoff_fhc_run2.root      run 12  xsec-ana-beamoff_rhc_run2.root
run 4  xsec-ana-beamoff_fhc_run4.root      run 13  xsec-ana-beamoff_rhc_run3.root
run 5  xsec-ana-beamoff_fhc_run5.root      run 14  xsec-ana-beamoff_rhc_run4.root
```

Every supporting macro instead uses the **single retired combined file** with one scalar:

| macro | line | EXT source |
|---|---|---|
| `total_xsec_counting.C` | 33 | `beamoff_run1Andrun3.root` × `sc_ext` |
| `sideband_compare.C` | 84 | `beamoff_run1Andrun3.root` × `sc_ext` |
| `cutflow_yields.C` | 35 | same |
| `cutflow_yields_1p.C` | 33 | same |

So the cut-flow table, the sideband table and the cut-and-count cross section rest on a
different cosmic model from the differential result they are presented alongside. EXT is
not a rounding term: from a live run it is **130.5 / 759.3 = 17%** of the FHC background
and **266.0 / 1624.4 = 16%** of the combined background. A 10–20% error in EXT moves the
cut-and-count σ by roughly 1.5–3%.

The per-run files exist, so the fix is mechanical: give the macros the same per-run list
and per-run trigger scaling the framework uses.

---

## 2. Combined-mode dirt is scaled to FHC exposure only

**Status: FIXED (9892a24).** Combined dirt 1.3 → 2.2 events; σ unchanged.

```cpp
double sc_dirt = (std::string(mode)=="rhc") ? 0.071666*0.65 : 0.092402*0.65;
```
`cutflow_yields.C:39`, `cutflow_yields_1p.C:35`, `sideband_compare.C:47-49`,
`total_xsec_counting.C:86-88`.

A two-branch ternary over a **three**-valued mode: `comb` falls through to the FHC branch.
EXT in the same macros is handled correctly (`12.0898 = 5.9313 + 6.1584`), so within one
function EXT gets combined exposure while dirt gets FHC-only — dirt is ~44% low for `comb`.

Confirmed directly in the live output: **combined dirt = 1.3 events, identical to FHC dirt
= 1.3**, where it should be roughly the sum.

Numerically this is negligible today (1.3 of 1624 background events, 0.08%), so **no
published number needs revising**. It is worth fixing because the pattern is silent and
would bite hard if the dirt sample ever grew.

Two related points on the same lines: the `0.65` dirt normalisation factor is asserted with
no derivation anywhere in the repo, and a single **FHC run-1** dirt file
(`dirt_fhc_mcc9_run1_v28_all_snapshot`) is used for RHC and combined alike.

---

## 3. `θ` generator predictions are not reproducible from committed code

**Status: FIXED (9892a24)** — plus a further problem found while testing: the
`*_newg4_final.root` inputs are **stale** w.r.t. the coarsened binning (ppi 6 / cosθπ 5 /
θμπ 7 vs the analysis's 5 / 4 / 5), so `make_fte` on them would regress three observables.
Added a staleness guard and `add_theta_fte()`, which derives θ from the cosθ histograms
inside an existing FTE file and reproduces all twelve on-disk files bit-for-bit.

`generator_predictions/newg4/make_fte.C:7` builds exactly five observables:

```cpp
const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
```

but every FTE file on disk contains **seven** — `thetamu_fte` and `thetapi_fte` as well.
`run_theta_rollout.sh:9` says they "were produced by exact bin reversal (sigma conserved)",
and the reversal is verifiably correct (θ is the exact bin-reversal of cosθ), but **no
committed script performs it**. Regenerating the FTE files today silently drops both θ
observables, and the θ configs would then fail or fall back.

Fix: add the two reversed observables to `make_fte.C` so the chain reproduces what is on
disk.

---

## 4. `XSEC_FORCE_REBUILD` is opt-in and the cache has no invalidation key

**Status: FIXED.** Every cached `total_` subfolder is now stamped with a digest of the
inputs that produced it — normalisation scheme tag, active file list, data POT and trigger
counts, and the per-file simulated POT (which this analysis rewrites in place via
`macros/set_summed_pot*.C`, so a config-only digest would have missed it). A cache whose key
is absent or does not match is discarded and rebuilt automatically, with the reason and both
keys printed. `XSEC_FORCE_REBUILD` still works as a manual override, but correctness no
longer depends on remembering it, so the 27 drivers that omitted it are safe by construction.

Verified end to end: an unkeyed cache invalidates ("carries no normalisation key"); an
unchanged run reuses ("normalisation key v2-distinct-sum-pot-bc6edad00c9530dd"); perturbing a
single beam-on trigger count invalidates with both keys shown
(`bc6edad00c9530dd` → `882b038b0266fcd2`); the override still forces a rebuild.

Bump `NORM_SCHEME_TAG` in `SystematicsCalculator.hh` whenever the normalisation changes.

> **Found while testing this:** the Makefile had `-include $(DEPS)` *before* its first real
> target, so a rule from a `.d` file became the default goal — a bare `make` built one object
> file and nothing else. Every plain `make` in this repo was a near-no-op that silently left
> stale binaries, which is indistinguishable from a successful build and compounds exactly the
> failure mode this item is about. (It is why my first attempt to test the cache key appeared
> to do nothing: the binary was 9 minutes older than the code.) Fixed with an explicit
> `.DEFAULT_GOAL := all`.

`SystematicsCalculator.cxx:166-178`: the POT-summed `total_` subfolder is reused whenever it
exists, unless the env var is set. There is no hash or version stamp tying the cache to the
normalisation code that produced it, so after the POT-convention fix a stale cache is
silently reused and yields old numbers with no warning. This is exactly the failure mode
that produced the "tune is 2× low" episode.

**27 of the 60 shell drivers call `UnfolderNuMI` without setting it**, including
`run_fhc5.sh`, `run_rhcfull.sh`, `run_comb.sh`, `run_full_chain.sh`, `reunfold_all.sh`,
`w_batch_fhc.sh`, `w_batch_rhc_comb.sh`.

Fix: stamp the cache with a normalisation-version key and invalidate automatically;
opt-in-by-env-var is not a safe default for a correctness-critical cache.

---

## 5. Missing `Flux` line silently falls back to the FHC value

**Status: FIXED (9892a24).** Sentinel + hard error in NuMI mode; `Flux` added to the 82
configs that lacked one (76 FHC, 5 RHC, 1 comb); all 151 verified numerically.

`CrossSectionExtractor.hh:190`:

```cpp
double flux_per_pot_ = 6.81159e-10;   // FHC
```

If a config has no `Flux` line the FHC normalisation is used with no warning. Five **RHC**
configs lack the line — `ccpi_xsec_config_numi_{costhmu,costhpi,pmu,ppi,thmupi}_rhc.txt` —
which would put those results 5.4% low.

Verified they are currently **unreferenced** by any script or macro, so nothing published is
affected. The eight detVar configs also lack it, but they are all FHC, so the default is
accidentally correct.

Fix: initialise to a sentinel and fail loudly when a NuMI config omits `Flux`. Relying on a
default that is right for one beam out of three is the same class of bug as #2.

The live configs are otherwise consistent: every config's flux matches its beam
(23× FHC, 12× comb, 12× RHC; the `6.608653` / `6.60865` pair differs by 5×10⁻⁷ — formatting
only, not a real discrepancy).

---

## 6. detVar normalisation is inconsistent with the POT fix and depends on hand-written POT

**Severity: medium — currently correct, structurally fragile.**

The recent distinct-sum POT fix was applied to the reweightable-MC path (`:1000`) and the
altCV path (`:886`), but **not** the detVar path, which still divides per file:

```cpp
temp_scale_factor = total_bnb_data_pot_ / file_pot;   // :894
```

The combined config has exactly **2 files per detVar type** (FHC+RHC) which are accumulated,
so the result is `D·(h₁/p₁ + h₂/p₂)` — correct only because each file's `summed_pot` was
hand-set to encode its mode's exposure (comment at `:837-843`). That hand-setting is
invisible from the configs and lives only inside the ROOT files.

Related: `macros/set_summed_pot_v2.C` is what *creates* the duplicate-POT convention — it
writes the run total into every file. The dual convention is therefore self-inflicted, not
inherited from production. Worth stating plainly in the note, since the distinct-sum fix is
otherwise hard to justify to a reader.

I checked the dedup tolerance (`1e-6` relative, `SystematicsCalculator.cxx:511`) against
every multi-file group in every live config: **no collisions**, and the distinct-sum
reproduces the run totals exactly (both conventions give 5.5185e21 for RHC run 3). But the
tolerance is load-bearing — two genuinely-different files whose POT agreed to 1e-6 would be
silently collapsed and the denominator under-counted. Evenly-split grid production makes
that plausible. Consider keying on file identity rather than POT value.

---

## 7. Smaller items

- **Build-sequence axes jump.** `dsigma_build.C:84` and `dsigma_build_1p.C:85` use
  `1.45*ymax`; `dsigma_current.C:82` and `dsigma_ccpi1p.C:83` use `1.35*ymax`. The build
  slides were specifically designed so the axes never move — but the final montage they
  build up to uses a 7% different y-range, so it jumps on the last transition.
- **39 debug `std::cout` statements** left in `SystematicsCalculator.cxx`
  (`"sys test after processing universes 5"` and similar), printed on every run.
- **`h_A_C.Write()` at `UnfolderNuMI.C:1150`** still fails with "the current directory is
  not associated with a file" on every unfold. Harmless — the new closure-file dump added
  today is the working one — but it prints an error that looks alarming in logs.
- **92 of 243 config files hard-code `/home/t2k/nowak/...`** absolute paths, so the repo is
  not relocatable.
- **55 of 60 shell drivers run without `set -e`.** Deliberate in the batch runners (one
  failed observable shouldn't kill the batch) but they also don't aggregate failures, so a
  partial failure is easy to miss — several already write a `.status` file, which is the
  right pattern to apply consistently.

---

## What I checked and found healthy

- Per-run POT scale factors are **identical across all four macros** that duplicate them
  (`0.14101 / 0.05085 / 0.07323 / 0.11560`, `0.06728 / 0.04478 / 0.08847`, `0.09066`) — the
  copy-paste has not drifted. Still worth centralising.
- The fiducial volume and argon target count in `total_xsec_counting.C:14-19` match
  `FiducialVolume.hh:36-50` and `Constants.hh` exactly, including the dead-slab subtraction
  and the deliberate per-nucleus (not per-nucleon) convention.
- Combined bookkeeping is self-consistent: `N_sig_gen(comb) = 5557.2 + 6114.1 = 11671.3`,
  EXT `12.0898 = 5.9313 + 6.1584`, efficiency stable at 17.5% across all three modes.
- The Run5 νe-contaminated detVar quarantine is holding — zero Run5 detVar lines in any of
  the 24 live `file_properties` files.
- The duplicate FV block in `Constants.hh:17-27` is commented out *and* numerically
  identical to the active one, so it is dead but not dangerous.
- `CC1mu1pi1p` inherits cleanly; its truth (`>0.3 GeV/c` true proton) and reco
  (`>0.3 GeV/c` from proton-hypothesis KE) thresholds are consistent, and the tightened
  `PROTON_LLR_CUT = 0.05` is documented with its cut-scan justification.
- The cut-and-count macro is honest about being a tautology under blinding
  (`σ = N_sig_gen/conv` by construction) and says so in its header.

---

## Status

**Fixed:** #1 (as corrected), #2, #3, #5 — commit `9892a24`; #4 plus the Makefile
default-goal bug found while testing it.
**Outstanding:** #6 (detVar POT) and the #7 items.

## Suggested order of work

1. Point the four supporting macros at the per-run EXT files (#1) — the only item that
   changes numbers already written up.
2. Add θ to `make_fte.C` (#3) — cheap, removes a reproducibility hole.
3. Make `Flux` mandatory (#5) and fix the `comb` dirt ternary (#2) — both are one-line
   guards against silent wrong answers.
4. Stamp the universe cache (#4).
5. Cosmetics: debug prints, the `h_A_C` write error, the y-range factor.
