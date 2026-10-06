# Plan: rewrite the analysis note as the current state, history to appendices

Written 2026-08-29 against `analysis_note.tex` @ 6061 lines / 91 pages (commit 465997b).
Line numbers below are from that revision.

---

## 0. The finding that determines the scope

**This is a re-derivation, not a restructure.**

None of the 18 current integrated cross sections appear anywhere in the note:

| | current values | occurrences in note |
|---|---|---|
| FHC5 | 0.9179 0.7946 0.8240 0.9304 0.9344 0.9231 | 0 of 6 |
| RHCFULL | 0.6279 0.7113 0.7788 0.8782 0.9173 0.8496 | 0 of 6 |
| COMB | 0.7797 0.8120 0.8680 0.9381 0.9526 0.8747 | 0 of 6 |

The results layer predates the beam-frame re-run and, for COMB, the detector-systematic
fix (80d613d). The `ppi2bin` table was already found stale for **all three**
configurations, not just COMB, and has been refreshed (465997b); the rest has not.

Consequence: moving sections around does not produce a correct note. Every results
table, every quoted sigma, and every claim derived from them has to be re-extracted
before or during the rewrite. Restructuring a document whose numbers are wrong is the
more expensive order of operations, because the prose has to be revisited twice.

**Recommended order: re-derive the numbers first, restructure second.**

---

## 1. Structural diagnosis

### 1.1 ~~Five sections sit after the Summary~~ — **WRONG, corrected 2026-08-30**

The note already has `\appendix` at :5367, so those five sections are appendices A-E,
not misplaced body sections. They are correctly placed. What follows in this
subsection was written before I checked for `\appendix` and is retained only to show
what was actually verified.

The real issue is narrower and survives: **the configuration-closure result
(`sec:combtruth`) sits inside Appendix C**, when it is a current, load-bearing finding
-- it is what identified and validated the combined detector-systematic fix. That one
subsection should be promoted to the body; A, B, D and E stay where they are.

The original (incorrect) diagnosis follows.

### 1.1-old Five sections after the Summary — 763 lines, 13% of the note

`Summary` is at :5220, and then the document continues:

| line | section | lines |
|---|---|---|
| 5298 | Identity closure: unfolding the prediction with itself | 84 |
| 5382 | Data/MC Normalisation | 90 |
| 5472 | NuWro Fake-Data Closure and Cross-Pipeline Validation | 344 |
| 5816 | Observable Binning | 34 |
| 5850 | Selection Cut Parameter Summary | 211 |

A reader who stops at the Summary misses the configuration-closure result, which is
currently the most consequential open item in the analysis. This is accretion: material
was appended as it was produced rather than placed.

### 1.2 The issues table presents history as status

`Known Issues and Required Follow-up` (:4797, 97 lines) has 19 rows:
**12 Resolved, 1 Understood, 2 Open**, plus `Not started`, `Data files needed`,
`Blocked`, `Dead`. Two-thirds of the table is a changelog. The two genuinely open rows
are buried among them.

### 1.3 ~920 lines of investigation narrative in the body

These read as lab notebook rather than result:

| line | subsection | lines | character |
|---|---|---|---|
| 2207 | The pion-momentum bias, and why $p_\pi$ is worst-conditioned | 659 | physics result + long derivation |
| 3244 | Survey of pion momentum estimators without a magnetic field | 163 | four approaches measured, all fail |
| 3145 | What limits $W_{\pi p}$, and what can be done about it | 99 | "It does not work, for a reason distinct from..." |
| 3046 | The $W$ binnings are finer than the resolution | 99 | superseded binning |

The *conclusions* belong in the body (they set real floors on the measurement). The
narrative of how they were reached does not.

### 1.4 Resolved review items as body subsections

`Data/MC Normalisation` (:5382) contains, as subsections: unfolding covariance
(~26% effect, resolved), EXT subtraction for fake data (resolved), opening-angle
binning "review item C4" (resolved), generator normalisation (resolved), and
**software trigger "verified to be a no-op"** (:5447) — a refuted hypothesis, retained
in full. All are history.

---

## 2. Proposed structure

**Body — the analysis as it now stands.** A reader should be able to go start to finish
and get the measurement, the method, and the honest limitations, without encountering a
superseded number or an investigation that went nowhere.

```
1  Introduction
2  Theory and prior measurements            (unchanged)
3  Detector                                 (unchanged)
4  Beam and flux                            (unchanged; beam-frame subsection is current)
5  Samples, signal definition, fiducial volume
6  Selection and cut-flow
7  Unfolding method                         (Wiener-SVD, A_C, response, resolution)
8  Results: FHC, RHC, combined              <-- FULLY RE-DERIVED
9  Proton-tagged subsample                  <-- FULLY RE-DERIVED
10 Systematics                              <-- re-derived; detVar section rewritten
11 Validation                               (identity closure, fake-data closure,
                                             configuration closure, cross-pipeline)
12 Limitations and open questions           <-- new, short, honest
13 Summary
```

**Appendices — everything that is true but historical.**

- **A. Superseded methods and refuted hypotheses.** Software-trigger no-op; the
  ~26% unfolding-covariance correction; EXT cancellation for fake data; generator
  normalisation; opening-angle C4; the five-bin and six-bin $W$ schemes; the
  fine-binned $p_\pi$ scheme; the detVar-files hypothesis and the missing-`apply_ac`
  hypothesis (both tested and refuted, 2026-08-28).
- **B. Investigations that set limits.** Pion-momentum estimator survey; what limits
  $W_{\pi p}$; the $p_\pi$ conditioning derivation. Keep the numbers, compress the
  narrative; the body cites the conclusion and points here.
- **C. Review responses.** The existing `NOTE_REVIEW_PLAN.md` material and the issues
  table's 13 resolved rows, as a dated changelog.
- **D. Reference tables.** Selection cut parameters (:5850, 211 lines), observable
  binning (:5816). These are lookup material, not narrative — they belong at the back.

---

## 3. Section-by-section disposition

| line | section | lines | disposition |
|---|---|---|---|
| 144-911 | Intro, theory, detector, beam/flux | 768 | **keep**, light edit. Beam-direction subsection (:718) is current. |
| 912-1228 | Samples, signal definition, FV | 317 | **keep**, verify POT table against current file lists |
| 1229-1702 | Selection, cut-flow, performance | 474 | **keep**; regenerate the data/MC and N-1 figures |
| 1703-2206 | Unfolding chain, response, resolution | 504 | **keep** — method, not results |
| 2207-2865 | $p_\pi$ bias (659L) | 659 | **split**: ~150L conclusion to body §7, remainder to **App. B** |
| 2866-3045 | Model comparison, proton-tagged $\chi^2$/$A_C$ | 180 | **keep**, re-derive all numbers |
| 3046-3243 | $W$ binnings, what limits $W_{\pi p}$ | 198 | **App. B**, cite conclusion in body |
| 3244-3406 | Pion momentum estimator survey | 163 | **App. B** |
| 3407-3475 | Cross-check methods | 69 | **keep** → body §11 |
| 3476-3716 | Full-exposure results | 241 | **rewrite entirely** — all numbers superseded |
| 3717-4096 | Angular variable, $\cos\theta_\mu$ vs $\theta_\mu$ | 380 | **keep** the $A_C$ rank argument; re-derive numbers. Note `fhc-thetamu-ac-is-rank1` finding must be reflected. |
| 4097-4321 | Sidebands, control regions | 225 | **keep**, re-derive |
| 4322-4765 | Proton-tagged subsample | 444 | **rewrite results**, keep method |
| 4766-4796 | Systematics | 31 | **expand** — currently thin for a systematics section; fold in the breakdown from :4971 |
| 4797-4893 | Known Issues table | 97 | **split**: 2 open rows → body §12; 13 resolved → **App. C** |
| 4894-5219 | Pipeline, binning opt, syst breakdown, flux anatomy, multi-run, nue-sibling | 326 | **redistribute**: syst breakdown → §10; flux anatomy → §10; multi-run closure → §11; nue comparison → §11; pipeline → **App. D** |
| 5220-5297 | Summary | 78 | **rewrite last** |
| 5298-5381 | Identity closure | 84 | **move up** → body §11 |
| 5382-5471 | Data/MC normalisation | 90 | **App. A** (all five subsections are resolved history) |
| 5472-5815 | NuWro closure, extraction bias, configuration closure | 344 | **move up** → body §11. Configuration closure (:5638) is current and important. |
| 5816-5849 | Observable binning | 34 | **App. D** |
| 5850-6061 | Selection cut parameters | 211 | **App. D** |

Net effect: body shrinks by roughly 1500 lines of narrative, appendices absorb it, and
the reading order stops doubling back.

---

## 4. What must be re-derived, not moved

1. **All 18 inclusive $\sigma_\mathrm{int}$** and the per-bin $d\sigma$ tables.
   Sources exist: `logs/systdump/<cfg>_<obs>.dump` and
   `/data/uboone/processed/closure_hists_xsec_<TAG>_<obs>.root`, both regenerated
   2026-08-29.
2. **`tab:sigint_all`** — the observable-vs-configuration sigma table.
3. **`tab:dagostini`** (:3468) — quotes WSVD means FHC 0.923 / RHC 0.836 / COMB 0.852.
   The current inclusive-family means are 0.887 / 0.794 / 0.748 pre-fix; COMB changes
   again post-fix. Requires a D'Agostini re-run to restate honestly.
4. **All proton-tagged results** — the seven observables were rebuilt 2026-08-27/28.
5. **Systematic breakdown** (:4971) — `detVar_total` for COMB changed 42.3% -> 24.6%,
   so every COMB systematic total moves.
6. **`ppi2bin` table** — DONE (465997b), as the worked example of what the rest needs.

---

## 5. Content that needs writing, not editing

- **§12 Limitations and open questions.** Does not exist. Should state plainly:
  - the unfolded-to-truth ratio is **1.13 in all three configurations** — a real,
    uniform ~13% Wiener-SVD normalisation offset on fake data, the same order as the
    physics differences being measured, documented but **not explained**;
  - the $A_C$ normalisation inconsistency ($A_C\cdot$truth $=0.79\times$truth while
    unfolded $=1.16\times$truth — both cannot hold if $\hat{x}=A_C x_\mathrm{true}$);
  - `costhmu`'s corrected COMB detVar sits slightly below both inputs (23.7% vs
    26.1/32.5) where the other five are bracketed;
  - FHC $\theta_\mu$ $A_C$ is rank-1, so "$\theta$ fixes $\cos\theta$" is by
    construction — already withdrawn, must not creep back.
- **§10 detVar subsection.** The per-mode double-counting bug and its fix, currently
  only in `sec:combtruth`, belongs in the systematics section proper.

---

## 6. Risks and prerequisites

- ~~Three figures still have no producer~~ **RESOLVED 2026-08-30.** `fw_bkgsub_*`,
  `fw_xsec_*` and `dsigma_build_1p_*` are referenced by **no** `.tex` file, note or
  slides. They are not producer-less stale content; they are dead files, and mislead
  nobody. My earlier claim that they "show superseded COMB numbers" was an inference
  from the `ppi2bin` precedent and was wrong.

  The audit that established this also found **82 of 237 figure basenames (35%) are
  orphaned** and **0 references are broken** (an earlier count of 16 "missing" was an
  artefact of listing only `.pdf` when many figures are `.png`). See
  `ORPHAN_FIGURES.txt`. Orphan cleanup is cosmetic and should not gate the rewrite.
- **D'Agostini numbers** need a re-run before `tab:dagostini` can be restated.
- Do not restructure and re-derive in the same pass. Re-derive, verify, commit; then
  restructure.

---

## 7. Suggested sequence

1. Re-extract all results into a machine-readable table; diff against the note; publish
   the diff. (Half a day. Largely scripted — the closure files already exist.)
2. ~~Check the three producer-less figure families~~ DONE — they are unreferenced;
   no action needed. See `ORPHAN_FIGURES.txt`.
3. Rewrite §8, §9, §10 from the re-extracted numbers.
4. Write §12 from the open items above.
5. Promote `sec:combtruth` from Appendix C into the body. Appendices A-E already
   exist and are correctly structured -- no wholesale move is needed (see §1.1).
6. Rewrite §1 and §13 last, once the body is settled.
7. Rebuild, check for undefined references, and diff the PDF page count as a sanity
   check on what moved.
