# Response to noteReview.txt and noteReview2.txt

Every claim below was checked against the source and the code before being accepted.
Verdicts are mine; where a reviewer is wrong or partly wrong I say so.

## A. Confirmed errors — fix now, independent of the re-run

### A1. N_Ar exponent typo  (both reviews)  **CORRECT — real error**
`analysis_note.tex:892` (Eq. 2) gives `N_Ar ~ 8.71e28`; the table at :1931 and
Sec. 9.1 at :3674 give `8.710e29`.

Verified against the code, not just internally. `CrossSectionExtractor.hh:517` calls
`num_Ar_targets_in_FV(FV, 675.1, 775.1)` with `FV = {10, 246, -101, 101, 10, 986}`
from `Constants.hh`, rho=1.3836 g/cm3, A=39.948, N_A=6.02214076e23:

    full FV                      V = 4.6528e7 cm3   N_Ar = 9.7046e29
    minus the 675.1-775.1 slab   V = 4.1761e7 cm3   N_Ar = 8.7103e29   <- used

So `8.710e29` is right and Eq. 2 is a one-character typo. Review 2's reasoning is
also right: `conv = 4.527e3` only closes with 1e29.

Fix: change the exponent in Eq. 2. Add the numerical derivation (FV box, dead-region
subtraction, density, molar mass) so the number is checkable on the page.

### A2. Off-axis angle quoted two ways  (review 2)  **CORRECT — and my fault**
:486 and :673 still say `110 mrad`; the new Sec. 2.2 says `134.6 mrad = 7.71 deg`
and footnotes that it supersedes the earlier value. I introduced the contradiction by
adding the superseding footnote instead of correcting the original.

The new value is the defensible one: it follows from the detector position in beam
coordinates (5502, 7259, 67270) cm, the same position the official NuMI beamline
geometry weights use. 7.7 deg is consistent with the ~8 deg usually quoted for
MicroBooNE; 110 mrad = 6.3 deg is not.

Fix: correct :486 and :673 to 135 mrad, keep one short footnote recording the change.

### A3. Uncertainty text contradicts Table 31  (review 2)  **CORRECT**
:3656 says FHC p_mu total ~27%, flux 23%, detector 7%, xsec 8%, data stat 8%.
Table 31 says total 32.9, flux 24.6, detector 12.3, xsec 9.2, data stat 11.6.

Table 31 is generated from the systematics dump and is authoritative; the prose is
stale. p_mu is frame-independent and its univmake is being reused, so these numbers
will not move in the re-run - safe to fix now.

Fix: rewrite the sentence from Table 31.

### A4. Response matrix defined two incompatible ways  (review 1, pt 3)  **CORRECT**
:1839-1841 says each column is normalised to unit sum, so
`A_ij = P(reco i | true j)`, and then that `eps_j = sum_i A_ij`. Both cannot hold: if
columns sum to unity then eps_j = 1 identically.

The reviewer's diagnosis is right. Two distinct objects are in play and the note
conflates them:

    M_ij = P(reco i | true j, selected)      sum_i M_ij = 1     (the plotted matrix)
    S_ij = eps_j * M_ij                      sum_i S_ij = eps_j (used in unfolding)

Resolved against the code. `UnfolderNuMI.C:259` passes
`syst.get_cv_smearceptance_matrix()` to the unfolder, and :405 computes the efficiency
as "column sums of smearceptance". `WienerSVDUnfolder.cxx` takes that matrix as
`smearcept` throughout. The `[DIAGDUMP]` block at :391 prints a "column-normalised
response diagonal", i.e. it normalises for the diagnostic - proof the stored matrix is
not column-normalised.

So the unfolder uses S, and of the two sentences in the note the *efficiency* one is
right and the *normalisation* one is wrong.

Fix: delete "Each column is normalised to unit sum so that A_ij = P(reco i | true j)",
define S_ij = eps_j M_ij with column sum eps_j as the matrix used in the unfolding, and
note that Fig. 18 plots the column-normalised M for legibility.

### A5. "Observed yield" in a blind analysis  (review 1, pt 8)  **CORRECT — important**
Abstract :123: "The observed yield is ~0.75x the GENIE v3 prediction". The analysis
runs entirely on Poisson-thrown fake data; no real beam-on data is unblinded. A reader
will take this as a measured 25% deficit.

Fix: "For the current fake-data extraction, the injected central-value sample
corresponds to ~0.75x the standalone GENIE v3 prediction in this phase space."
Sweep the note for every other "observed"/"data deficit" phrasing.

### A6. Sec. 11 contradicts the systematics claim  (review 1, pt 7)  **CORRECT**
:3651 says every result carries full statistical and systematic uncertainties;
:4740 says "Full systematic uncertainty propagation has not yet been performed".

Both are true of different things - the framework propagates the covariance, the
custom pipeline imports it rather than generating native multisims - but the note
never says so. Rewrite Sec. 11 to state exactly that, and drop the blanket sentence.

### A7. Minor  (review 2)  **CORRECT, cosmetic**
- RHC POT written `1.108e21` at :98 and `11.082e20` at :694/:830. Pick one.
- "Resolved" rows in the Known Issues table whose impact text still reads as pending
  (:4788 "re-enable for real data", :4794 "full custom propagation pending"). Split
  the column into status and follow-up.
- A changelog note inside a table cell ("Changed 2026-08-17, see 8.8") - move to the
  caption or drop.

## B. Correct, but resolved by the re-run rather than by editing

The reviews object at length to superseded material sitting beside current material -
Sec. 8.10 "Withdrawn pending re-extraction" followed by pre-fix arguments, Sec. 9.2's
retracted peak-sharpness explanation, Sec. 10.2's pre-fix A_C warning contradicted by
its own table note.

They are right that it reads as contradictory. But those retraction boxes are
deliberate: they mark results that the beam-frame correction invalidated, pending the
re-run now in progress. Deleting them now would leave the note asserting numbers known
to be wrong. Deleting them *before* replacements exist is the worse failure.

Plan: when the re-run lands, replace each box with the final post-beta result and move
the historical material to a "Superseded studies" appendix. Not before.

This covers review 1 pts 4, 6 (partly) and review 2's "ghost sections" and
"narrative whiplash" items.

## C. Where I disagree, or would go further

### C1. "Do not present five-bin p_pi next to the final result" - agree
Reviewer 1's framing advice is sound and I would adopt it, with one addition: the same
argument now applies to *four* more observables. The TKI re-derivation on the
beta-corrected data (configs/ccpi1p_TKI_binning.README) shows six bins reach only
25-34% migration diagonals; delta_pT, delta_phiT and p_n support two bins, and
delta_alphaT meets the criterion at no binning and is retained by explicit decision.

So the note needs a single, general statement of the principle - "these are
response-limited measurements; bin counts are set by the migration matrix, not by
model separation" - covering p_pi, W_pipr and all four TKI observables, rather than
arguing it once per observable.

### C2. "The upper p_pi bin is an integral, not a shape point" - agree, and it
generalises. Same wording should apply to the upper TKI bins.

### C3. p_pi threshold vs the [0.113, 0.175] bin  (review 2)  **partly wrong**
The reviewer reads the [0.113, 0.175] row as an added analysis bin that Sec. 14.1
fails to validate. It is not: :1063 and the table at :1078 are a *comparison with the
BNB analysis threshold*, showing why we cut at 0.175 rather than 0.113 - the resolution
is 56-73% below 0.175. No such bin is measured. No fix needed beyond making the
table caption say plainly that these rows are the rejected region.

### C4. Compression advice - agree in principle, defer
Both reviews want the p_pi saturation story and the A_C caveat consolidated. Right,
but this is a large structural edit that will conflict with the re-run rewrite. Do it
in the same pass as B.

## Order of work

Now, independent of the re-run:  A1, A2, A3, A5, A7  (factual and blinding fixes)
Now, also settled:               A4  (unfolder uses the smearceptance matrix)
Needs the re-run:                A6, B, C1, C2, C4
