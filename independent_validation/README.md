# Independent validation of the NuMI CC1π extraction

Cross-checks of the `xsec_analyzer` framework extraction against a fully
framework-independent analysis of the same Run-1 NuMI FHC final state
(`/home/t2k/nowak/MicroBooNE/XSecCCPip`, standalone selection + Wiener-SVD /
D'Agostini unfolding). Supports §6.7 of `report/NuMIInternalNoteCC1pi.md`.

## Provenance
`WienerSVD.h`, `sel_run/ccpi_selection.C`, `FV_new.h`, `CCPiConfig.h` are copies
of the standalone analysis at `XSecCCPip/`. `sel_run/ccpi_selection.C` carries two
local edits for the efficiency/purity cross-check (marked in-file): skip the
zero-weight Run2/4/5 MC and NuWro slots, and a `SEL_FRAC` uniform event-fraction
cap (ratios preserved). Set `SEL_FRAC=1` for the full exact run.

## Contents
- `fw_inputs_indep.C` — runs the independent `RunWienerSVD_FW` on the framework's
  own numuMC inputs (algorithm-equivalence check).
- `selfclose_indep.C` — standalone self-closure per observable.
- `sel_run/ccpi_selection.C` — the independent selection, for the eff/purity check.

## Key results (2026-07-25)
- **Wiener-SVD algorithm:** framework `WienerSVDUnfolder.cxx` is identical to the
  independent `RunWienerSVD_FW`; the framework's exact matrices reproduce its
  result bin-for-bin. **No framework unfolding bug.**
- **numuMC self-closure:** x̂/truth = 0.94–1.06 across all five observables × three
  regularisations — the extraction is unbiased once the fake-data sample matches
  the response. The ~1.3× in the detVar-CV closure was the 1.29× POT offset.
- **Selection:** signal definition identical to the framework; efficiency 16.5%
  (framework 16%), purity 54.9% (framework 60%). The old 12.1% report value is stale.
