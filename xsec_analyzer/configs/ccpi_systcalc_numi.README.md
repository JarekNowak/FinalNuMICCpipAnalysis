# ccpi_systcalc_numi.conf — what this config does and does not cover

**The `.conf` format has no comment syntax.** The parser is a bare
`while (config_file >> name >> type)` loop, so any `#` line is read as a
name/type pair and throws (this is exactly what broke `configs/systcalc.conf`,
whose orphan beamline token block still makes it unparseable). Hence this
separate file.

## Removed: per-category detector variations and dirt normalization

The following entries were removed because they contributed **identically
zero** to the total uncertainty:

    detVarNumupizero MCFullCorrCategory 0.2 2
    detVarNumuOther  MCFullCorrCategory 0.2 3
    detVarNCpizero   MCFullCorrCategory 0.2 4
    detVarNCOther    MCFullCorrCategory 0.2 5
    detVarOutFV      MCFullCorrCategory 0.2 6
    detVarNumu       sum 5 <the above>
    dirtNorm         MCFullCorrCategory 1.0 10

### Why they were zero

`MCFullCorrCategory` matches a true bin by exact string equality against
`"category == N"`. Every CCpi bin config declares its background as a single
`"!CC1mu1piXp_MC_Signal"` true bin, so nothing ever matched and the covariance
came out zero in every bin — silently, with no diagnostic.

### Why restoring the numbers as-is would have been worse

Those category indices are not CC1mu1piXp categories. They map one-to-one,
names and numbers both, onto the **NuMI CC1e** scheme in
`EventCategoriesNuMICC1e.hh`:

| entry | | NuMICC1e |
|---|---|---|
| `detVarNumupizero` 2 | ↔ | `kNuMuCCPi0 = 2` |
| `detVarNumuOther` 3 | ↔ | `kNuMuCCOther = 3` |
| `detVarNCpizero` 4 | ↔ | `kNCPi0 = 4` |
| `detVarNCOther` 5 | ↔ | `kNCOther = 5` |
| `detVarOutFV` 6 | ↔ | `kOOFV = 6` |

Under CC1mu1piXp's scheme (`EventCategoriesXp.hh`) those same numbers mean
`kOOFV`, `kUnknown`, `kNuECC`, `kNC`, `kOther`. Had the match ever worked, a
"νμ CC π⁰" detector variation would have been applied to out-of-FV events.
The config was copied from the CC1e analysis and never re-mapped.

`dirtNorm 10` is wrong twice over: no scheme defines category 10, and dirt is
not an event category at all — it is a *sample type* (`kDirtMC`, two entries in
`file_properties_numi.txt`).

### Verified impact of removal: none

Running `bin/Unfolder` on `ccpi_Run1_costhmu_univmake.root` before and after
this change produced **152 bit-identical histograms, worst relative difference
0**. The only difference is that seven all-zero placeholder histograms are no
longer written.

## What is therefore NOT in the quoted uncertainty

1. **Background detector-variation normalization.** `detVar_total` now covers
   only the eight `DV` sample-based variations (LY, recombination, SCE, wire
   modification). There is no additional normalization term on the background
   prediction by category.
2. **Dirt normalization.** No dirt uncertainty of any kind is included.

Both were absent from every result produced before this change too — removing
the entries did not lose anything that was previously counted. It only stopped
the config from claiming coverage it did not provide.

## Restoring them properly

**Per-category detector variations** require, in order:

1. `categorize_event()` in `CC1mu1piXp.cxx` to actually break out the
   background. It currently returns only `kUnknown`, `kOOFV`, `kNC`, `kNuECC`
   and `kNumuCC_sig`; `kEXT` and `kOther` are commented out, and — critically —
   **real data and the dominant νμ CC in-FV background both fall into
   `kUnknown`**. The νμ CC breakdown lines are commented out.
2. The `MakeConfig` background-bin regression fixed: `MakeConfig.cxx:290` loops
   over `background_index`, but `TutorialBinScheme` sets neither
   `background_index` nor `CATEGORY` (both empty-initialized in
   `BinSchemeBase.hh:46`), so regenerating any bin config currently emits *no*
   background true bins at all.
3. Bin configs regenerated with per-category background true bins.
4. These entries re-added with CC1mu1piXp category numbers.

**Dirt normalization** cannot be expressed in this config format today.
`MCFullCorr` applies its fraction to the *entire* CV prediction, so
`dirtNorm MCFullCorr 1.0` would put 100% on everything. The options are a
sample-selective covariance type, or an `MCFullCorr` scaled to the dirt
fraction of the selected sample (fully correlated across bins, approximate).

## Also missing: NuMI beamline flux systematics

`Configs/systcalc_numi.conf` (capital C, unused by any code path) holds an
intact block of `flux_Horn_2kA`, `flux_Horn1_x_3mm`, `flux_Beam_spot_*` etc.
Those weights are **absent from the processed ntuples** — only
`weight_ppfx_all` is present. `instructions_numi.txt:26` documents
`AddBeamlineGeometryWeights` as the tool that injects them, but no `run_*.sh`
invokes it. Restoring the config block alone would do nothing; that stage has
to go back into the pipeline ahead of `ProcessNTuples`.
