# How to Run the Full xsec_analyzer Chain

The framework runs as a 3-stage pipeline. Each stage has a self-documenting
wrapper script that sources `setup_xsec_analyzer.sh`, builds the needed binary
if missing, runs preflight checks, and prints what to run next.

**Prerequisite:** ROOT must already be in your `PATH` before you start (the
scripts error out otherwise).

## Full chain (run from `xsec_analyzer/`)

```bash
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer

# Stage 1 — ProcessNTuples: apply CC1mu1piXp selection to raw PeLEE ntuples
./run_process_ntuples.sh              # add -j 4 to process 4 files in parallel

# Stage 2 — UniverseMaker: build systematic-universe histograms
./run_universe_maker.sh

# Stage 3 — Unfolder: extract unfolded dsigma/dpmu with full systematics
./run_unfolder.sh
```

## What each stage does and its key I/O

| Stage | Script | Reads | Writes |
|-------|--------|-------|--------|
| 1 | `run_process_ntuples.sh` | `configs/files_to_process_numi.txt` | `xsec-ana-*.root` -> `OUT_DIR` (default `/data/uboone/processed`) |
| 2 | `run_universe_maker.sh` | `configs/file_properties_numi.txt` + `configs/ccpi_pmu_bin_config.txt` | `ccpi_Run1_pmu_univmake.root` |
| 3 | `run_unfolder.sh` | `configs/ccpi_xsec_config_numi.txt` + `configs/ccpi_pmu_slice_config.txt` (+ univ file) | `unfold_output/ccpi_Run1_pmu_xsec.root` + diagnostic PDFs |

## Common overrides (env vars)

```bash
# Stage 1: where processed files go, which selection, which input list
OUT_DIR=/somewhere SELECTION=CC1mu1piXp ./run_process_ntuples.sh -j 4

# Stage 2: pick a different bin scheme / output
BIN_CONFIG=configs/foo.txt OUT=/path/univ.root ./run_universe_maker.sh

# Stage 3: skip the diagnostic-plot step, or swap configs
PLOTS=0 ./run_unfolder.sh
XSEC=configs/foo.txt SLICE=configs/bar.txt OUT=/path/out.root ./run_unfolder.sh
```

## Notes / gotchas

- **Order matters and is enforced:** Stage 2 checks that every processed input
  from the FPM file exists (else tells you to run Stage 1); Stage 3 checks the
  `UnivFile` referenced in the xsec config exists (else tells you to run Stage 2).
- Stage 3's `Unfolder` produces the measurement ROOT file; `UnfolderNuMI` (run
  after, unless `PLOTS=0`) produces the diagnostic step1-4 PDFs into
  `unfold_output/`. The plots step is non-fatal if it fails.
- There is also `run_nuwro_observables.sh` for NuWro generator-level observables
  — a separate generator-comparison path, not part of the core measurement chain.
- There is no single master "do-everything" wrapper — run the three scripts in
  sequence.
