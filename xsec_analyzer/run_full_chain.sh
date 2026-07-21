#!/usr/bin/env bash
# run_full_chain.sh
# One-shot driver for the complete NuMI CC1mu1piXp analysis chain, end to end:
#
#   Stage 1  ProcessNTuples   raw PeLEE ntuples  -> processed stv_tree files
#   Stage 2  univmake         processed files    -> per-observable universe file
#   Stage 3  Unfolder(+NuMI)  universe file       -> per-observable cross section
#
# Stage 1 runs ONCE and produces all processed files (they are
# observable-independent). Stages 2 and 3 run once per observable. This is the
# run required after the momentum-branch corrections, because the corrected
# reco branches (candidate_muon_mom_reco, the range-momentum pion estimator) do
# not exist in the previously processed files.
#
# Everything is logged, per stage and per observable, under logs/full_chain_<ts>/.
# Each stage is preflighted; the run aborts on the first hard failure.
#
# Usage:
#   ./run_full_chain.sh                     # full chain, all observables
#   ./run_full_chain.sh -j 6                # 6 parallel ProcessNTuples jobs
#   ./run_full_chain.sh --start-at 2        # skip stage 1 (reuse processed files)
#   ./run_full_chain.sh --obs pmu,ppi       # only these observables in stages 2-3
#   ./run_full_chain.sh --no-plots          # skip UnfolderNuMI diagnostic plots
#
# Env overrides honoured by the underlying scripts also work here
# (OUT_DIR for stage 1, etc.).

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"

# ── Options ────────────────────────────────────────────────────────────────
JOBS=4
START_AT=1
OBS_LIST=(pmu ppi costhmu costhpi thmupi)
PLOTS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)     JOBS="$2"; shift 2 ;;
    --start-at)    START_AT="$2"; shift 2 ;;
    --obs)         IFS=',' read -ra OBS_LIST <<< "$2"; shift 2 ;;
    --no-plots)    PLOTS=0; shift ;;
    -h|--help)     sed -n '2,30p' "$0" | sed 's/^# \?//'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# TS is passed in from the caller (the script cannot call date under some
# harnesses); fall back to a fixed label if unset.
TS="${FULL_CHAIN_TS:-run}"
LOGDIR="logs/full_chain_${TS}"
mkdir -p "$LOGDIR"

# ── Environment ────────────────────────────────────────────────────────────
set +u; source ./setup_xsec_analyzer.sh; set -u
command -v root &>/dev/null || { echo "ERROR: ROOT not in PATH" >&2; exit 1; }

banner() { echo; echo "======== $* ========"; }
fail()   { echo "FATAL: $*" >&2; exit 1; }

echo "NuMI CC1mu1piXp full chain"
echo "  log dir     : $LOGDIR"
echo "  stage 1 jobs: $JOBS"
echo "  start at    : stage $START_AT"
echo "  observables : ${OBS_LIST[*]}"
echo "  plots       : $PLOTS"

# ── Stage 1: ProcessNTuples ─────────────────────────────────────────────────
if [[ "$START_AT" -le 1 ]]; then
  banner "STAGE 1  ProcessNTuples  (-j $JOBS)"
  if ./run_process_ntuples.sh -j "$JOBS" > "$LOGDIR/stage1_process.log" 2>&1; then
    echo "  stage 1 OK   -> $LOGDIR/stage1_process.log"
  else
    tail -20 "$LOGDIR/stage1_process.log" >&2
    fail "stage 1 (ProcessNTuples) failed; see $LOGDIR/stage1_process.log"
  fi
else
  echo "  (skipping stage 1)"
fi

# ── Stage 2: univmake, per observable ───────────────────────────────────────
if [[ "$START_AT" -le 2 ]]; then
  banner "STAGE 2  univmake  (${#OBS_LIST[@]} observables)"
  for o in "${OBS_LIST[@]}"; do
    bin_cfg="configs/ccpi_${o}_bin_config.txt"
    out="/data/uboone/processed/ccpi_Run1_${o}_univmake.root"
    [[ -f "$bin_cfg" ]] || fail "missing bin config: $bin_cfg"
    echo "  [$o] univmake ..."
    if BIN_CONFIG="$bin_cfg" OUT="$out" \
         ./run_universe_maker.sh > "$LOGDIR/stage2_univ_${o}.log" 2>&1; then
      echo "        OK   -> $out"
    else
      tail -20 "$LOGDIR/stage2_univ_${o}.log" >&2
      fail "stage 2 univmake for $o failed; see $LOGDIR/stage2_univ_${o}.log"
    fi
  done
else
  echo "  (skipping stage 2)"
fi

# ── Stage 3: Unfolder (+ UnfolderNuMI), per observable ──────────────────────
if [[ "$START_AT" -le 3 ]]; then
  banner "STAGE 3  Unfolder  (${#OBS_LIST[@]} observables)"
  for o in "${OBS_LIST[@]}"; do
    xsec_cfg="configs/ccpi_xsec_config_numi_${o}.txt"
    slice_cfg="configs/ccpi_${o}_slice_config.txt"
    out="unfold_output/ccpi_Run1_${o}_xsec.root"
    [[ -f "$xsec_cfg" ]]  || fail "missing xsec config: $xsec_cfg"
    [[ -f "$slice_cfg" ]] || fail "missing slice config: $slice_cfg"
    echo "  [$o] Unfolder ..."
    if XSEC="$xsec_cfg" SLICE="$slice_cfg" OUT="$out" PLOTS="$PLOTS" \
         ./run_unfolder.sh > "$LOGDIR/stage3_unfold_${o}.log" 2>&1; then
      # Surface the zeroed-category warning summary if it fired.
      if grep -q "SYSTEMATICS WARNING SUMMARY" "$LOGDIR/stage3_unfold_${o}.log"; then
        echo "        OK  (with systematics warning -- see log)"
      else
        echo "        OK   -> $out"
      fi
    else
      tail -20 "$LOGDIR/stage3_unfold_${o}.log" >&2
      fail "stage 3 Unfolder for $o failed; see $LOGDIR/stage3_unfold_${o}.log"
    fi
  done
else
  echo "  (skipping stage 3)"
fi

banner "DONE"
echo "Cross-section files:"
for o in "${OBS_LIST[@]}"; do
  f="unfold_output/ccpi_Run1_${o}_xsec.root"
  [[ -f "$f" ]] && echo "  $o  ->  $f ($(du -h "$f" | cut -f1))" || echo "  $o  ->  (missing)"
done
echo
echo "Logs in $LOGDIR"
