#!/bin/bash
# run_nuwro_pred.sh — run inside the SL7 container:
#   apptainer exec -B /cvmfs -B /home -B /data <img> bash run_nuwro_pred.sh
# Generates the NuWro NuMI numu-on-Ar sample, applies the CC1mu1piXp signal,
# and writes the differential-cross-section prediction ROOT file.
HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
cd "$HERE"
echo "[start] $(pwd)"

# --- NuWro environment ------------------------------------------------------
# Do NOT run these under `set -e`: the UPS setup scripts return non-zero from
# internal probes even on success, which would abort the whole script before
# any output. Enable strict mode only once the environment is ready.
set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
if [[ -z "$NUWRO_FQ_DIR" ]]; then
  echo "[fatal] NuWro UPS setup failed (NUWRO_FQ_DIR unset)"; exit 1
fi
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
mkdir -p "$HOME/generators/.nuwro_run/bin"
ln -sfn "$NUWRO_FQ_DIR/bin/nuwro" "$HOME/generators/.nuwro_run/bin/nuwro"
ln -sfn "$NUWRO/data" "$HOME/generators/.nuwro_run/data"
export PATH="$HOME/generators/.nuwro_run/bin:$PATH"

set -e   # environment ready; strict from here on
echo "[nuwro] $(command -v nuwro)"
echo "[root]  $(command -v root-config) $(root-config --version)"

# --- 1. generate ------------------------------------------------------------
OUT_NUWRO=nuwro/out_numu_numi.root
if [[ -s "$OUT_NUWRO" && "${SKIP_GEN:-0}" == "1" ]]; then
  echo "[gen] reusing existing $OUT_NUWRO"
else
  echo "[gen] running NuWro (this is the slow step) ..."
  nuwro -i nuwro/params_numu_numi.txt -o "$OUT_NUWRO"
fi

# --- 2. build analyzer ------------------------------------------------------
# The `event` class dictionary lives in event1.so (next to the nuwro binary),
# built by NuWro with rootcint; link it so event::Class()/vtable resolve.
echo "[build] nuwro_cc1pi ..."
EVENT_SO="$NUWRO_FQ_DIR/bin/event1.so"
[[ -f "$EVENT_SO" ]] || { echo "[fatal] event1.so not found at $EVENT_SO"; exit 1; }
g++ analysis/nuwro_cc1pi.cc -o analysis/nuwro_cc1pi \
    -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$EVENT_SO"

# --- 3. analyze -------------------------------------------------------------
echo "[analyze] applying CC1mu1piXp signal ..."
./analysis/nuwro_cc1pi "$OUT_NUWRO" out/nuwro_cc1pi_pred.root

echo "[done] prediction -> out/nuwro_cc1pi_pred.root"
