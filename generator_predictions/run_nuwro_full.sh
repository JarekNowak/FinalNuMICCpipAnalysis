#!/bin/bash
# run_nuwro_full.sh — inside the SL7 container. Produces the flux-combined,
# per-nucleus NuWro CC1mu1piXp prediction from separate numu and numubar runs.
HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
cd "$HERE"
echo "[start] $(pwd)"

set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof
[[ -z "$NUWRO_FQ_DIR" ]] && { echo "[fatal] NuWro setup failed"; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
mkdir -p "$HOME/generators/.nuwro_run/bin"
ln -sfn "$NUWRO_FQ_DIR/bin/nuwro" "$HOME/generators/.nuwro_run/bin/nuwro"
ln -sfn "$NUWRO/data" "$HOME/generators/.nuwro_run/data"
export PATH="$HOME/generators/.nuwro_run/bin:$PATH"
set -e

EVENT_SO="$NUWRO_FQ_DIR/bin/event1.so"
echo "[build] nuwro_cc1pi ..."
g++ analysis/nuwro_cc1pi.cc -o analysis/nuwro_cc1pi \
    -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$EVENT_SO"

# numu: reuse the existing 2M sample if present
if [[ ! -s nuwro/out_numu_numi.root ]]; then
  echo "[gen] numu ..."; nuwro -i nuwro/params_numu_numi.txt -o nuwro/out_numu_numi.root
fi
echo "[analyze] numu ..."; ./analysis/nuwro_cc1pi nuwro/out_numu_numi.root out/nuwro_numu_pred.root

# numubar
if [[ ! -s nuwro/out_numubar_numi.root ]]; then
  echo "[gen] numubar ..."; nuwro -i nuwro/params_numubar_numi.txt -o nuwro/out_numubar_numi.root
fi
echo "[analyze] numubar ..."; ./analysis/nuwro_cc1pi nuwro/out_numubar_numi.root out/nuwro_numubar_pred.root

echo "[combine] flux-weighted, per-nucleus ..."
root -l -b -q "analysis/combine_nuwro.C(\"out/nuwro_numu_pred.root\",\"out/nuwro_numubar_pred.root\",\"out/nuwro_cc1pi_final.root\")"
echo "[done] out/nuwro_cc1pi_final.root"
