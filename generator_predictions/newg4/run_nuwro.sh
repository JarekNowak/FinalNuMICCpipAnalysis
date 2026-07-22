HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
cd "$HERE/newg4"; echo "[nuwro-newg4] start"
set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup nuwro v21_09_1 -q e20:prof; [[ -z "$NUWRO_FQ_DIR" ]] && { echo fatal; exit 1; }
export NUWRO="$NUWRO_FQ_DIR/nuwro-nuwro_21.09.1"
mkdir -p "$HOME/generators/.nuwro_run/bin"
ln -sfn "$NUWRO_FQ_DIR/bin/nuwro" "$HOME/generators/.nuwro_run/bin/nuwro"
ln -sfn "$NUWRO/data" "$HOME/generators/.nuwro_run/data"
export PATH="$HOME/generators/.nuwro_run/bin:$PATH"
set -e
g++ "$HERE/analysis/nuwro_cc1pi.cc" -o nuwro_cc1pi -I"$NUWRO/src" $(root-config --cflags --libs) -lEG "$NUWRO_FQ_DIR/bin/event1.so"
[[ -s out_numu.root ]]    || { echo "[gen] numu";    nuwro -i params_numu.txt    -o out_numu.root; }
./nuwro_cc1pi out_numu.root    pred_numu.root
[[ -s out_numubar.root ]] || { echo "[gen] numubar"; nuwro -i params_numubar.txt -o out_numubar.root; }
./nuwro_cc1pi out_numubar.root pred_numubar.root
root -l -b -q "combine_newg4.C(\"pred_numu.root\",\"pred_numubar.root\",\"nuwro_newg4_final.root\")"
echo "[nuwro-newg4] done"
