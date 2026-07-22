HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/gibuu
export GIBUU_ROOT="$HOME/generators/gibuu"; export BUUINPUT="$GIBUU_ROOT/buuinput"
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup gcc v9_3_0; export PATH="$GIBUU_ROOT/release/objects:$PATH"
command -v GiBUU.x || { echo "no GiBUU.x"; exit 1; }
for tag in numu numubar; do
  mkdir -p "$HERE/run_$tag"; cd "$HERE/run_$tag"
  echo "[gibuu] $tag ..."; GiBUU.x < "$HERE/job_$tag.job" > gibuu_$tag.log 2>&1
  root -l -b -q "$HERE/gibuu_cc1pi.C(\"FinalEvents.dat\",20,\"$HERE/pred_$tag.root\")" 2>&1 | grep -iE "CC1pi|wrote"
done
cd "$HERE"; root -l -b -q "../combine_newg4.C(\"pred_numu.root\",\"pred_numubar.root\",\"gibuu_newg4_final.root\")" 2>&1 | grep integral
echo "[gibuu] done"
