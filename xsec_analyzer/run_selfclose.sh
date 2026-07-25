#!/usr/bin/env bash
# run_selfclose.sh — numuMC self-closure: numuMC in the onBNB (fake-data) slot,
# so D = R*x_true by construction (same sample builds the response and the data).
# Isolates the pure regularization/statistics unfolding bias from the detVar-CV
# fake-data response mismatch. Rebuilds the univmakes, then reports x_hat/truth
# per observable for identity / second-deriv / D'Agostini via ac_diagnostic.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$PWD/lib:${LD_LIBRARY_PATH:-}"
PROC=/data/uboone/processed; LOG=../logs; mkdir -p "$LOG"
OBS=(pmu costhmu costhpi ppi thmupi)

echo "########## STAGE A: univmake (parallel, numuMC self-closure) ##########"
pids=()
for o in "${OBS[@]}"; do
  FPM=configs/file_properties_numi_selfclose.txt \
  BIN_CONFIG=configs/ccpi_${o}_bin_config.txt \
  OUT=$PROC/ccpi_Run1_${o}_selfclose_univmake.root \
  ./run_universe_maker.sh > "$LOG/selfclose_${o}_univ.log" 2>&1 &
  pids+=($!); echo "  launched univmake[$o] pid=${pids[-1]}"
done
fail=0; i=0
for o in "${OBS[@]}"; do
  if wait "${pids[$i]}"; then echo "  univmake[$o] OK"; else echo "  univmake[$o] FAILED"; fail=1; fi
  i=$((i+1))
done

echo "########## STAGE B: ac_diagnostic per obs/regularization ##########"
printf "%-9s %-14s %10s %10s\n" "obs" "reg" "xhat/truth" "A_C_diag" | tee $LOG/selfclose_results.txt
for o in "${OBS[@]}"; do
  univ=$PROC/ccpi_Run1_${o}_selfclose_univmake.root
  [[ -s "$univ" ]] || { echo "  [skip] $o"; continue; }
  grep -q 'univ FakeData' configs/ccpi_xsec_config_numi_${o}.txt || echo 'Prediction FakeData "Fakedata" univ FakeData' >> configs/ccpi_xsec_config_numi_${o}.txt
  base=configs/ccpi_xsec_config_numi_${o}_selfclose.txt
  sed -e "s#ccpi_Run1_${o}_univmake.root#ccpi_Run1_${o}_selfclose_univmake.root#" \
      -e "s#^FPFile configs/file_properties_numi.txt#FPFile configs/file_properties_numi_selfclose.txt#" \
      configs/ccpi_xsec_config_numi_${o}.txt > "$base"
  for spec in "WienerSVD 1 identity" "WienerSVD 1 second-deriv" "DAgostini iter 4"; do
    sed "s/^Unfold .*/Unfold $spec/" "$base" > /tmp/sc_cfg_${o}.txt
    out=$(./bin/ac_diagnostic /tmp/sc_cfg_${o}.txt 2>/dev/null || /tmp/ac_diagnostic /tmp/sc_cfg_${o}.txt 2>/dev/null)
    xt=$(echo "$out"|grep "xhat / truth"|grep -oE "= [0-9.]+"|head -1|tr -d '= ')
    diag=$(echo "$out"|grep -oE "diagonal mean = [0-9.]+"|grep -oE "[0-9.]+$")
    printf "%-9s %-14s %10s %10s\n" "$o" "${spec##* }" "$xt" "$diag" | tee -a $LOG/selfclose_results.txt
  done
done
echo "########## SELFCLOSE DONE (fail=$fail) ##########"
