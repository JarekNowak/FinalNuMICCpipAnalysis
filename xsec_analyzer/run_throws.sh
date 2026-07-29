#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; N=15; RES=../logs/throw_results.txt; : > "$RES"
OBS=(pmu ppi costhmu costhpi thmupi)
for t in $(seq 1 $N); do
  echo "===== THROW $t / $N ====="
  root.exe -l -b -q "macros/throw_cv.C($t)" 2>/dev/null | grep seed
  pids=()
  for o in "${OBS[@]}"; do
    FPM=configs/file_properties_numi_cv.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
      OUT=$PROC/ccpi_Run1_${o}_cv_univmake.root ./run_universe_maker.sh >/dev/null 2>&1 &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p"; done
  for o in "${OBS[@]}"; do
    line=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_cv.txt "$o" 2>/dev/null | grep "flux RECO")
    unf=$(echo "$line" | sed -nE 's/.*UNFOLDED:[^|]*bin-avg= *([0-9.]+)%.*/\1/p')
    amp=$(echo "$line" | sed -nE 's/.*amp\(diag\)=([0-9.]+)x.*/\1/p')
    stat=$(echo "$line" | sed -nE 's/.*stat\(unf\)= *([0-9.]+)%.*/\1/p')
    echo "$t $o ${unf:-NA} ${amp:-NA} ${stat:-NA}" >> "$RES"
  done
  echo "  throw $t: $(grep "^$t " "$RES" | awk '{printf "%s=%s ",$2,$4}')"
done
echo "===== THROWS DONE ====="
# summary: mean +/- std of unfolded flux and amp per observable
python3 - "$RES" <<'PY'
import sys,statistics as st
rows=[l.split() for l in open(sys.argv[1]) if l.strip()]
obs=["pmu","ppi","costhmu","costhpi","thmupi"]
print(f"{'obs':8} {'N':>3} {'unfl_flux mean±sd':>20} {'amp mean±sd':>16} {'stat mean':>10}")
for o in obs:
  u=[float(r[2]) for r in rows if r[1]==o and r[2]!='NA']
  a=[float(r[3]) for r in rows if r[1]==o and r[3]!='NA']
  s=[float(r[4]) for r in rows if r[1]==o and r[4]!='NA']
  if u: print(f"{o:8} {len(u):>3} {st.mean(u):>8.1f} ± {st.pstdev(u):>4.1f}%      {st.mean(a):>5.2f} ± {st.pstdev(a):>4.2f}x    {st.mean(s):>6.1f}%")
PY
echo "===== SUMMARY DONE ====="
