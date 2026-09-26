#!/usr/bin/env bash
# dagostini_050.sh -- D'Agostini (iterative Bayesian) cross-check of the inclusive release after the binnings
# of the 0.50 migration criterion were adopted (2026-09-26): the five inclusive observables in FHC, RHC and
# combined, 1, 2 and 4 iterations, on the live universe files, fake data and flux of the Wiener-SVD
# extractions. Each configuration is derived from the live Wiener-SVD config by changing only the Unfold
# line, so it cannot drift from the release (the old *_dag configs carried stale prediction binnings).
# Run after promote_050.sh.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
PROC=/data/uboone/processed; LOG=../logs/dag050; mkdir -p "$LOG"
SCR=$PROC/unfold_scratch/dag050; mkdir -p "$SCR"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
run_cfg() {
  local cfg=$1 t=${TAG[$1]} W=$SCR/$1; mkdir -p "$W/unfold_output"
  ln -sfn "$REPO/configs" "$W/configs"; ln -sfn "$REPO/bin" "$W/bin"; ln -sfn "$REPO/lib" "$W/lib"; cd "$W"
  for it in 1 2 4; do
    for o in pmu ppi3bin costhmu costhpi thmupi; do
      if [ $o = ppi3bin ]; then sc=configs/ccpi_ppi_slice_config_3bin.txt; else sc=configs/ccpi_${o}_slice_config_opt.txt; fi
      c="$W/dag_${o}_it${it}.txt"
      sed -E "s|^Unfold .*|Unfold DAgostini iter ${it}|" "$REPO/configs/ccpi_xsec_config_numi_${o}_${cfg}.txt" > "$c"
      "$REPO/bin/UnfolderNuMI" "$c" "$sc" "$PROC/xsec_${t}dag${it}_${o}.root" > "$REPO/$LOG/${cfg}_${o}_it${it}.log" 2>&1 || echo "  FAIL $cfg $o it$it"
    done
  done
  echo "  $cfg done $(date +%H:%M)"
}
for cfg in fhc5 rhcfull comb; do run_cfg $cfg & done; wait
cd "$REPO"
python3 - <<'PYEOF' | tee "$LOG/summary.txt"
import uproot
P='/data/uboone/processed/'
def ig(f, h):
    try:
        x = uproot.open(f)[h]; return float((x.values() * (x.axis().edges()[1:] - x.axis().edges()[:-1])).sum())
    except Exception: return float('nan')
print("Wiener-SVD against D'Agostini, integrated unfolded sigma (1e-38 cm^2/Ar)")
for c in ('FHC5', 'RHCFULL', 'COMB'):
    print(f"[{c}] {'obs':8s} {'WSVD':>7s} | {'DAg i1':>7s} {'DAg i2':>7s} {'DAg i4':>7s} | i4/WSVD")
    for o in ('pmu', 'ppi3bin', 'costhmu', 'costhpi', 'thmupi'):
        w = ig(f'{P}closure_hists_xsec_{c}_{o}.root', 'h_unfolded_nuwro')
        d = [ig(f'{P}closure_hists_xsec_{c}dag{i}_{o}.root', 'h_unfolded_nuwro') for i in (1, 2, 4)]
        print(f"      {o:8s} {w:7.3f} | {d[0]:7.3f} {d[1]:7.3f} {d[2]:7.3f} | {d[2]/w:6.3f}")
PYEOF
echo "==== DAGOSTINI 050 DONE $(date) ===="
