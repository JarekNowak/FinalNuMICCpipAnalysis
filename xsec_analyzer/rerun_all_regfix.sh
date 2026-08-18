#!/usr/bin/env bash
# rerun_all_regfix.sh -- regenerate every extraction with the width-aware Wiener-SVD
# regularisation (commit 8383556). Unfolds only: the univmake files are stage 2 and are
# not affected by the fix, so nothing needs rebuilding there.
#
# Writes to the CANONICAL output names so the closure files the figures and note read are
# the ones refreshed. The previous set is backed up in /data/uboone/processed/pre_regfix_backup.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; LOG=../logs/regfix; mkdir -p "$LOG"
SUM=$LOG/summary.txt; : > "$SUM"

declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
ok=0; fail=0

run() { # xsec_config slice_config output_stem label
  local xc=$1 sc=$2 out=$3 lbl=$4
  if [ ! -f "$xc" ] || [ ! -f "$sc" ]; then echo "  SKIP $lbl (missing config)"; return; fi
  bin/UnfolderNuMI "$xc" "$sc" "$PROC/${out}.root" > "$LOG/${lbl}.log" 2>&1
  local w=$(grep -c "REGWIDTH] supplied" "$LOG/${lbl}.log")
  local s=$(grep -oE "SYSTDUMP\] sigma_int [0-9.eE+-]+" "$LOG/${lbl}.log" | awk '{print $3}')
  local c=$(grep -oE "truth: .* = [0-9.]+/[0-9]+ bins, p-value = [0-9.]+" "$LOG/${lbl}.log" | tail -1 \
            | sed -E 's/.*= ([0-9.]+\/[0-9]+) bins, p-value = ([0-9.]+)/\1 p=\2/')
  local i=$(grep -oE "IDENTITY\] .*relative [0-9.eE+-]+" "$LOG/${lbl}.log" | awk '{print $NF}')
  if [ -n "$s" ]; then ok=$((ok+1)); else fail=$((fail+1)); fi
  printf "%-26s widths=%s sigma_int=%-10s closure=%-22s identity=%s\n" \
         "$lbl" "${w:-0}" "${s:-FAIL}" "${c:-n/a}" "${i:-n/a}" | tee -a "$SUM"
  rm -f "$PROC/${out}.root"
}

echo "==== INCLUSIVE ====" | tee -a "$SUM"
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}
  for o in pmu costhmu costhpi thmupi thetamu; do
    run configs/ccpi_xsec_config_numi_${o}_${cfg}.txt configs/ccpi_${o}_slice_config_opt.txt \
        xsec_${t}_${o} "${cfg}_${o}"
  done
  run configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt configs/ccpi_ppi_slice_config_2bin.txt \
      xsec_${t}_ppi2bin "${cfg}_ppi2bin"
done

echo "==== PROTON-TAGGED ====" | tee -a "$SUM"
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}
  for k in Wpipr Whad dpt dalphat dphit pn; do
    run configs/ccpi1p_xsec_config_numi_${k}_${cfg}.txt configs/ccpi1p_${k}_slice_config.txt \
        xsec_ccpi1p_${t}_${k} "1p_${cfg}_${k}"
  done
  run configs/ccpi1p_xsec_config_numi_ppi2bin_${cfg}.txt configs/ccpi1p_ppi_slice_config_2bin.txt \
      xsec_ccpi1p_${t}_ppi2bin "1p_${cfg}_ppi2bin"
done

echo "==== DONE: $ok succeeded, $fail failed ====" | tee -a "$SUM"
