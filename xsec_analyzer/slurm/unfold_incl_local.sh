#!/bin/bash
# Re-unfold the 18 inclusive extractions to produce their covariance matrices and refresh
# their A_C. Run locally: each unfold takes well under a minute, so 18 of them is faster
# than queueing, and this keeps the cluster free.
#
# Each unfold still gets its own working directory: UnfolderNuMI writes its matrix dumps
# to a FIXED relative path, so sharing one would silently overwrite results.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RELEASE="$(cd "$REPO/../report" && pwd)/data_release"
PROC=/data/uboone/processed
SCRATCH=$PROC/unfold_scratch
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
ok=0; fail=0
for cfg in fhc5 rhcfull comb; do
  for obs in pmu ppi costhmu costhpi thmupi thetamu; do
    # p_pi is measured in the two-bin scheme; the five-bin configs are withdrawn.
    if [ "$obs" = ppi ]; then
      XSEC="configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt"
      SLICE="configs/ccpi_ppi_slice_config_2bin.txt"; OUTOBS=ppi2bin
    else
      XSEC="configs/ccpi_xsec_config_numi_${obs}_${cfg}.txt"
      SLICE="configs/ccpi_${obs}_slice_config_opt.txt"; OUTOBS=$obs
    fi
    for f in "$REPO/$XSEC" "$REPO/$SLICE"; do
      [ -f "$f" ] || { echo "  SKIP $cfg $obs (missing $(basename $f))"; continue 2; }
    done
    W="$SCRATCH/incl_${TAG[$cfg]}_${OUTOBS}.$$"
    rm -rf "$W"; mkdir -p "$W/unfold_output"
    ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
    # The closure sidecar name is derived from this output filename
    # (closure_hists_<name>), so it must be the real one -- passing /dev/null silently
    # leaves the sidecar, and hence the released A_C, un-refreshed.
    ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XSEC" "$SLICE" \
        "$PROC/xsec_${TAG[$cfg]}_${OUTOBS}.root" ) > "$W/run.log" 2>&1
    rc=$?
    n=$(ls "$W"/unfold_output/mat_table_cov_*.txt 2>/dev/null | wc -l)
    if [ $rc -ne 0 ] || [ "$n" -eq 0 ]; then
      echo "  FAILED $cfg $OUTOBS (rc=$rc, $n matrices)"; fail=$((fail+1)); rm -rf "$W"; continue
    fi
    OUTD="$RELEASE/cov/incl_${TAG[$cfg]}_${OUTOBS}"; mkdir -p "$OUTD"
    for f in "$W"/unfold_output/mat_table_cov_*.txt; do
      cp -p "$f" "$OUTD/$(basename "$f" | sed 's/^mat_table_//')"
    done
    for extra in add_smear unfolding err_prop; do
      [ -f "$W/unfold_output/mat_table_${extra}.txt" ] && cp -p "$W/unfold_output/mat_table_${extra}.txt" "$OUTD/${extra}.txt"
    done
    [ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$OUTD/unfolded_signal.txt"
    sig=$(grep 'SYSTDUMP] sigma_int' "$W/run.log" | tail -1 | awk '{print $3}')
    echo "  OK $cfg $OUTOBS  sigma_int=$sig  ($n cov matrices)"
    ok=$((ok+1)); rm -rf "$W"
  done
done
echo "  ---- $ok ok, $fail failed ----"
