#!/bin/bash
# harvest_proposed_local.sh -- re-unfold and harvest covariance matrices into the data release for the
# results proposed for approval on 2026-09-27 that were outside the release index: the proton polar
# angle theta_p and the pion-proton opening angle theta_pi_p (FHC, RHC, combined) and the two
# two-dimensional results in the combined configuration, theta_p x delta p_T (proton-tagged) and
# cos theta_pi x cos theta_mu (inclusive). Same per-extraction body as unfold_1p_local.sh; the universe
# files are unchanged, so the unfolded results and sidecars are reproduced, not changed.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RELEASE="$(cd "$REPO/../report" && pwd)/data_release"
PROC=/data/uboone/processed; SCRATCH=$PROC/unfold_scratch
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
ok=0; fail=0
# tag | xsec config | slice config | output root
while IFS='|' read -r TAGN XSEC SLICE OUT; do
  [ -z "$TAGN" ] && continue
  W="$SCRATCH/prop_${TAGN}.$$"; rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XSEC" "$SLICE" "$OUT" ) > "$W/run.log" 2>&1
  rc=$?; n=$(ls "$W"/unfold_output/mat_table_cov_*.txt 2>/dev/null | wc -l)
  if [ $rc -ne 0 ] || [ "$n" -eq 0 ]; then echo "  FAILED $TAGN (rc=$rc, $n)"; cp "$W/run.log" "$REPO/../logs/fdfix/harvest_proposed_${TAGN}.log"; fail=$((fail+1)); rm -rf "$W"; continue; fi
  OUTD="$RELEASE/cov/$TAGN"; mkdir -p "$OUTD"
  for f in "$W"/unfold_output/mat_table_cov_*.txt; do cp -p "$f" "$OUTD/$(basename "$f" | sed 's/^mat_table_//')"; done
  for extra in add_smear unfolding err_prop; do [ -f "$W/unfold_output/mat_table_${extra}.txt" ] && cp -p "$W/unfold_output/mat_table_${extra}.txt" "$OUTD/${extra}.txt"; done
  [ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$OUTD/unfolded_signal.txt"
  cp "$W/run.log" "$REPO/../logs/fdfix/harvest_proposed_${TAGN}.log"
  echo "  OK $TAGN ($n cov matrices)"; ok=$((ok+1)); rm -rf "$W"
done <<LIST
1p_FHC5_thetap|configs/ccpi1p_xsec_config_numi_thetap_fhc5.txt|configs/ccpi1p_thetap_slice_config.txt|$PROC/xsec_ccpi1p_FHC5_thetap.root
1p_RHCFULL_thetap|configs/ccpi1p_xsec_config_numi_thetap_rhcfull.txt|configs/ccpi1p_thetap_slice_config.txt|$PROC/xsec_ccpi1p_RHCFULL_thetap.root
1p_COMB_thetap|configs/ccpi1p_xsec_config_numi_thetap_comb.txt|configs/ccpi1p_thetap_slice_config.txt|$PROC/xsec_ccpi1p_COMB_thetap.root
1p_FHC5_thpipr|configs/ccpi1p_xsec_config_numi_thpipr_fhc5.txt|configs/ccpi1p_thpipr_slice_config.txt|$PROC/xsec_ccpi1p_FHC5_thpipr.root
1p_RHCFULL_thpipr|configs/ccpi1p_xsec_config_numi_thpipr_rhcfull.txt|configs/ccpi1p_thpipr_slice_config.txt|$PROC/xsec_ccpi1p_RHCFULL_thpipr.root
1p_COMB_thpipr|configs/ccpi1p_xsec_config_numi_thpipr_comb.txt|configs/ccpi1p_thpipr_slice_config.txt|$PROC/xsec_ccpi1p_COMB_thpipr.root
1p_COMB_thetap_dpt2d|configs/d2/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb.txt|configs/d2/ccpi1p_thetap_dpt_slice_config_2d.txt|$PROC/rebuild_2d/xsec_ccpi1p_COMB_thetap_dpt.root
incl_COMB_costhpi_costhmu2d|configs/d2/ccpi_xsec_config_numi_costhpi_costhmu_2d_comb.txt|configs/d2/ccpi_costhpi_costhmu_slice_config_2d.txt|$PROC/rebuild_2d/xsec_ccpi_COMB_costhpi_costhmu.root
LIST
echo "  ---- $ok ok, $fail failed ----"
