#!/usr/bin/env bash
# relabel_regen_1p.sh -- proton-tagged (CC1mu1pi1p) counterpart of relabel_regen.sh.
# Regenerates the wstep1_reco_*, wstep2_bsub_* and wsmear_* document figures after the
# "NuMI Data" -> "NuMI fake data" relabelling. FHC only, which is what the documents show.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; FIG=../report/figures; UO=unfold_output; LOG=../logs/relabel
mkdir -p "$LOG"
n=0; fail=0
# The TKI observables are published in the COARSENED two-bin scheme; the fine-binned
# configs still exist but their univmakes were deliberately not rebuilt, so unfolding
# them fails and the step figures silently go missing. Draw from the published config
# while keeping the historical output filename.
for k in Wpipr Whad dpt dalphat dphit pn; do
  case "$k" in
    dpt|dalphat|dphit|pn) src="${k}2bin"; sc=configs/ccpi1p_${k}_slice_config_2bin.txt ;;
    *)                    src="$k";       sc=configs/ccpi1p_${k}_slice_config.txt ;;
  esac
  xc=configs/ccpi1p_xsec_config_numi_${src}_fhc5.txt
  [ -f "$xc" ] && [ -f "$sc" ] || { echo "  MISSING config for $k"; fail=$((fail+1)); continue; }
  rm -f "$UO"/plot_step[1-4]_*.pdf
  bin/UnfolderNuMI "$xc" "$sc" "$PROC/relabel_1p_${k}.root" > "$LOG/1p_${k}.log" 2>&1
  [ -f "$UO/plot_step1_reco_spectrum.pdf" ]    && { cp "$UO/plot_step1_reco_spectrum.pdf"    "$FIG/wstep1_reco_${k}.pdf"; n=$((n+1)); } || { echo "  !! no step1 $k"; fail=$((fail+1)); }
  [ -f "$UO/plot_step2_bkgd_subtraction.pdf" ] && { cp "$UO/plot_step2_bkgd_subtraction.pdf" "$FIG/wstep2_bsub_${k}.pdf"; n=$((n+1)); } || { echo "  !! no step2 $k"; fail=$((fail+1)); }
  [ -f "$UO/plot_step3_smearing_matrix.pdf" ]  && { cp "$UO/plot_step3_smearing_matrix.pdf"  "$FIG/wsmear_${k}.pdf";     n=$((n+1)); } || { echo "  !! no step3 $k"; fail=$((fail+1)); }
  echo "  [$(date +%H:%M)] 1p $k done"
  rm -f "$PROC/relabel_1p_${k}.root"
done
echo "==== proton-tagged figures written: $n, failures: $fail ===="
