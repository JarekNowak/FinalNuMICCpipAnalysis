#!/usr/bin/env bash
# relabel_regen.sh -- regenerate every document figure that UnfolderNuMI draws, after the
# "NuMI Data" -> "NuMI fake data" relabelling. One unfold per (config, observable); each
# unfold writes fixed plot_stepN_*.pdf names which are copied to the document figure names.
#
# p_pi uses the ADOPTED TWO-BIN configs (ppi2bin + slice_config_2bin), not the superseded
# five-bin ones, so these figures stay consistent with the rest of the analysis.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; FIG=../report/figures; UO=unfold_output; LOG=../logs/relabel
mkdir -p "$LOG"

declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
declare -A PFX=( [fhc5]=""   [rhcfull]="rhc_"  [comb]="comb_" )
declare -A SFX=( [fhc5]=""   [rhcfull]="_RHCFULL" [comb]="_COMB" )

n=0; fail=0
for cfg in fhc5 rhcfull comb; do
  t=${TAG[$cfg]}; p=${PFX[$cfg]}; s=${SFX[$cfg]}
  for o in pmu ppi costhmu costhpi thmupi thetamu; do

    # p_pi is extracted at the two-bin scheme
    if [ "$o" = ppi ]; then
      xc=configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt
      sc=configs/ccpi_ppi_slice_config_2bin.txt
    else
      xc=configs/ccpi_xsec_config_numi_${o}_${cfg}.txt
      sc=configs/ccpi_${o}_slice_config_opt.txt
    fi
    [ -f "$xc" ] && [ -f "$sc" ] || { echo "  MISSING config for $cfg $o"; fail=$((fail+1)); continue; }

    rm -f "$UO"/plot_step[1-4]_*.pdf
    bin/UnfolderNuMI "$xc" "$sc" "$PROC/relabel_${t}_${o}.root" > "$LOG/${cfg}_${o}.log" 2>&1

    # step1 -> the reco-spectrum figures
    if [ -f "$UO/plot_step1_reco_spectrum.pdf" ]; then
      cp "$UO/plot_step1_reco_spectrum.pdf" "$FIG/step1_reco_${o}${s}.pdf"; n=$((n+1))
      # fw_reco_* exists for the five cross-section observables, not for thetamu
      [ "$o" != thetamu ] && { cp "$UO/plot_step1_reco_spectrum.pdf" "$FIG/fw_reco_${p}${o}.pdf"; n=$((n+1)); }
    else echo "  !! no step1 for $cfg $o"; fail=$((fail+1)); fi

    # step2 -> background-subtracted
    if [ -f "$UO/plot_step2_bkgd_subtraction.pdf" ]; then
      cp "$UO/plot_step2_bkgd_subtraction.pdf" "$FIG/step2_bsub_${o}${s}.pdf"; n=$((n+1))
    else echo "  !! no step2 for $cfg $o"; fail=$((fail+1)); fi

    # response and efficiency are FHC-only in the documents, and not kept for thetamu
    if [ "$cfg" = fhc5 ] && [ "$o" != thetamu ]; then
      [ -f "$UO/plot_step3_smearing_matrix.pdf" ] && { cp "$UO/plot_step3_smearing_matrix.pdf" "$FIG/fw_resp_${o}.pdf"; n=$((n+1)); }
      [ -f "$UO/plot_step4_efficiency.pdf" ]      && { cp "$UO/plot_step4_efficiency.pdf"      "$FIG/fw_eff_${o}.pdf";  n=$((n+1)); }
    fi
    echo "  [$(date +%H:%M)] $cfg $o done"
    rm -f "$PROC/relabel_${t}_${o}.root"
  done
done
echo "==== figures written: $n, failures: $fail ===="
