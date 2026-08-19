#!/usr/bin/env bash
# Regenerate the per-config/observable reco-spectrum figures (data vs MC stack, now
# Okabe-Ito) the note references as fw_reco[_<cfg>_]<obs>. UnfolderNuMI writes a fixed
# plot_step1_reco_spectrum.pdf, so re-run each unfold and copy it to the target name.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed; FIG=../report/figures
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
declare -A PFX=( [fhc5]="" [rhcfull]="rhc_" [comb]="comb_" )
n=0
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}; p=${PFX[$cfg]}
  for o in pmu ppi costhmu costhpi thmupi thetamu; do
    bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_${cfg}.txt configs/ccpi_${o}_slice_config_opt.txt \
      $PROC/xsec_${t}_${o}.root > /dev/null 2>&1
    if [ -f unfold_output/plot_step1_reco_spectrum.pdf ]; then
      cp unfold_output/plot_step1_reco_spectrum.pdf "$FIG/fw_reco_${p}${o}.pdf"; n=$((n+1))
    else echo "  MISSING plot_step1 for $cfg $o"; fi
  done
done
echo "reco spectra written: $n/18"
ls "$FIG"/fw_reco_*.pdf 2>/dev/null | wc -l
echo "########## RECO SPECTRA DONE ##########"
