#!/usr/bin/env bash
# regen_figs.sh — regenerate ALL per-run note figures + the systdump tables from the
# completed per-run univmakes. One UnfolderNuMI pass per config x observable (reads the
# per-run univmake named in each xsec config -> identical to the batch result), capturing:
#   - stdout [SYSTDUMP] lines -> logs/systdump/<cfg>_<obs>.dump   (for systbreak_fig.C)
#   - closure_hists_xsec_<TAG>_<obs>.root                         (for dsigma_current.C)
#   - plot_step1 -> fw_reco[_<pfx>]<obs>  (all configs)
#   - plot_step3 -> fw_resp_<obs>, plot_step4 -> fw_eff_<obs>     (FHC only, per note)
# then draws the dsigma montages (3) + systbreak figures (3) + cutflow yields (3).
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
PROC=/data/uboone/processed; FIG=../report/figures; DUMP=../logs/systdump; UO=unfold_output
mkdir -p "$DUMP" "$UO"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
declare -A PFX=( [fhc5]="" [rhcfull]="rhc_" [comb]="comb_" )
declare -A GT=(  [fhc5]=newg4 [rhcfull]=rhc [comb]=comb )
obs=(pmu ppi costhmu costhpi thmupi thetamu)

echo "==== REGEN FIGS START $(date) ===="
nr=0
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}; p=${PFX[$cfg]}
  for o in "${obs[@]}"; do
    xc=configs/ccpi_xsec_config_numi_${o}_${cfg}.txt
    sc=configs/ccpi_${o}_slice_config_opt.txt
    [ -f "$xc" ] || { echo "  MISSING xsec config $xc"; continue; }
    bin/UnfolderNuMI "$xc" "$sc" $PROC/xsec_${t}_${o}.root > "$DUMP/${cfg}_${o}.raw" 2>&1
    grep '\[SYSTDUMP\]' "$DUMP/${cfg}_${o}.raw" > "$DUMP/${cfg}_${o}.dump"
    [ -f "$UO/plot_step1_reco_spectrum.pdf" ] && { cp "$UO/plot_step1_reco_spectrum.pdf" "$FIG/fw_reco_${p}${o}.pdf"; nr=$((nr+1)); }
    if [ "$cfg" = "fhc5" ]; then
      [ -f "$UO/plot_step3_smearing_matrix.pdf" ] && cp "$UO/plot_step3_smearing_matrix.pdf" "$FIG/fw_resp_${o}.pdf"
      [ -f "$UO/plot_step4_efficiency.pdf" ]      && cp "$UO/plot_step4_efficiency.pdf"      "$FIG/fw_eff_${o}.pdf"
    fi
    echo "  [$(date +%H:%M)] $cfg $o  sigma_int=$(awk '$2=="sigma_int"{print $3}' "$DUMP/${cfg}_${o}.dump")"
  done
done
echo "reco spectra copied: $nr/15"

echo "---- dsigma montages ----"
for cfg in FHC5 RHCFULL COMB; do
  gt=newg4; [ "$cfg" = RHCFULL ] && gt=rhc; [ "$cfg" = COMB ] && gt=comb
  root.exe -l -b -q "macros/dsigma_current.C(\"$cfg\",\"$gt\")" > "$DUMP/dsigma_${cfg}.log" 2>&1
  [ -f "$UO/dsigma_${cfg}.pdf" ] && cp "$UO/dsigma_${cfg}.pdf" "$FIG/dsigma_${cfg}.pdf" \
    && echo "  dsigma $cfg ok" || echo "  dsigma $cfg FAIL"
done
echo "---- systbreak figures ----"
for cfg in fhc5 rhcfull comb; do
  root.exe -l -b -q "macros/systbreak_fig.C(\"$cfg\")" > "$DUMP/systbreak_${cfg}.log" 2>&1 \
    && echo "  systbreak $cfg ok" || echo "  systbreak $cfg FAIL"
done
echo "---- cutflow yields ----"
for m in fhc rhc comb; do
  root.exe -l -b -q "macros/cutflow_yields.C(\"$m\")" > "$DUMP/cutflow_${m}.log" 2>&1
  [ -f "$UO/cutflow_yields_${m}.pdf" ] && cp "$UO/cutflow_yields_${m}.pdf" "$FIG/cutflow_yields_${m}.pdf" \
    && echo "  cutflow $m ok" || echo "  cutflow $m FAIL"
done
echo "==== REGEN FIGS DONE $(date) ===="
