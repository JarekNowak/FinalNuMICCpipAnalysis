#!/usr/bin/env bash
# regen_fdfix.sh -- re-run all 51 release extractions after the fake-data-TRUTH weighting
# fix (SystematicsCalculator: an external throw's truth is now the unweighted true
# distribution, as its reco side already was). Unfolded spectra and covariances are
# unchanged; the closure sidecars (h_fakedata_truth), closure chi2 and truth curves are
# not. Writes the RELEASE names (xsec_<TAG>_<obs>, xsec_ccpi1p_<TAG>_<obs>) so every
# downstream artefact is regenerated from them; the previous sidecars are in
# /data/uboone/processed/closure_backup_20260902/.
# Three configurations run in parallel, each in its own scratch directory, because
# UnfolderNuMI writes fixed plot_stepN_*.pdf names into ./unfold_output.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
cd "$REPO"; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
PROC=/data/uboone/processed; FIG=$REPO/../report/figures; DUMP=$REPO/../logs/systdump; LOG=$REPO/../logs/fdfix
SCR=${FDFIX_SCRATCH:-/tmp/fdfix_work}
mkdir -p "$LOG" "$DUMP"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
declare -A PFX=( [fhc5]=""   [rhcfull]="rhc_"  [comb]="comb_" )
declare -A SFX=( [fhc5]=""   [rhcfull]="_RHCFULL" [comb]="_COMB" )

run_cfg() {
  local cfg=$1 t=${TAG[$1]} p=${PFX[$1]} s=${SFX[$1]}
  local W=$SCR/$cfg; mkdir -p "$W/unfold_output"
  ln -sfn "$REPO/configs" "$W/configs"; ln -sfn "$REPO/bin" "$W/bin"; ln -sfn "$REPO/lib" "$W/lib"
  ln -sfn "$REPO/booster_decision_tree" "$W/booster_decision_tree" 2>/dev/null || true
  cd "$W"
  # ---- inclusive (6) ----
  for o in pmu ppi costhmu costhpi thmupi thetamu; do
    if [ "$o" = ppi ]; then xc=configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt; sc=configs/ccpi_ppi_slice_config_2bin.txt; o2=ppi2bin
    else xc=configs/ccpi_xsec_config_numi_${o}_${cfg}.txt; sc=configs/ccpi_${o}_slice_config_opt.txt; o2=$o; fi
    rm -f unfold_output/plot_step[1-4]_*.pdf
    bin/UnfolderNuMI "$xc" "$sc" "$PROC/xsec_${t}_${o2}.root" > "$LOG/${cfg}_${o}.raw" 2>&1 || echo "  FAIL $cfg $o"
    cp "$LOG/${cfg}_${o}.raw" "$DUMP/${cfg}_${o}.raw"; grep '\[SYSTDUMP\]' "$LOG/${cfg}_${o}.raw" > "$DUMP/${cfg}_${o}.dump"
    [ -f unfold_output/plot_step1_reco_spectrum.pdf ] && { cp unfold_output/plot_step1_reco_spectrum.pdf "$FIG/step1_reco_${o}${s}.pdf"; [ "$o" != thetamu ] && cp unfold_output/plot_step1_reco_spectrum.pdf "$FIG/fw_reco_${p}${o}.pdf"; }
    [ -f unfold_output/plot_step2_bkgd_subtraction.pdf ] && cp unfold_output/plot_step2_bkgd_subtraction.pdf "$FIG/step2_bsub_${o}${s}.pdf"
    if [ "$cfg" = fhc5 ] && [ "$o" != thetamu ] && [ "$o" != ppi ]; then
      [ -f unfold_output/plot_step3_smearing_matrix.pdf ] && cp unfold_output/plot_step3_smearing_matrix.pdf "$FIG/fw_resp_${o}.pdf"
      [ -f unfold_output/plot_step4_efficiency.pdf ]      && cp unfold_output/plot_step4_efficiency.pdf      "$FIG/fw_eff_${o}.pdf"
    fi
    echo "  [$(date +%H:%M)] $cfg $o  truth-chi2: $(grep -m1 '^truth:' "$LOG/${cfg}_${o}.raw")"
  done
  # ---- proton-tagged (11) ----
  for o in pmu ppi2bin costhmu costhpi thmupi Wpipr Whad dpt2bin dphit2bin dalphat2bin pn2bin; do
    xc=configs/ccpi1p_xsec_config_numi_${o}_${cfg}.txt
    case "$o" in
      dpt2bin|dphit2bin|dalphat2bin|pn2bin) sc=configs/ccpi1p_${o%2bin}_slice_config_2bin.txt ;;
      ppi2bin) sc=configs/ccpi1p_ppi_slice_config_2bin.txt ;;
      *) sc=configs/ccpi1p_${o}_slice_config.txt ;;
    esac
    rm -f unfold_output/plot_step[1-4]_*.pdf
    bin/UnfolderNuMI "$xc" "$sc" "$PROC/xsec_ccpi1p_${t}_${o}.root" > "$LOG/1p_${cfg}_${o}.raw" 2>&1 || echo "  FAIL 1p $cfg $o"
    cp "$LOG/1p_${cfg}_${o}.raw" "$DUMP/1p_${cfg}_${o}.raw"; grep '\[SYSTDUMP\]' "$LOG/1p_${cfg}_${o}.raw" > "$DUMP/1p_${cfg}_${o}.dump"
    if [ "$cfg" = fhc5 ]; then
      k=${o%2bin}
      case "$o" in Wpipr|Whad|dpt2bin|dalphat2bin|dphit2bin|pn2bin)
        [ -f unfold_output/plot_step1_reco_spectrum.pdf ]    && cp unfold_output/plot_step1_reco_spectrum.pdf    "$FIG/wstep1_reco_${k}.pdf"
        [ -f unfold_output/plot_step2_bkgd_subtraction.pdf ] && cp unfold_output/plot_step2_bkgd_subtraction.pdf "$FIG/wstep2_bsub_${k}.pdf"
        [ -f unfold_output/plot_step3_smearing_matrix.pdf ]  && cp unfold_output/plot_step3_smearing_matrix.pdf  "$FIG/wsmear_${k}.pdf" ;;
      esac
    fi
    echo "  [$(date +%H:%M)] 1p $cfg $o  truth-chi2: $(grep -m1 '^truth:' "$LOG/1p_${cfg}_${o}.raw")"
  done
}

echo "==== FDFIX REGEN START $(date) ===="
for cfg in fhc5 rhcfull comb; do run_cfg "$cfg" > "$LOG/stream_${cfg}.log" 2>&1 & done
wait
echo "==== extractions done $(date) ===="
cd "$REPO"
echo "---- dsigma montages ----"
for cfg in FHC5 RHCFULL COMB; do
  gt=newg4; [ "$cfg" = RHCFULL ] && gt=rhc; [ "$cfg" = COMB ] && gt=comb
  root.exe -l -b -q "macros/dsigma_current.C(\"$cfg\",\"$gt\")" > "$LOG/dsigma_${cfg}.log" 2>&1
  [ -f "unfold_output/dsigma_${cfg}.pdf" ] && cp "unfold_output/dsigma_${cfg}.pdf" "$FIG/dsigma_${cfg}.pdf" && echo "  dsigma $cfg ok" || echo "  dsigma $cfg FAIL"
done
for set in wtki wtki_noW incl; do
  root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"FHC5\",\"$set\")" > "$LOG/dsigma_ccpi1p_${set}.log" 2>&1
done
cp unfold_output/dsigma_ccpi1p_FHC5.pdf "$FIG/" 2>/dev/null; cp unfold_output/dsigma_ccpi1p_noW_FHC5.pdf "$FIG/" 2>/dev/null; cp unfold_output/dsigma_ccpi1p_incl_FHC5.pdf "$FIG/" 2>/dev/null
echo "---- ppi2bin + systbreak figures ----"
for cfg in FHC5 RHCFULL COMB; do root.exe -l -b -q "macros/ppi2bin_figs.C(\"$cfg\")" > "$LOG/ppi2bin_${cfg}.log" 2>&1 || echo "  ppi2bin $cfg FAIL"; done
for cfg in fhc5 rhcfull comb; do root.exe -l -b -q "macros/systbreak_fig.C(\"$cfg\")" > "$LOG/systbreak_${cfg}.log" 2>&1 || echo "  systbreak $cfg FAIL"; done
echo "---- data release ----"
root.exe -l -b -q 'macros/export_curves.C("../report/data_release")'   > "$LOG/export_curves.log"   2>&1 || echo "  export_curves FAIL"
root.exe -l -b -q 'macros/export_matrices.C("../report/data_release")' > "$LOG/export_matrices.log" 2>&1 || echo "  export_matrices FAIL"
echo "==== FDFIX REGEN DONE $(date) ===="
