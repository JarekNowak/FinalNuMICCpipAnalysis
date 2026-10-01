#!/usr/bin/env bash
# run_crdirtfix.sh -- regenerate the control-region outputs with the per-run dirt normalisation
# (dirt_perrun.h / bkgfit/dirt.py, 2026-10-01). Logs in logs/crdirtfix/; bkgfit templates, JSON and
# figures go to separate locations so the released ones stay in place until review. The beam-on
# figures (sideband_data_*.pdf) and frozen-binning figures (sideband_plots_*.pdf) are overwritten in
# report/figures/ (tracked in git).
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$PWD"
L=../logs/crdirtfix
R="nice -n 10 root -l -b -q"
streamA(){  # control-region yield tables (released: logs/note_2026-09-24, logs/sidebands_perrun)
  for m in fhc rhc comb; do
    $R "macros/sideband_compare.C(\"$m\",\"fake\",\"/data/uboone/processed/sb/\",\"CC1mu1piXp\",false)" > $L/incl_$m.log 2>&1
    $R "macros/sideband_compare.C(\"$m\",\"fakestack\",\"/data/uboone/processed/sb/\",\"CC1mu1piXp\",false)" > $L/fakestack_$m.log 2>&1
    $R "macros/sideband_compare.C(\"$m\",\"fake\",\"/data/uboone/processed/w/\",\"CC1mu1pi1p\",true)" > $L/1p_tag_$m.log 2>&1
  done; echo "A done $(date +%H:%M)"; }
streamB(){  # figures and per-run beam-on table
  for m in fhc rhc; do $R "macros/sb_plots.C(\"$m\")" > $L/sb_plots_$m.log 2>&1; done
  for m in fhc rhc; do for c in false true; do $R "macros/sb_plots_data.C(\"$m\",$c)" >> $L/sb_plots_data_$m.log 2>&1; done; done
  $R macros/sb_perrun_data.C > $L/sb_perrun_data.log 2>&1
  for m in fhc rhc; do $R "macros/pi0final/sb_pi0_final_yields.C(\"$m\")" > $L/pi0final_$m.log 2>&1; done
  echo "B done $(date +%H:%M)"; }
streamC(){  # background-fit templates, studies, figures; far sidebands; Run-1 trigger record
  export BKGFIT_TEMPLATES=/data/uboone/processed/bkgfit_dirtfix/templates.root BKGFIT_LOG=$PWD/../logs/crdirtfix/bkgfit BKGFIT_FIG=$PWD/../logs/crdirtfix/bkgfit_figs
  mkdir -p "$BKGFIT_LOG" "$BKGFIT_FIG"
  $R 'macros/bkgfit_templates.C+("/data/uboone/processed/bkgfit_dirtfix")' > $L/templates.log 2>&1 || { echo "templates FAILED"; return 1; }
  export OMP_NUM_THREADS=4
  for s in composition sensitivity "robustness 30" "validate 300" data; do
    nice -n 10 python3 -m bkgfit.studies $s > $L/bkgfit_$(echo $s | cut -d' ' -f1).log 2>&1 || echo "study $s FAILED"
  done
  nice -n 10 python3 -m bkgfit.figures > $L/bkgfit_figures.log 2>&1 || echo "figures FAILED"
  nice -n 10 python3 -m bkgfit.farsb > $L/farsb.log 2>&1 || echo "farsb FAILED"
  nice -n 10 python3 -m bkgfit.run1_trigger > $L/run1_trigger.log 2>&1 || echo "run1_trigger FAILED"
  echo "C done $(date +%H:%M)"; }
echo "start $(date +%H:%M)"
streamA & streamB & streamC & wait
echo "CRDIRTFIX_DONE $(date +%H:%M)"
