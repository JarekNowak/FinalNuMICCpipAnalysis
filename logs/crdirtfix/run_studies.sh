#!/usr/bin/env bash
# re-run of the bkgfit studies of run_crdirtfix.sh after linking detvar_matched.npz next to the templates
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export XSEC_ANALYZER_DIR="$PWD" OMP_NUM_THREADS=4
export BKGFIT_TEMPLATES=/data/uboone/processed/bkgfit_dirtfix/templates.root BKGFIT_LOG=$PWD/../logs/crdirtfix/bkgfit BKGFIT_FIG=$PWD/../logs/crdirtfix/bkgfit_figs
L=../logs/crdirtfix
for s in sensitivity "robustness 30" "validate 300" data; do
  nice -n 10 python3 -m bkgfit.studies $s > $L/bkgfit_$(echo $s | cut -d' ' -f1).log 2>&1 || echo "study $s FAILED"
done
nice -n 10 python3 -m bkgfit.figures > $L/bkgfit_figures.log 2>&1 || echo "figures FAILED"
echo "STUDIES_DONE $(date +%H:%M)"
