#!/usr/bin/env bash
# regen_release.sh -- the downstream chain of promote_run1fix.sh without the promotion: re-unfold all
# 51 extractions from the live universe files, harvest covariances, export the release, regenerate
# the figures and the counting; then the cut-flow yields from the cf/ and cf_1p/ reprocesses (full CV
# weight) whose logs feed report/tools/cutflow_table.py and cutflow_1p_table.py.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; S=${FDFIX_SCRATCH:-$PROC/unfold_scratch/fdfix_work}
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
mkdir -p "$S" ../logs/fdfix
mkdir -p "$S" ../logs/fdfix
FDFIX_SCRATCH=$S ./regen_fdfix.sh > ../logs/fdfix/regen_run1fix.log 2>&1
bash slurm/unfold_incl_local.sh > ../logs/fdfix/harvest_incl_run1fix.log 2>&1 &
bash slurm/unfold_1p_local.sh   > ../logs/fdfix/harvest_1p_run1fix.log 2>&1 &
wait
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
L=../logs/fdfix; FIG=../report/figures
root.exe -l -b -q 'macros/export_curves.C("../report/data_release")' > $L/export_curves_run1fix.log 2>&1
root.exe -l -b -q 'macros/export_matrices.C("../report/data_release")' > $L/export_matrices_run1fix.log 2>&1
for t in FHC5 RHCFULL COMB; do cp -p $PROC/closure_hists_xsec_${t}_ppi2bin.root $PROC/closure_hists_xsec_${t}_ppi.root; root.exe -l -b -q "macros/ppi2bin_figs.C(\"$t\")" > $L/ppi2bin_run1fix_$t.log 2>&1; done
for set in wtki wtki_noW incl; do root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"FHC5\",\"$set\")" > $L/ccpi1p_run1fix_$set.log 2>&1; done
for c in RHCFULL COMB; do for set in wtki wtki_noW incl; do root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"$c\",\"$set\")" > $L/ccpi1p_run1fix_${set}_$c.log 2>&1; done; done
root.exe -l -b -q 'macros/dsigma_build.C("FHC5")' > $L/dsigma_build_run1fix.log 2>&1; root.exe -l -b -q 'macros/dsigma_build_1p.C("FHC5")' > $L/dsigma_build_1p_run1fix.log 2>&1
for c in fhc5 rhcfull comb; do root.exe -l -b -q "macros/systbreak_fig.C(\"$c\")" > $L/systbreak_run1fix_$c.log 2>&1; done
for m in fhc rhc comb; do root.exe -l -b -q "macros/cutflow_yields.C(\"$m\")" > $L/cutflow_run1fix_$m.log 2>&1; cp unfold_output/cutflow_yields_$m.pdf $FIG/; done
root -l -b -q 'macros/total_xsec_counting.C' > ../logs/counting_review7.log 2>&1
cp -p unfold_output/dsigma_ccpi1p_*.pdf unfold_output/dsigma_build_FHC5_*.pdf unfold_output/dsigma_build_1p_FHC5_*.pdf $FIG/ 2>/dev/null
for m in fhc rhc comb; do
  root.exe -l -b -q "macros/cutflow_yields.C(\"$m\",\"$PROC/cf/\",true)" > $L/cutflow_incl_$m.log 2>&1; cp unfold_output/cutflow_yields_$m.pdf $FIG/
  root.exe -l -b -q "macros/cutflow_yields_1p.C(\"$m\")" > $L/cutflow_1p_$m.log 2>&1
done
python3 check_staleness.py | sed -n 1,12p
echo "==== REGEN RELEASE DONE $(date) ===="
