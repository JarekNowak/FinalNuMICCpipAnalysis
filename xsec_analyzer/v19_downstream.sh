#!/usr/bin/env bash
# v19_downstream.sh -- version 1.9 (MCS calibration): everything downstream of the extractions, after braggfix_regen.sh
# has re-unfolded the release from the configs that point at the calibrated universes (/data/uboone/processed/mcscal/univ):
#   1. regression: the forty results that do not use the muon momentum must equal version 1.8 bit for bit
#   2. the MCS-calibration term into the data release (report/tools/mcs_cal_release.py)
#   3. closure summary, results tables, 2D export, extraction index
#   4. counting scale factors (= official prediction-total fraction of p_mu) and the counting log
#   5. generated tables and figures, normalisation manifest, D'Agostini cross-check
#   6. alternative-model evaluation and ensemble tables
# Requires the calibration term per result (scripts/mcs_cal_compare.py on mcscal/nominal_out and mcscal/var).
#   ./v19_downstream.sh > ../logs/v19/downstream.log 2>&1
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer; cd "$REPO"
L=logs/v19; mkdir -p "$L"
CAL="incl_FHC5_pmu incl_RHCFULL_pmu incl_COMB_pmu 1p_COMB_thetap_dpt2d"
for c in FHC5 RHCFULL COMB; do for o in pmu Whad dpt2bin dalphat2bin pn2bin; do CAL="$CAL 1p_${c}_$o"; done; done
echo "==== 1. regression against version 1.8 $(date)"
n=0; bad=0
for d in report/data_release/cov/*/; do
  t=$(basename "$d"); case " $CAL " in *" $t "*) continue;; esac; n=$((n+1))
  for f in unfolded_signal.txt cov_total.txt; do
    cmp -s "$d$f" <(git show "release/v1.8:$d$f") || { bad=$((bad+1)); echo "  CHANGED $t/$f"; }
  done
done
echo "  $n unaffected results, $bad changed files"
echo "==== 2. MCS-calibration term $(date)"
python3 report/tools/mcs_cal_release.py
echo "==== 3. closure, results, exports $(date)"
( cd report && root -l -b -q tools/closure_tables.C ) > $L/closure_tables.log 2>&1; tail -1 $L/closure_tables.log
python3 report/tools/make_results_tsv.py report/results/closure_summary.tsv > $L/make_results.log 2>&1; tail -1 $L/make_results.log
python3 report/tools/export_2d.py > $L/export_2d.log 2>&1; tail -1 $L/export_2d.log
python3 report/tools/index_extractions.py > $L/index_extractions.log 2>&1; tail -1 $L/index_extractions.log
echo "==== 4. counting $(date)"
python3 - <<'PYEOF'
import re, sys
sys.path.insert(0, 'report/tools')
import official_cov as oc
sf = [[round(oc.summary(f, c, 'pmu')['PredTotal'] / 100, 3) for c in ('FHC5', 'RHCFULL', 'COMB')] for f in ('incl', '1p')]
p = 'xsec_analyzer/macros/total_xsec_counting.C'; s = open(p).read()
s = re.sub(r'const double sf\[2\]\[3\] = \{\{[^}]*\},\{[^}]*\}\};',
           'const double sf[2][3] = {{%.3f,%.3f,%.3f},{%.3f,%.3f,%.3f}};' % tuple(sf[0] + sf[1]), s)
s = re.sub(r'per-configuration result tables\): incl [0-9./]+ %, 1p [0-9./]+ %',
           'per-configuration result tables): incl %.1f/%.1f/%.1f %%, 1p %.1f/%.1f/%.1f %%' % tuple(100 * v for v in sf[0] + sf[1]), s)
open(p, 'w').write(s); print('  counting scale factors', sf)
PYEOF
( cd xsec_analyzer && set +u && source ./setup_xsec_analyzer.sh >/dev/null 2>&1 && set -u && root -l -b -q 'macros/total_xsec_counting.C' > ../logs/counting_review7.log 2>&1 )
command grep "sigma_count" logs/counting_review7.log
echo "==== 5. tables and figures $(date)"
python3 report/tools/rebuild_tables.py report/results/closure_summary.tsv logs/counting_review7.log > $L/rebuild_tables.log 2>&1; tail -1 $L/rebuild_tables.log | cut -c1-200
for t in ppi2bin_tables ac_rowsums_table infocontent binning_table xsec2d_figs newobs_xsec_fig; do
  python3 report/tools/$t.py > $L/$t.log 2>&1; echo "  $t rc=$?"
done
cp -p report/figures/bkgfit/newobs_xsec.pdf report/figures/
python3 xsec_analyzer/norm_manifest.py > $L/norm_manifest.log 2>&1; tail -1 $L/norm_manifest.log
( cd xsec_analyzer && bash dagostini_050.sh > ../$L/dagostini.log 2>&1 ); tail -2 $L/dagostini.log
echo "==== 6. validation summaries $(date)"
for a in "FHC5 incl" "RHCFULL incl" "FHC5 1p"; do set -- $a; python3 report/tools/altmodel_eval.py $1 $2 > $L/altmodel_$1_$2.log 2>&1; done
python3 report/tools/altmodel_eval.py V1 > $L/altmodel_V1.log 2>&1
( for a in "FHC5 incl" "RHCFULL incl" "FHC5 1p"; do set -- $a; python3 report/tools/altmodel_eval.py $1 $2 perbin; done ) > $L/altmodel_perbin.log 2>&1
for s in ensemble_tables_v1 ensemble_tables_1p ensemble_tables; do python3 report/tools/$s.py > $L/$s.log 2>&1; echo "  $s rc=$?"; done
echo "==== V19 DOWNSTREAM DONE $(date)"
