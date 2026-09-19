#!/usr/bin/env bash
# promote_run1fix.sh -- after slurm_rebuild.sbatch (job 3351878, rebuild_run1fix.list) reports
# 34/34 OK: move the FHC5/COMB universe files rebuilt on the corrected FHC Run-1 exposure
# (2.192e20 POT / 5,748,692 triggers) from $PROC/rebuild into the live tree (old ones backed
# up), then run the full downstream chain exactly as promote_extfix.sh / regen_costhpi_release.sh
# do. RHC is untouched; the one-bin totals follow in promote_run1fix_totals.sh.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; RB=$PROC/rebuild; JOB=${RUN1FIX_JOB:-3351878}; BK=$PROC/univmake_backup_job$JOB
S=${FDFIX_SCRATCH:-$PROC/unfold_scratch/fdfix_work}
n_ok=$(grep -l '^OK ' rebuild_${JOB}_*.out 2>/dev/null | wc -l)
[ "$n_ok" -eq 34 ] || { echo "only $n_ok/34 OK -- not promoting"; exit 1; }
mkdir -p "$BK"; n=0
while read -r FAM CFG OBS; do
  case $FAM in 1p) PFX=ccpi1p;; *) PFX=ccpi;; esac
  case $CFG in fhc5) T=FHC5;; comb) T=COMB;; *) echo "unexpected cfg $CFG"; exit 1;; esac
  u="$RB/${PFX}_${T}_${OBS}_univmake.root"
  [ -s "$u" ] && [ "$(stat -c %Y "$u")" -gt "$(date -d '2026-09-19 10:00' +%s)" ] || { echo "missing or stale $u"; exit 1; }
  [ -f "$PROC/$(basename $u)" ] && mv "$PROC/$(basename $u)" "$BK/"; mv "$u" "$PROC/"; n=$((n+1))
done < slurm/rebuild_run1fix.list
echo "promoted $n univmakes (old in $BK)"
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
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
for c in RHCFULL COMB; do root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"$c\",\"incl\")" > $L/ccpi1p_run1fix_incl_$c.log 2>&1; done
root.exe -l -b -q 'macros/dsigma_build.C("FHC5")' > $L/dsigma_build_run1fix.log 2>&1; root.exe -l -b -q 'macros/dsigma_build_1p.C("FHC5")' > $L/dsigma_build_1p_run1fix.log 2>&1
for c in fhc5 rhcfull comb; do root.exe -l -b -q "macros/systbreak_fig.C(\"$c\")" > $L/systbreak_run1fix_$c.log 2>&1; done
for m in fhc rhc comb; do root.exe -l -b -q "macros/cutflow_yields.C(\"$m\")" > $L/cutflow_run1fix_$m.log 2>&1; cp unfold_output/cutflow_yields_$m.pdf $FIG/; done
root -l -b -q 'macros/total_xsec_counting.C' > ../logs/counting_review7.log 2>&1
cp -p unfold_output/dsigma_ccpi1p_FHC5.pdf unfold_output/dsigma_ccpi1p_noW_FHC5.pdf unfold_output/dsigma_ccpi1p_incl_*.pdf unfold_output/dsigma_build_FHC5_*.pdf unfold_output/dsigma_build_1p_FHC5_*.pdf $FIG/ 2>/dev/null
python3 check_staleness.py | sed -n 1,6p
echo "==== RUN1FIX PROMOTE+REGEN DONE $(date) ===="
