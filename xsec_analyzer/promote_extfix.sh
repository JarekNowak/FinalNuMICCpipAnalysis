#!/usr/bin/env bash
# promote_extfix.sh -- after slurm_extfix.sbatch reports 51/51 OK: move the per-run-EXT
# universes from $PROC/rebuild_extfix into the live tree (old ones backed up), then run
# the full downstream chain (re-unfold 51, harvest covariances, export release, figures).
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; RB=$PROC/rebuild_extfix; BK=$PROC/univmake_pooledext_backup_20260902
S=${FDFIX_SCRATCH:-/tmp/fdfix_work}
n_ok=$(grep -l '^OK ' slurm/extfix_*.out 2>/dev/null | wc -l)
[ "$n_ok" -eq 51 ] || { echo "only $n_ok/51 OK -- not promoting"; exit 1; }
mkdir -p "$BK"
for u in "$RB"/*_univmake.root; do b=$(basename "$u"); [ -f "$PROC/$b" ] && mv "$PROC/$b" "$BK/$b"; mv "$u" "$PROC/$b"; done
echo "promoted $(ls $PROC/*_univmake.root | wc -l) univmakes (old in $BK)"
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
FDFIX_SCRATCH=$S ./regen_fdfix.sh > ../logs/fdfix/regen_perrun.log 2>&1
bash slurm/unfold_incl_local.sh > ../logs/fdfix/harvest_incl3.log 2>&1 &
bash slurm/unfold_1p_local.sh   > ../logs/fdfix/harvest_1p3.log 2>&1 &
wait
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
L=../logs/fdfix; FIG=../report/figures
root.exe -l -b -q 'macros/export_curves.C("../report/data_release")' > $L/export_curves4.log 2>&1
root.exe -l -b -q 'macros/export_matrices.C("../report/data_release")' > $L/export_matrices4.log 2>&1
for t in FHC5 RHCFULL COMB; do cp -p $PROC/closure_hists_xsec_${t}_ppi2bin.root $PROC/closure_hists_xsec_${t}_ppi.root; root.exe -l -b -q "macros/ppi2bin_figs.C(\"$t\")" > $L/ppi2bin4_$t.log 2>&1; gt=newg4; [ $t = RHCFULL ] && gt=rhc; [ $t = COMB ] && gt=comb; root.exe -l -b -q "macros/dsigma_current.C(\"$t\",\"$gt\")" > $L/dsigma4_$t.log 2>&1; cp unfold_output/dsigma_$t.pdf $FIG/; done
for set in wtki wtki_noW incl; do root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"FHC5\",\"$set\")" > $L/ccpi1p4_$set.log 2>&1; done
for c in RHCFULL COMB; do root.exe -l -b -q "macros/dsigma_ccpi1p.C(\"$c\",\"incl\")" > $L/ccpi1p4_incl_$c.log 2>&1; done
root.exe -l -b -q 'macros/dsigma_build.C("FHC5")' > $L/dsigma_build4.log 2>&1; root.exe -l -b -q 'macros/dsigma_build_1p.C("FHC5")' > $L/dsigma_build_1p4.log 2>&1
for c in fhc5 rhcfull comb; do root.exe -l -b -q "macros/systbreak_fig.C(\"$c\")" > $L/systbreak4_$c.log 2>&1; done
for m in fhc rhc comb; do root.exe -l -b -q "macros/cutflow_yields.C(\"$m\")" > $L/cutflow4_$m.log 2>&1; cp unfold_output/cutflow_yields_$m.pdf $FIG/; done
root -l -b -q 'macros/total_xsec_counting.C' > ../logs/counting_review5.log 2>&1
cp -p unfold_output/dsigma_ccpi1p_FHC5.pdf unfold_output/dsigma_ccpi1p_noW_FHC5.pdf unfold_output/dsigma_ccpi1p_incl_*.pdf unfold_output/dsigma_build_FHC5_*.pdf unfold_output/dsigma_build_1p_FHC5_*.pdf $FIG/ 2>/dev/null
python3 check_staleness.py | sed -n 1,6p
echo "==== PROMOTE+REGEN DONE $(date) ===="
