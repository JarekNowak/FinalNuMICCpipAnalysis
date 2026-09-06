#!/usr/bin/env bash
# regen_costhpi_release.sh -- adopt the forward-split five-bin cos(theta_pi) binning (2026-09-06):
# install the regenerated generator predictions (FHC/RHC/combined FTE files; the other observables
# were verified identical to the previous files), promote the three rebuilt cos(theta_pi) universe
# files (old ones backed up), then run the full downstream chain exactly as promote_extfix.sh does.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; RB=$PROC/rebuild_extfix; BK=$PROC/univmake_backup_costhpi4bin_20260906
S=${FDFIX_SCRATCH:-/tmp/fdfix_work}; GP=../generator_predictions
mkdir -p "$BK" "$GP/newg4/fte_backup_costhpi4bin"
for g in genie gibuu neut nuwro; do for pair in "newg4:c5" "rhc:c5_rhc" "comb:c5_comb"; do
  old="$GP/newg4/${g}_${pair%%:*}_fte.root"; new="$GP/costhpi5/${g}_${pair##*:}_fte.root"
  [ -f "$new" ] || { echo "missing $new"; exit 1; }; cp -p "$old" "$GP/newg4/fte_backup_costhpi4bin/"; cp "$new" "$old"; done; done
echo "installed 12 prediction files (old in newg4/fte_backup_costhpi4bin)"
for T in FHC5 RHCFULL COMB; do u="$RB/ccpi_${T}_costhpi_univmake.root"; [ -s "$u" ] || { echo "missing $u"; exit 1; }
  mv "$PROC/ccpi_${T}_costhpi_univmake.root" "$BK/"; mv "$u" "$PROC/"; done
echo "promoted 3 cos(theta_pi) univmakes (old in $BK)"
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
echo "==== COSTHPI5 RELEASE DONE $(date) ===="
