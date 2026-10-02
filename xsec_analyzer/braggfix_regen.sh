#!/usr/bin/env bash
# braggfix_regen.sh -- downstream chain after braggfix_promote_univ.sh (release rebuilt without the
# Bragg-pion cut): the 51 extractions with covariances, exports, figures, counting and cut-flows
# (regen_release.sh; the cut-flows read processed/cf and processed/cf_1p, links into the promoted trees
# after braggfix_cf_links.sh), the six one-bin totals and their export, and the proposed results
# (theta_p, theta_pi_p and the two combined 2D results) with their covariances.
#   ./braggfix_regen.sh > ../logs/braggfix/regen.log 2>&1
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; DUMP=../logs/systdump; LOG=../logs/fdfix
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
[ -L "$PROC/cf/xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root" ] || { echo "run braggfix_cf_links.sh first"; exit 1; }
echo "==== one-bin totals $(date)"
S=$PROC/unfold_scratch/totals_braggfix; rm -rf "$S"; mkdir -p "$S/unfold_output"
ln -sfn "$REPO/configs" "$S/configs"; ln -sfn "$REPO/bin" "$S/bin"; ln -sfn "$REPO/lib" "$S/lib"
for FAM in incl 1p; do for CFG in fhc5 rhcfull comb; do
  case $FAM in 1p) PFX=ccpi1p; SL=configs/ccpi1p_total_slice_config.txt; DP=1p_; FP=ccpi1p_;; *) PFX=ccpi; SL=configs/ccpi_total_slice_config_opt.txt; DP=""; FP="";; esac
  case $CFG in fhc5) T=FHC5;; rhcfull) T=RHCFULL;; comb) T=COMB;; esac
  ( cd "$S" && rm -f unfold_output/plot_step[1-4]_*.pdf && bin/UnfolderNuMI "configs/${PFX}_xsec_config_numi_total_${CFG}.txt" "$SL" "$PROC/xsec_${FP}${T}_total.root" ) > "$LOG/${DP}${CFG}_total.raw" 2>&1 || echo "  FAIL $FAM $CFG total"
  cp "$LOG/${DP}${CFG}_total.raw" "$DUMP/${DP}${CFG}_total.raw"; grep '\[SYSTDUMP\]' "$LOG/${DP}${CFG}_total.raw" > "$DUMP/${DP}${CFG}_total.dump"
  echo "  $FAM $CFG total sigma_int=$(awk '$2=="sigma_int"{print $3}' "$DUMP/${DP}${CFG}_total.dump")"
done; done
rm -rf "$S"
python3 ../report/tools/total_export.py
echo "==== proposed results $(date)"
bash slurm/harvest_proposed_local.sh
# regen_release.sh last: its exports and figures must postdate the totals and the proposed results
echo "==== regen_release $(date)"
bash regen_release.sh > $LOG/regen_release_braggfix.log 2>&1; tail -n 14 $LOG/regen_release_braggfix.log
echo "==== BRAGGFIX REGEN DONE $(date) ===="
