#!/usr/bin/env bash
# promote_run1fix_totals.sh -- after slurm_rebuild.sbatch (rebuild_run1fix_total.list) reports
# 4/4 OK: promote the four one-bin total universe files (incl/1p x FHC5/COMB) rebuilt on the
# corrected FHC Run-1 exposure, re-unfold them into the live sidecars and systdump logs, and
# re-export data_release/total_xsec.tsv. RHC totals are untouched.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; RB=$PROC/rebuild; JOB=${RUN1FIX_TOTAL_JOB:?set RUN1FIX_TOTAL_JOB}; BK=$PROC/univmake_backup_job$JOB
DUMP=../logs/systdump; LOG=../logs/fdfix; S=$PROC/unfold_scratch/totals_run1fix
n_ok=$(grep -l '^OK ' rebuild_${JOB}_*.out 2>/dev/null | wc -l)
[ "$n_ok" -eq 4 ] || { echo "only $n_ok/4 OK -- not promoting"; exit 1; }
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
mkdir -p "$BK" "$S/unfold_output"; ln -sfn "$REPO/configs" "$S/configs"; ln -sfn "$REPO/bin" "$S/bin"; ln -sfn "$REPO/lib" "$S/lib"
while read -r FAM CFG OBS; do
  case $FAM in 1p) PFX=ccpi1p; SL=configs/ccpi1p_total_slice_config.txt; DP=1p_; FP=ccpi1p_;; *) PFX=ccpi; SL=configs/ccpi_total_slice_config_opt.txt; DP=""; FP="";; esac
  case $CFG in fhc5) T=FHC5;; comb) T=COMB;; *) echo "unexpected cfg $CFG"; exit 1;; esac
  u="$RB/${PFX}_${T}_total_univmake.root"
  [ -s "$u" ] && [ "$(stat -c %Y "$u")" -gt "$(date -d '2026-09-19 18:00' +%s)" ] || { echo "missing or stale $u"; exit 1; }
  [ -f "$PROC/$(basename $u)" ] && mv "$PROC/$(basename $u)" "$BK/"; mv "$u" "$PROC/"
  ( cd "$S" && rm -f unfold_output/plot_step[1-4]_*.pdf && bin/UnfolderNuMI "configs/${PFX}_xsec_config_numi_total_${CFG}.txt" "$SL" "$PROC/xsec_${FP}${T}_total.root" ) > "$LOG/${DP}${CFG}_total.raw" 2>&1 || echo "  FAIL $FAM $CFG total"
  cp "$LOG/${DP}${CFG}_total.raw" "$DUMP/${DP}${CFG}_total.raw"; grep '\[SYSTDUMP\]' "$LOG/${DP}${CFG}_total.raw" > "$DUMP/${DP}${CFG}_total.dump"
  echo "  $FAM $CFG total sigma_int=$(awk '$2=="sigma_int"{print $3}' "$DUMP/${DP}${CFG}_total.dump")"
done < slurm/rebuild_run1fix_total.list
rm -rf "$S"
python3 ../report/tools/total_export.py
echo "==== RUN1FIX TOTALS DONE $(date) ===="
