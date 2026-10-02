#!/usr/bin/env bash
# braggfix_cf_links.sh -- after braggfix_promote_trees.sh: the cut-flow inputs processed/cf (inclusive) and
# processed/cf_1p (proton-tagged) become links into the promoted trees. They were separate reprocesses only
# because the live trees predated the full-CV-weight cut-flow counters (2026-09-02); the promoted trees are
# made by the same code from the same raw files, so their content is what a cf reprocess would produce
# (and the cf_1p dirt, previously processed with CC1mu1piXp, now carries the CC1mu1pi1p counters). The old
# directories move to processed/pre_braggfix/.
set -euo pipefail
PROC=/data/uboone/processed; BK=$PROC/pre_braggfix
[ -d "$BK" ] || { echo "run braggfix_promote_trees.sh first"; exit 1; }
declare -A ALIAS=(
  [numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana]=dirt_fhc_run4c
  [numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana]=dirt_fhc_run4d
  [run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns]=dirt_fhc_run5
  [neutrinoselection_filt_run3b_dirt_overlay]=dirt_rhc_run3b
  [run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana]=dirt_rhc_run4a
  [numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana]=dirt_rhc_run4b
  [prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot]=prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot )
for pair in "cf:$PROC" "cf_1p:$PROC/w"; do
  IFS=: read -r d tree <<< "$pair"
  [ -e "$BK/$d" ] && { echo "REFUSING: $BK/$d exists"; exit 1; }
  names=$(cd "$PROC/$d" && ls xsec-ana-*.root)
  mv "$PROC/$d" "$BK/$d"; mkdir "$PROC/$d"
  for f in $names; do
    tag=${f#xsec-ana-}; tag=${tag%.root}; tgt="$tree/xsec-ana-${ALIAS[$tag]:-$tag}.root"
    [ -f "$tgt" ] || { echo "missing target $tgt for $d/$f"; exit 1; }
    ln -s "$tgt" "$PROC/$d/$f"
  done
  echo "$d: $(ls "$PROC/$d" | wc -l) links into $tree"
done
