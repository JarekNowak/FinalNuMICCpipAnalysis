#!/usr/bin/env bash
# dvfrac_impact.sh -- size of the fractional detector-variation treatment (SystematicsCalculator "DVfrac",
# configs/ccpi_systcalc_numi_dvfrac.conf) on every released single-pion extraction, without touching the release.
#
# Per product: the release universe file is COPIED to $OUT/univ (the extraction opens it in update mode), and the
# release extraction configuration is run twice with this checkout's UnfolderNuMI, once as released ("base", detVar
# as DV) and once with the detVar lines as DVfrac ("dvfrac"). Outputs, closure sidecars, covariance tables and logs
# go to $OUT/{base,dvfrac}/<tag>/ only; nothing under report/, logs/ or the release names in /data/uboone/processed
# is written. "base" against the release SYSTDUMP logs checks that the comparison is like for like.
#
#   bash scripts/dvfrac_impact.sh [NPAR]        (default 4 products in parallel; each unfold takes about a minute)
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROC=/data/uboone/processed
OUT=${DVFRAC_OUT:-$PROC/dvfrac_study}
NPAR=${1:-4}
cd "$REPO"; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO" REPO OUT
mkdir -p "$OUT/univ" "$OUT/cfg" "$OUT/base" "$OUT/dvfrac" "$OUT/work"

# tag | xsec config | slice config (the release lists: regen_fdfix.sh, braggfix_regen.sh, harvest_proposed_local.sh)
LIST=$OUT/products.txt; : > "$LIST"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
for cfg in fhc5 rhcfull comb; do
  t=${TAG[$cfg]}
  for o in pmu ppi costhmu costhpi thmupi; do
    if [ "$o" = ppi ]; then echo "incl_${t}_ppi3bin|configs/ccpi_xsec_config_numi_ppi3bin_${cfg}.txt|configs/ccpi_ppi_slice_config_3bin.txt"
    else echo "incl_${t}_${o}|configs/ccpi_xsec_config_numi_${o}_${cfg}.txt|configs/ccpi_${o}_slice_config_opt.txt"; fi
  done
  echo "incl_${t}_total|configs/ccpi_xsec_config_numi_total_${cfg}.txt|configs/ccpi_total_slice_config_opt.txt"
  for o in pmu ppi2bin costhmu costhpi thmupi Wpipr Whad dpt2bin dphit3bin dalphat2bin pn2bin pp thetap thpipr; do
    case "$o" in
      dpt2bin|dalphat2bin|pn2bin) sc=configs/ccpi1p_${o%2bin}_slice_config_2bin.txt ;;
      dphit3bin) sc=configs/ccpi1p_dphit_slice_config_3bin.txt ;;
      ppi2bin) sc=configs/ccpi1p_ppi_slice_config_2bin.txt ;;
      *) sc=configs/ccpi1p_${o}_slice_config.txt ;;
    esac
    echo "1p_${t}_${o}|configs/ccpi1p_xsec_config_numi_${o}_${cfg}.txt|$sc"
  done
  echo "1p_${t}_total|configs/ccpi1p_xsec_config_numi_total_${cfg}.txt|configs/ccpi1p_total_slice_config.txt"
done >> "$LIST"
echo "1p_COMB_thetap_dpt2d|configs/d2/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb.txt|configs/d2/ccpi1p_thetap_dpt_slice_config_2d.txt" >> "$LIST"
echo "incl_COMB_costhpi_costhmu2d|configs/d2/ccpi_xsec_config_numi_costhpi_costhmu_2d_comb.txt|configs/d2/ccpi_costhpi_costhmu_slice_config_2d.txt" >> "$LIST"
# DVFRAC_ONLY=<regex> restricts the run to matching tags (a test on one product first)
if [ -n "${DVFRAC_ONLY:-}" ]; then grep -E "^(${DVFRAC_ONLY})\|" "$LIST" > "$LIST.sel"; mv "$LIST.sel" "$LIST"; fi

# configurations and universe copies
while IFS='|' read -r tag xc sc; do
  for f in "$xc" "$sc"; do [ -f "$f" ] || { echo "missing $f ($tag)"; exit 1; }; done
  univ=$(awk '$1=="UnivFile"{print $2}' "$xc"); copy="$OUT/univ/$(basename "$univ")"
  if [ ! -f "$copy" ] || [ "$(stat -c %s "$copy")" != "$(stat -c %s "$univ")" ]; then cp -p "$univ" "$copy.tmp" && mv "$copy.tmp" "$copy"; fi
  sed "s|^UnivFile .*|UnivFile $copy|" "$xc" > "$OUT/cfg/${tag}_base.txt"
  sed -e "s|^UnivFile .*|UnivFile $copy|" -e "s|^SystFile .*|SystFile configs/ccpi_systcalc_numi_dvfrac.conf|" "$xc" > "$OUT/cfg/${tag}_dvfrac.txt"
done < "$LIST"
echo "$(wc -l < "$LIST") products, universe copies in $OUT/univ ($(du -sh "$OUT/univ" | cut -f1))"

one() {
  local tag=$1 sc=$2 mode
  for mode in base dvfrac; do
    local W="$OUT/work/${tag}_${mode}" D="$OUT/$mode/$tag"
    rm -rf "$W" "$D"; mkdir -p "$W/unfold_output" "$D"
    ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
    ln -s "$REPO/booster_decision_tree" "$W/booster_decision_tree"
    ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$OUT/cfg/${tag}_${mode}.txt" "$sc" "$D/xsec_${tag}.root" ) > "$D/run.log" 2>&1
    local rc=$?
    for f in "$W"/unfold_output/mat_table_cov_*.txt; do [ -f "$f" ] && cp -p "$f" "$D/$(basename "$f" | sed 's/^mat_table_//')"; done
    [ -f "$W/unfold_output/vec_table_unfolded_signal.txt" ] && cp -p "$W/unfold_output/vec_table_unfolded_signal.txt" "$D/unfolded_signal.txt"
    grep '\[SYSTDUMP\]' "$D/run.log" > "$D/systdump.txt"
    local cache; cache=$(grep -c "Reusing cached universes" "$D/run.log")
    echo "$(date +%H:%M) $tag $mode rc=$rc cov=$(ls "$D"/cov_*.txt 2>/dev/null | wc -l) cache_reused=$cache"
    rm -rf "$W"
  done
}
export -f one
echo "==== DVFRAC IMPACT START $(date) ===="
awk -F'|' '{print $1, $3}' "$LIST" | xargs -P "$NPAR" -n 2 bash -c 'one "$0" "$1"'
echo "==== DVFRAC IMPACT DONE $(date) ===="
