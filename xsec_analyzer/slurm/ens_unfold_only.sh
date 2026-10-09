#!/usr/bin/env bash
# ens_unfold_only.sh TAG NPAR CFG:OBS [CFG:OBS ...] -- re-unfold every member of the given release ensembles on its
# EXISTING universe file, for a change that touches only the covariance or the unfolding (version 1.8: detector
# variations as DVfrac). One UnfolderNuMI per member, in its kept work directory ens/statcov_<CFG>_<OBS>_t<T> and with
# the member's output name ens/xsec_<CFG>_<OBS>_t<T>.root, refreshes both what the member job wrote (the closure
# sidecars in ens/, read by report/tools/ensemble_stat_pulls.py) and what slurm/ens_statcov.sh wrote (the DataStats
# and total covariance tables). The previous statcov directory, sidecars and log of each member are first moved to
# ens/backup_<TAG>/. CFG and OBS as in slurm/ens_statcov.sh.
#   bash slurm/ens_unfold_only.sh v1.7 6 fhc5:pmu fhc5:costhmu ... > ../logs/v18/ens_unfold.log 2>&1
set -u
TAG=${1:?backup tag}; NPAR=${2:?parallel}; shift 2
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; PROC=/data/uboone/processed; E=$PROC/ens
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO" REPO E TAG
mkdir -p "$E/backup_$TAG"
slice() {
  local CFG=$1 OBS=$2
  if [[ "$CFG" == 1p_* ]]; then
    case "$OBS" in *2bin) echo configs/ccpi1p_${OBS%2bin}_slice_config_2bin.txt;; *3bin) echo configs/ccpi1p_${OBS%3bin}_slice_config_3bin.txt;; *_2d) echo configs/d2/ccpi1p_${OBS%_2d}_slice_config_2d.txt;; total) echo configs/ccpi1p_total_slice_config.txt;; *) echo configs/ccpi1p_${OBS}_slice_config.txt;; esac
  else
    case "$OBS" in ppi2bin) echo configs/ccpi_ppi_slice_config_2bin.txt;; ppi3bin) echo configs/ccpi_ppi_slice_config_3bin.txt;; total) echo configs/ccpi_total_slice_config_opt.txt;; *_2d) echo configs/d2/ccpi_${OBS%_2d}_slice_config_2d.txt;; *) echo configs/ccpi_${OBS}_slice_config_opt.txt;; esac
  fi
}
member() {
  local N=$1 SL=$2 T=$3
  local XC="$E/xsec_${N}_t${T}.txt" W="$E/statcov_${N}_t${T}" B="$E/backup_$TAG"
  [ -f "$XC" ] || return 0
  local U; U=$(awk '$1=="UnivFile"{print $2}' "$XC"); [ -s "$U" ] || { echo "MISSING universe $N t$T"; return 0; }
  [ -d "$W" ] && [ ! -e "$B/statcov_${N}_t${T}" ] && mv "$W" "$B/"
  for f in "$E/closure_hists_xsec_${N}_t${T}.root" "$E/closure_hists_all_xsec_${N}_t${T}.root" "$E/unfold_${N}_t${T}.log"; do
    [ -f "$f" ] && [ ! -f "$B/$(basename "$f")" ] && cp -p "$f" "$B/"
  done
  rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && nice -n 10 "$REPO/bin/UnfolderNuMI" "$XC" "$REPO/$SL" "$E/xsec_${N}_t${T}.root" ) > "$W/run.log" 2>&1
  local rc=$?
  cp -p "$W/run.log" "$E/unfold_${N}_t${T}.log"
  if [ $rc -ne 0 ] || [ ! -f "$W/unfold_output/mat_table_cov_DataStats.txt" ]; then echo "FAIL $N t$T rc=$rc"
  elif ! grep -q "Reusing cached universes" "$W/run.log"; then echo "REBUILT-CACHE $N t$T"; fi
}
export -f member
for pair in "$@"; do
  CFG=${pair%%:*}; OBS=${pair#*:}; N="${CFG}_${OBS}"; SL=$(slice "$CFG" "$OBS")
  [ -f "$REPO/$SL" ] || { echo "missing slice config $SL"; continue; }
  echo "== $N ($SL) $(date +%H:%M)"
  seq 1 100 | xargs -P "$NPAR" -I{} bash -c "member '$N' '$SL' {}"
  n=$(ls "$E"/statcov_${N}_t*/unfold_output/mat_table_cov_DataStats.txt 2>/dev/null | wc -l)
  echo "   $N: $n members with tables $(date +%H:%M)"
done
echo "==== ENS UNFOLD DONE $(date)"
