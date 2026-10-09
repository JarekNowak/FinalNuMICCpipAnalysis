#!/bin/bash
# unfold_only_manifest.sh -- the unfolding half of slurm/slurm_extract.sbatch for every line of a manifest
# ("XSEC_CONFIG BIN_CONFIG SLICE_CONFIG"), on the EXISTING universe files: for a change that touches only the
# covariance or the unfolding (version 1.8: detector variations as DVfrac), where slurm_extract.sbatch would rebuild
# identical universes for hours. Same outputs as slurm_extract.sbatch (xsec_<B>.root with its closure sidecars,
# unfold_<B>.log, plots_<B>/ next to the universe file); the previous sidecars and log of each line are first copied
# to <dir>/backup_<TAG>/. Runs locally, NPAR lines at a time (each unfold takes seconds to a minute).
#   bash slurm/unfold_only_manifest.sh slurm/braggfix_altmodel_manifest.list v1.7 [NPAR]
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$REPO"
MANIFEST=${1:?manifest}; TAG=${2:?backup tag}; NPAR=${3:-6}
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO" REPO TAG
one() {
  local XC=$1 SLICE=$2
  local U; U=$(awk '$1=="UnivFile"{print $2}' "$XC")
  [ -s "$U" ] || { echo "MISSING universe $U ($XC)"; return 1; }
  local D B W; D=$(dirname "$U"); B=$(basename "$U" _univmake.root); W="$D/work_${B}"
  mkdir -p "$D/backup_${TAG}"
  for f in "$D"/closure_hists_xsec_"$B".root "$D"/closure_hists_all_xsec_"$B".root "$D/unfold_${B}.log"; do
    [ -f "$f" ] && [ ! -f "$D/backup_${TAG}/$(basename "$f")" ] && cp -p "$f" "$D/backup_${TAG}/"
  done
  rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && nice -n 10 "$REPO/bin/UnfolderNuMI" "$REPO/$XC" "$REPO/$SLICE" "$D/xsec_${B}.root" ) > "$D/unfold_${B}.log" 2>&1
  local rc=$?; mkdir -p "$D/plots_${B}"; cp "$W"/unfold_output/*.pdf "$D/plots_${B}/" 2>/dev/null; rm -rf "$W"
  local cache; cache=$(grep -c "Reusing cached universes" "$D/unfold_${B}.log")
  [ $rc -eq 0 ] && echo "OK $B cache_reused=$cache sigma_int=$(grep 'SYSTDUMP] sigma_int' "$D/unfold_${B}.log" | tail -1 | awk '{print $3}')" \
    || echo "FAILED $B rc=$rc"
}
export -f one
awk 'NF>=3 && $1!~/^#/ {print $1, $3}' "$MANIFEST" | xargs -P "$NPAR" -n 2 bash -c 'one "$0" "$1"'
