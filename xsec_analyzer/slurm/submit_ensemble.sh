#!/usr/bin/env bash
# Fake-data ensemble: N independent Poisson throws through the full chain, to measure
# pull widths and interval coverage. The note's closure uses ONE throw and therefore
# cannot say whether the quoted uncertainties are correctly sized.
#   ./submit_ensemble.sh [-n N] [-o OBS] [--dry-run]
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPTDIR/.." && pwd)"
PROC=/data/uboone/processed; N=32; OBS=pmu; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    -n) N="$2"; shift 2 ;; -o) OBS="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;; *) echo "unknown: $1" >&2; exit 1 ;;
  esac
done
for f in "$REPO/configs/ccpi_${OBS}_bin_config_opt.txt" \
         "$REPO/configs/ccpi_${OBS}_slice_config_opt.txt" \
         "$REPO/configs/ccpi_xsec_config_numi_${OBS}_fhc5.txt" \
         "$REPO/macros/throw_ensemble_fhc.C"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
# the released fake data must be backed up, since the ensemble writes alongside it
[ -d "$PROC/ensemble_backup" ] || { echo "ERROR: no ensemble_backup; refusing" >&2; exit 1; }
mkdir -p "$PROC/ens"
echo "  throws: $N   observable: $OBS   one whole node each"
args=( --array="0-$((N-1))" --partition=main
       --export="ALL,REPO=${REPO},PROC=${PROC},OBS=${OBS}"
       --chdir="$SCRIPTDIR" "$SCRIPTDIR/slurm_ensemble.sbatch" )
if [ "$DRY" -eq 1 ]; then echo "[dry-run] sbatch ${args[*]}"; exit 0; fi
sbatch "${args[@]}"
