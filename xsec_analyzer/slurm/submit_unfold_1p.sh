#!/bin/bash
# Unfold the proton-tagged extractions from the detVar-carrying universes, and harvest
# the covariance and additional-smearing matrices for the data release.
#
# Run this only after submit_univmake_1p.sh has completed. It refuses to submit for any
# extraction whose univmake is missing or older than the file_properties file that
# introduced the detector variations, because unfolding one of those would silently
# reproduce the detVar-free result this whole exercise exists to fix.
#
# Usage:
#   ./submit_unfold_1p.sh              # all ready extractions
#   ./submit_unfold_1p.sh --dry-run
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPTDIR/.." && pwd)"
RELEASE="$(cd "$REPO/../report" && pwd)/data_release"
PROC=/data/uboone/processed
PARTITION=main; DRYRUN=0

while [ $# -gt 0 ]; do
  case "$1" in
    -p) PARTITION="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

SRC="$SCRIPTDIR/univmake_1p_manifest.list"
[ -f "$SRC" ] || { echo "ERROR: $SRC not found (run submit_univmake_1p.sh first)" >&2; exit 1; }

declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
MANIFEST="$SCRIPTDIR/unfold_1p_manifest.list"
: > "$MANIFEST"
notready=0
while read -r cfg obs; do
  u="$PROC/ccpi1p_${TAG[$cfg]}_${obs}_univmake.root"
  fpm="$REPO/configs/file_properties_numi_${cfg}_w.txt"
  if [ ! -s "$u" ]; then
    echo "  NOT READY $cfg $obs (univmake missing)"; notready=$((notready+1)); continue
  fi
  if [ "$fpm" -nt "$u" ]; then
    echo "  NOT READY $cfg $obs (univmake predates the detVar config)"; notready=$((notready+1)); continue
  fi
  echo "$cfg $obs" >> "$MANIFEST"
done < "$SRC"

N=$(wc -l < "$MANIFEST")
echo "  ready: $N   not ready: $notready"
[ "$N" -gt 0 ] || { echo "ERROR: nothing ready to unfold" >&2; exit 1; }
if [ "$notready" -gt 0 ]; then
  echo "  (submitting only the ready ones; re-run this script for the rest)"
fi

mkdir -p "$RELEASE/cov"
args=( --array="0-$((N-1))" --partition="$PARTITION"
       --export="ALL,MANIFEST=${MANIFEST},REPO=${REPO},PROC=${PROC},RELEASE=${RELEASE}"
       --chdir="$SCRIPTDIR"
       "$SCRIPTDIR/slurm_unfold_1p.sbatch" )
if [ "$DRYRUN" -eq 1 ]; then echo "[dry-run] sbatch ${args[*]}"; exit 0; fi
sbatch "${args[@]}"
