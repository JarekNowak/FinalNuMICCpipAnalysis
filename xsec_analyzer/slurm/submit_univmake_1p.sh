#!/bin/bash
# Rebuild the proton-tagged univmakes with detector variations included, on SLURM.
#
# The proton-tagged detVar samples were processed into /data/uboone/processed/w/ but were
# never listed in the _w file_properties, and ccpi1p_systcalc_numi.conf had no DV block,
# so all 33 proton-tagged extractions carry no detector systematic at all. This rebuilds
# the universes so they do.
#
# Only the univmakes actually used by published results are built (33), not every one
# defined (49) -- the remainder are superseded fine-binned variants.
#
# Usage:
#   ./submit_univmake_1p.sh                # all 33
#   ./submit_univmake_1p.sh -c fhc5        # one configuration
#   ./submit_univmake_1p.sh -n 2           # first N (test)
#   ./submit_univmake_1p.sh --dry-run
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPTDIR/.." && pwd)"
PROC=/data/uboone/processed
PARTITION=main; DRYRUN=0; ONLYCFG=""; NTASK=0

while [ $# -gt 0 ]; do
  case "$1" in
    -c) ONLYCFG="$2"; shift 2 ;;
    -n) NTASK="$2"; shift 2 ;;
    -p) PARTITION="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

# Derive the work list from the published results rather than from a hand-kept list, so
# it cannot drift out of step with what is actually reported.
RESULTS="$REPO/../report/current_results.tsv"
[ -f "$RESULTS" ] || { echo "ERROR: $RESULTS not found" >&2; exit 1; }

MANIFEST="$SCRIPTDIR/univmake_1p_manifest.list"
# The published "ppi" row is the FIVE-bin scheme, which the analysis note withdraws
# (four of five bins fall below the 68% diagonal criterion). A complete two-bin chain
# exists -- ccpi1p_ppi_bin_config_2bin.txt with the ppi2bin xsec and slice configs -- so
# the rebuild targets that rather than reproducing a binning we do not intend to publish.
awk -F'\t' 'NR>1 && $1=="1p" {obs=$2; if (obs=="ppi") obs="ppi2bin"; print $3, obs}' \
  "$RESULTS" | sort -u > "$MANIFEST"
[ -n "$ONLYCFG" ] && { grep "^${ONLYCFG} " "$MANIFEST" > "$MANIFEST.t" && mv "$MANIFEST.t" "$MANIFEST"; }
[ "$NTASK" -gt 0 ] && { head -n "$NTASK" "$MANIFEST" > "$MANIFEST.t" && mv "$MANIFEST.t" "$MANIFEST"; }

# Refuse to submit if any task's inputs are missing -- a whole-node job that dies on a
# missing config has cost a node-hour to discover something checkable in a second.
missing=0
while read -r cfg obs; do
  case "$obs" in
    *2bin) bin="$REPO/configs/ccpi1p_${obs%2bin}_bin_config_2bin.txt" ;;
    *)     bin="$REPO/configs/ccpi1p_${obs}_bin_config.txt" ;;
  esac
  for f in "$REPO/configs/file_properties_numi_${cfg}_w.txt" "$bin"; do
    [ -f "$f" ] || { echo "  MISSING $f (for $cfg $obs)"; missing=1; }
  done
done < "$MANIFEST"
[ "$missing" -eq 0 ] || { echo "ERROR: inputs missing, not submitting" >&2; exit 1; }

# And that the configs actually carry detVar -- otherwise this rebuild achieves nothing.
for cfg in $(awk '{print $1}' "$MANIFEST" | sort -u); do
  n=$(grep -vE '^[[:space:]]*#' "$REPO/configs/file_properties_numi_${cfg}_w.txt" \
      | awk '$3 ~ /^detVar/' | wc -l)
  [ "$n" -gt 0 ] || { echo "ERROR: $cfg _w config has no detVar lines" >&2; exit 1; }
  echo "  $cfg: $n detVar entries"
done
grep -q '^detVar' "$REPO/configs/ccpi1p_systcalc_numi.conf" \
  || { echo "ERROR: ccpi1p_systcalc_numi.conf has no DV block" >&2; exit 1; }

N=$(wc -l < "$MANIFEST")
[ "$N" -gt 0 ] || { echo "ERROR: manifest empty" >&2; exit 1; }
echo "  tasks: $N   partition: $PARTITION   one whole node each"
sed 's/^/    /' "$MANIFEST"

args=( --array="0-$((N-1))" --partition="$PARTITION"
       --export="ALL,MANIFEST=${MANIFEST},REPO=${REPO},PROC=${PROC}"
       --chdir="$SCRIPTDIR"
       "$SCRIPTDIR/slurm_univmake_1p.sbatch" )
if [ "$DRYRUN" -eq 1 ]; then echo "[dry-run] sbatch ${args[*]}"; exit 0; fi
sbatch "${args[@]}"
