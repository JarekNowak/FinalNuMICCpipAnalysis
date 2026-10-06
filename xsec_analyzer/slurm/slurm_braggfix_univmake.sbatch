#!/bin/bash
# SLURM array: ONE universe file of the release, rebuilt after the Bragg-pion cut was dropped (commit
# 1fb6be9), with the exact bin configuration recovered from the live file
# (slurm/make_braggfix_univ_manifest.py) and the release file-properties list, which by then points at the
# promoted processed trees. Output under the live file name in $PROC/rebuild_braggfix/; nothing live is
# touched (promotion is a separate step). Checks that the reco and true bin specifications of the new
# file equal those of the live one.
#   sbatch --chdir=$PWD/slurm --array=0-64 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/braggfix_univ_manifest.list \
#          slurm/slurm_braggfix_univmake.sbatch
#SBATCH --job-name=bfunivmake
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; PROC="${PROC:-/data/uboone/processed}"; MANIFEST="${MANIFEST:?}"
RB="$PROC/rebuild_braggfix"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r LIVE FPM BINS <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${BINS:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"; export UNIV_THREADS="${UNIV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
mkdir -p "$RB"
for f in "$FPM" ${BINS//;/ }; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
grep -q 'ext_run' "$FPM" || { echo "ERROR: $FPM is not the per-run EXT configuration" >&2; exit 1; }
[[ "$BINS" == *";"* ]] && { echo "ERROR: several directories per file are not handled" >&2; exit 1; }
U="$RB/$(basename "$LIVE")"; rm -f "$U"
echo "host=$(hostname) live=$LIVE fpm=$FPM bin=$BINS"
./bin/univmake "$FPM" "$BINS" "$U" "$FPM" > "$RB/$(basename "$LIVE" .root).log" 2>&1
[ -s "$U" ] || { echo "FAILED univmake $LIVE" >&2; exit 1; }
python3 - "$LIVE" "$U" <<'PY' || exit 1
import sys, uproot
a, b = uproot.open(sys.argv[1]), uproot.open(sys.argv[2])
ok = True
for n in {k.split(';')[0] for k in b.keys(recursive=False)}:
    for s in ('reco_bin_spec', 'true_bin_spec', 'sel_for_categ'):
        if a[n][s] != b[n][s]: print('BIN SPEC DIFFERS', n, s); ok = False
print('bin specs', 'identical' if ok else 'DIFFER'); sys.exit(0 if ok else 1)
PY
echo "OK $(basename "$U") $(stat -c %s "$U")"
