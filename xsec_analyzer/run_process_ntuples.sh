#!/usr/bin/env bash
# run_process_ntuples.sh
# Stage 1 of the xsec_analyzer chain: run ProcessNTuples over the raw PeLEE
# ntuples listed in configs/files_to_process_numi.txt, applying the CC1mu1piXp
# selection, and write the processed xsec-ana-*.root files to OUT_DIR.
#
# These outputs are what configs/file_properties_numi.txt points at; the
# normalisation (POT / EXT triggers) is verified there against the main
# analysis (selection/ccpi_selection.C Scale[]).
#
# Usage:
#   ./run_process_ntuples.sh                       # full run (build if needed)
#   ./run_process_ntuples.sh -j 4                  # 4 files in parallel
#   OUT_DIR=/somewhere ./run_process_ntuples.sh    # override output directory
#
# Env overrides:  OUT_DIR  SELECTION  NTUPLE_LIST

set -euo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$HERE"

OUT_DIR="${OUT_DIR:-/data/uboone/processed}"
SELECTION="${SELECTION:-CC1mu1piXp}"
NTUPLE_LIST="${NTUPLE_LIST:-configs/files_to_process_numi.txt}"
JOBS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs) JOBS="$2"; shift 2 ;;
    -h|--help) sed -n '2,18p' "$0" | sed 's/^# \?//'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# ── Environment ────────────────────────────────────────────────────────────
# setup_xsec_analyzer.sh exports XSEC_ANALYZER_DIR, PATH (+bin), ROOT_INCLUDE_PATH
# and LD_LIBRARY_PATH (+lib).  It assumes ROOT is already in your environment.
# It appends to vars that may be unset, so relax nounset while sourcing it.
set +u
source ./setup_xsec_analyzer.sh
set -u

if ! command -v root &>/dev/null; then
  echo "ERROR: ROOT not found in PATH. Source your ROOT setup first." >&2
  exit 1
fi

# ── Build if the binary is missing ─────────────────────────────────────────
if [[ ! -x bin/ProcessNTuples ]]; then
  echo "[build] bin/ProcessNTuples not found — running make ..."
  make -j"$(nproc)" bin/ProcessNTuples
fi

# ── Preflight: list exists, every input file exists ────────────────────────
[[ -f "$NTUPLE_LIST" ]] || { echo "ERROR: $NTUPLE_LIST not found" >&2; exit 1; }

mkdir -p "$OUT_DIR"
echo "[info] selection : $SELECTION"
echo "[info] ntuple list: $NTUPLE_LIST"
echo "[info] output dir : $OUT_DIR"

missing=0
while read -r path _type; do
  [[ -z "${path:-}" || "$path" == \#* ]] && continue
  [[ -f "$path" ]] || { echo "  MISSING INPUT: $path" >&2; missing=1; }
done < "$NTUPLE_LIST"
[[ $missing -eq 0 ]] || { echo "ERROR: missing input files (see above)." >&2; exit 1; }

# ── Process ────────────────────────────────────────────────────────────────
# One ProcessNTuples call per input line; output name matches the convention
# used downstream by file_properties_numi.txt: xsec-ana-<basename>.
run_one() {
  local path="$1" ftype="$2"
  local out="${OUT_DIR}/xsec-ana-$(basename "$path")"
  echo "[proc] $(basename "$path")  ($ftype)  ->  $out"
  ProcessNTuples "$path" "$ftype" "$SELECTION" "$out"
}
export -f run_one
export OUT_DIR SELECTION

mapfile -t LINES < <(grep -vE '^\s*#|^\s*$' "$NTUPLE_LIST")

if [[ "$JOBS" -gt 1 ]]; then
  printf '%s\n' "${LINES[@]}" \
    | xargs -P "$JOBS" -I{} bash -c 'run_one $0 $1' {}
else
  for line in "${LINES[@]}"; do
    # shellcheck disable=SC2086
    run_one $line
  done
fi

echo
echo "[done] Processed $(ls -1 "${OUT_DIR}"/xsec-ana-*.root 2>/dev/null | wc -l) file(s) in ${OUT_DIR}"
echo "       Next: UniverseMaker with configs/file_properties_numi.txt + a bin config."
