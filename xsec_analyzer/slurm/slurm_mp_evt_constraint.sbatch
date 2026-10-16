#!/bin/bash
#SBATCH --job-name=mpevtcon
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=04:00:00
# SLURM array task: one point of the working-point scan of the two-pion event classifier on the constrained total
# (scripts/mp_evt.py constraint). Line SLURM_ARRAY_TASK_ID + 1 of MANIFEST: MP_EVT_DIR|tag|feature set|threshold.
# Output: report/planning/multipion/phase4_constraint_<tag>[_<set>]_t<threshold>.{md,json}; summary with
# scripts/mp_evt_wp_scan_summary.py.
#   sbatch --array=0-41 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_evt_constraint_manifest.list slurm/slurm_mp_evt_constraint.sbatch
set -u
REPO="${REPO:?}"; MANIFEST="${MANIFEST:?}"
IFS='|' read -r DIR TAG SET T <<< "$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$MANIFEST")"
[ -n "${T:-}" ] || { echo "ERROR: no work on line $(( SLURM_ARRAY_TASK_ID + 1 ))" >&2; exit 1; }
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
echo "host=$(hostname) dir=$DIR tag=$TAG set=$SET threshold=$T"
MP_EVT_DIR="$DIR" python3 scripts/mp_evt.py constraint --tag "$TAG" --set "$SET" --threshold "$T"; rc=$?
echo "DONE $TAG $SET $T rc=$rc"; exit $rc
