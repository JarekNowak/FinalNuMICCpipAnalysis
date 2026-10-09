#!/bin/bash
#SBATCH --job-name=mcscalf
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500M
#SBATCH --time=10:00:00
#SBATCH --output=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mcscal/logs/%x_%A_%a.out
# SLURM array task: one line "in|out|map" of slurm/mcscal_files.list (scripts/mcs_cal_setup.py): a copy of a processed
# file with the MCS momentum of the simulated exiting muon candidate calibrated (macros/mcs_cal_apply.C). The macro is
# compiled once beforehand and loaded as a library, so concurrent tasks never recompile it:
#   root -l -b -e 'gSystem->AddIncludePath("-I'$PWD'/include")' -q -e '.L macros/mcs_cal_apply.C+'
#   mkdir -p /data/uboone/processed/mcscal/logs
#   sbatch --array=0-211%25 --export=ALL,REPO=$PWD slurm/slurm_mcscal_files.sbatch
set -u
REPO="${REPO:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r IN OUT MAP <<< "$(sed -n "${line}p" "$REPO/slurm/mcscal_files.list")"
[ -n "${MAP:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
[ -s "$OUT" ] && { echo "SKIP existing $OUT"; exit 0; }
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
[ -f macros/mcs_cal_apply_C.so ] || { echo "ERROR: compile macros/mcs_cal_apply.C first" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
echo "host=$(hostname) in=$IN map=$MAP"
root -l -b -e 'gSystem->Load("macros/mcs_cal_apply_C.so")' -q -e "mcs_cal_apply(\"$IN\",\"$OUT\",\"$MAP\")" 2>&1 | grep -v '^$'
[ -s "$OUT" ] || { echo "FAIL no output" >&2; exit 1; }
nin=$(root -l -b -q -e "TFile f(\"$IN\"); cout << ((TTree*)f.Get(\"stv_tree\"))->GetEntries() << endl;" 2>/dev/null | tail -1)
nout=$(root -l -b -q -e "TFile f(\"$OUT\"); cout << ((TTree*)f.Get(\"stv_tree\"))->GetEntries() << endl;" 2>/dev/null | tail -1)
[ "$nin" = "$nout" ] || { echo "FAIL entries $nin -> $nout" >&2; command rm -f -- "$OUT"; exit 1; }
echo "DONE $OUT entries=$nout"
