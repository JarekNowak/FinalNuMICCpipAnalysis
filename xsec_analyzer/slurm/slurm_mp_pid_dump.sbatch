#!/bin/bash
# SLURM array: Phase 1 of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md. Dump the primary track-like
# PFParticles of ONE raw overlay into the particle-classifier training sample
# (macros/mp_pid/dump_pid_tracks.C) under /data/uboone/processed/mp_pid. Overlays only.
#   sbatch --array=0-13 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_pid_manifest.list \
#          slurm/slurm_mp_pid_dump.sbatch
# The Run-4 detector-variation samples go to a separate directory, outside the training glob:
#   sbatch --array=0-26 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/mp_pid_detvar_manifest.list,OUT=/data/uboone/processed/mp_pid_detvar \
#          slurm/slurm_mp_pid_dump.sbatch
# Compile the macro once before submitting (root -l -b -q 'macros/mp_pid/dump_pid_tracks.C+'), so the
# array members do not race to rebuild it.
#SBATCH --job-name=mppid
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"; OUT="${OUT:-/data/uboone/processed/mp_pid}"; MANIFEST="${MANIFEST:?}"
line=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS='|' read -r RAW PERIOD HORN TAG <<< "$(sed -n "${line}p" "$MANIFEST")"
[ -n "${TAG:-}" ] || { echo "ERROR: no work on line $line" >&2; exit 1; }
case "$RAW" in *beamon*|*beam_on*|*onbeam*) echo "REFUSED: beam-on input $RAW"; exit 2;; esac
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
mkdir -p "$OUT"
OUTF="$OUT/pidtrk-$TAG.root"
root -l -b -q "macros/mp_pid/dump_pid_tracks.C+(\"$RAW\",\"$OUTF\",$PERIOD,$HORN,0)" 2>&1 | grep -v 'TClass::Init' | tail -5
NT=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$OUTF\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"trk\"):nullptr;printf(\"%lld\n\",t?t->GetEntries():-1LL);" 2>/dev/null | tail -1)
echo "DONE $TAG period=$PERIOD horn=$HORN tracks=$NT"
[ "${NT:--1}" -gt 0 ]
