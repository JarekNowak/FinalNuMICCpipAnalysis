#!/bin/bash
# SLURM: one variant of the multi-pion particle classifier for the detector-robustness study (follow-up of Phase 2
# of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md): scripts/mp_pid_train.py --unmatched --common with the variant's
# options, then its paired detector sensitivity on the odd events of the Run-4 variations (scripts/mp_pid_detsens.py).
# Reports in report/planning/multipion/detrobust/, the model in /data/uboone/processed/mp_pid/.
#   sbatch --export=ALL,REPO=$PWD,TAG=aug50,OPTS="--aug 0.5" slurm/slurm_mp_pid_train.sbatch
#SBATCH --job-name=mppidtr
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.out
set -u
REPO="${REPO:?}"; TAG="${TAG:?}"; OPTS="${OPTS:-}"
cd "$REPO" || exit 1
python3 scripts/mp_pid_train.py --unmatched --common --tag "$TAG" $OPTS || exit 1
python3 scripts/mp_pid_detsens.py model /data/uboone/processed/mp_pid/model_xgb_unmatched_common_${TAG}.json
