#!/bin/bash
#SBATCH --job-name=softpion
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/data/uboone/processed/mp2pi_2f/softpion/logs/%x_%A_%a.out
#SBATCH --error=/data/uboone/processed/mp2pi_2f/softpion/logs/%x_%A_%a.out
# SLURM array task: one overlay of the two-fold production for scripts/mp2pi_softpion.py (per-universe sums of the
# signal and the sub-threshold class by softest-pion momentum), into /data/uboone/processed/mp2pi_2f/softpion/.
#   mkdir -p /data/uboone/processed/mp2pi_2f/softpion/logs
#   sbatch --array=0-$(( $(python3 scripts/mp2pi_softpion.py n) - 1 )) --export=ALL,REPO=$PWD slurm/slurm_mp2pi_softpion.sbatch
set -u
REPO="${REPO:?}"; cd "$REPO" || exit 1
python3 scripts/mp2pi_softpion.py build1 "$SLURM_ARRAY_TASK_ID"
