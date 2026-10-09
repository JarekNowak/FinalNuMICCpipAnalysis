#!/bin/bash
#SBATCH --job-name=dvfluxw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500M
#SBATCH --time=01:30:00
#SBATCH --output=/data/uboone/processed/mp2pi/detvar_fluxw/logs/fluxw_%a.log
# SLURM array task: reweighted copy of one Run-4 detector-variation sample (line SLURM_ARRAY_TASK_ID of
# slurm/mp2pi_detvar_manifest.list) with the new / old NuMI flux weight, scripts/detvar_fluxw.py.
#   mkdir -p /data/uboone/processed/mp2pi/detvar_fluxw/logs
#   sbatch --array=0-17 --export=ALL,REPO=$PWD slurm/slurm_detvar_fluxw.sbatch
set -u
REPO="${REPO:?}"
cd "$REPO" || exit 1
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID"
python3 scripts/detvar_fluxw.py "$SLURM_ARRAY_TASK_ID"
