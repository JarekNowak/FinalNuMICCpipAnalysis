#!/bin/bash
#SBATCH --job-name=mp2fluxw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2500M
#SBATCH --time=12:00:00
#SBATCH --output=/data/uboone/processed/mp2pi/detvar_fluxw/logs/univ.log
# scripts/mp2pi_fluxw_univ.sh on a batch node (univmake on the 18 reweighted detector-variation samples, merge,
# extractions).   sbatch slurm/slurm_mp2pi_fluxw.sbatch
cd /home/t2k/nowak/MicroBooNE/mp_phase1_wt/xsec_analyzer || exit 1
echo "host=$(hostname)"
UNIV_THREADS=4 bash scripts/mp2pi_fluxw_univ.sh
