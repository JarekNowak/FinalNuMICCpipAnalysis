#!/bin/bash
#SBATCH --job-name=ensunf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500M
#SBATCH --time=06:00:00
#SBATCH --output=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/v18/ens_unfold_%a.log
# SLURM array task: slurm/ens_unfold_only.sh for ONE ensemble (CFG:OBS = word SLURM_ARRAY_TASK_ID+1 of PAIRS),
# members one after another.
#   sbatch --array=0-15 --export=ALL,PAIRS="fhc5:pmu ...",BTAG=v1.7 slurm/slurm_ens_unfold_only.sbatch
set -u
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO" || exit 1
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
read -r -a P <<< "${PAIRS:?}"
echo "host=$(hostname) pair=${P[$SLURM_ARRAY_TASK_ID]}"
bash slurm/ens_unfold_only.sh "${BTAG:?}" 1 "${P[$SLURM_ARRAY_TASK_ID]}"
