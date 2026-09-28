#!/usr/bin/env bash
# ens_v1.sh -- fake-data ensembles (100 members each) for the results proposed for approval on 2026-09-27
# that had none (prerequisite V1): proton-tagged theta_p and theta_pi_p in FHC, RHC and combined, and the two
# combined-only two-dimensional results, theta_p x delta p_T (proton-tagged) and cos theta_pi x cos theta_mu
# (inclusive). Same member jobs as ens_050.sh (each member rebuilds its universes from its own throw; throws
# are shared and created once, macros/throw_guard.h); the two-dimensional configs are read from configs/d2.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO/slurm"
PROC=/data/uboone/processed; LOGF=rebuild_jobs_2026-09-19.txt
RES="--cpus-per-task=2 --time=16:00:00"
for spec in "comb thetap configs 40" "comb thpipr configs 40" "fhc5 thetap configs 20" "fhc5 thpipr configs 20" "rhcfull thetap configs 20" "rhcfull thpipr configs 20" "comb thetap_dpt_2d configs/d2 40"; do
  set -- $spec
  J=$(sbatch --parsable $RES --array=0-99%$4 --export=ALL,REPO=$REPO,PROC=$PROC,CONF=$3,UNIV_THREADS=2,OBS=$2,CFG=$1 slurm_ensemble_1p.sbatch)
  echo "$J ensemble 1p $1 $2 (100 members, V1) $(date)" >> $LOGF; echo "submitted $J 1p $1 $2"
done
J=$(sbatch --parsable $RES --array=0-99%40 --export=ALL,REPO=$REPO,PROC=$PROC,CONF=configs/d2,UNIV_THREADS=2,OBS=costhpi_costhmu_2d,CFG=comb slurm_ensemble_cfg.sbatch)
echo "$J ensemble incl comb costhpi_costhmu_2d (100 members, V1) $(date)" >> $LOGF; echo "submitted $J incl comb costhpi_costhmu_2d"
