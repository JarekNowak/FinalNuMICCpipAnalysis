#!/usr/bin/env bash
# ens_050.sh -- fake-data ensembles (100 members each) for the binnings adopted under the 0.50 migration
# criterion (2026-09-26): FHC p_mu (six bins), FHC and RHC cos theta_mu (eight bins), combined p_pi (three
# regions) and the proton-tagged FHC W_had (two regions). The one-bin totals and the proton-tagged two-bin
# delta p_T and p_n keep their binnings and their ensembles. Each member rebuilds its universes from its own
# throw and needs none of the release universe files, so the members read the staged configs
# (CONF=configs/adopt050, identical to what promote_050.sh installs) and can run before the promotion.
# Outputs of the previous binnings with the same tags are moved to ens/pre050_20260926/ first. The arrays
# share each member's throw, which is created once and never replaced (macros/throw_guard.h, ens/throws/).
# Two CPUs and two universe threads per member: a four-thread build keeps only about 2.6 cores busy, so on a
# full cluster two-CPU members finish more work per hour. The combined array is submitted first with the
# larger throttle because its members take twice as long.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO/slurm"
PROC=/data/uboone/processed; E=$PROC/ens; A=$E/pre050_20260926; mkdir -p "$A"
# once per tag (marker file): a second pass would move a later run's files onto the archived ones
for tag in fhc5_pmu fhc5_costhmu rhcfull_costhmu 1p_fhc5_Whad; do
  [ -e "$A/.archived_$tag" ] && { echo "already archived: $tag"; continue; }
  n=0; for f in "$E"/*_${tag}_t[0-9]*; do [ -e "$f" ] || continue; mv -n "$f" "$A/"; n=$((n+1)); done
  touch "$A/.archived_$tag"; echo "archived $n files of $tag"
done
LOGF=rebuild_jobs_2026-09-19.txt; CONF=configs/adopt050
RES="--cpus-per-task=2 --time=16:00:00"; EXP="ALL,REPO=$REPO,PROC=$PROC,CONF=$CONF,UNIV_THREADS=2"
for spec in "comb ppi3bin 50" "fhc5 pmu 20" "fhc5 costhmu 20" "rhcfull costhmu 20"; do
  set -- $spec
  J=$(sbatch --parsable $RES --array=0-99%$3 --export=$EXP,OBS=$2,CFG=$1 slurm_ensemble_cfg.sbatch)
  echo "$J ensemble $1 $2 (100 members, binnings of the 0.50 criterion) $(date)" >> $LOGF; echo "submitted $J $1 $2"
done
J=$(sbatch --parsable $RES --array=0-99%20 --export=$EXP,OBS=Whad,CFG=fhc5 slurm_ensemble_1p.sbatch)
echo "$J ensemble 1p fhc5 Whad (100 members, two regions under the 0.50 criterion) $(date)" >> $LOGF; echo "submitted $J 1p fhc5 Whad"
