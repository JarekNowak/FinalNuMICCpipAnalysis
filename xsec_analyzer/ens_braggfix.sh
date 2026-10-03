#!/usr/bin/env bash
# ens_braggfix.sh -- the seventeen fake-data ensembles of the release (100 members each), re-run on the
# selection without the Bragg-pion requirement (version 1.7, 2026-10-03). Inclusive: FHC p_mu, FHC and RHC
# cos theta_mu, combined p_pi (three regions), the FHC and RHC one-bin totals and the combined
# cos theta_pi x cos theta_mu. Proton-tagged: FHC W_had and delta p_T, combined p_n, theta_p and theta_pi_p in
# FHC, RHC and combined, and the combined theta_p x delta p_T. Same member jobs as ens_050.sh and ens_v1.sh;
# the bin, slice and xsec configs are those of the release (bin configs checked byte-identical to the
# configs/braggfix/ set the release universes were built from).
#
# The previous outputs and throws are moved to ens/pre_braggfix/ and w/ens/pre_braggfix/ first. Every existing
# throw is older than the reprocessed trees, so throw_guard.h would replace it by rename while other members
# read it (the NFS race of 2026-09-26); with empty throw directories every throw is created once by link(2).
# Two CPUs and 5 GB per member (peak resident memory of the previous 1300 members: 2.8 GB), so eight members
# fit on a 16-core, 60 GB node. Combined arrays are submitted first: their members take about twice as long.
set -uo pipefail
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO/slurm"
PROC=/data/uboone/processed; E=$PROC/ens; A=$E/pre_braggfix; AW=$PROC/w/ens/pre_braggfix
squeue -h -u "$USER" -o '%j' | grep -q -E '^ens' && { echo "ensemble jobs are still queued or running"; exit 1; }
mkdir -p "$A" "$AW"
if [ ! -e "$A/.archived_throws" ]; then
  mv -n "$E/throws" "$A/throws" && mkdir "$E/throws" || { echo "FAILED to archive $E/throws"; exit 1; }
  mv -n "$PROC/w/ens/throws" "$AW/throws" && mkdir "$PROC/w/ens/throws" || { echo "FAILED to archive w/ens/throws"; exit 1; }
  touch "$A/.archived_throws"; echo "archived the throws"
fi
INCL="comb:ppi3bin:configs comb:costhpi_costhmu_2d:configs/d2 fhc5:pmu:configs fhc5:costhmu:configs rhcfull:costhmu:configs fhc5:total:configs rhcfull:total:configs"
P1="comb:pn2bin:configs comb:thetap:configs comb:thpipr:configs comb:thetap_dpt_2d:configs/d2 fhc5:Whad:configs fhc5:dpt2bin:configs fhc5:thetap:configs fhc5:thpipr:configs rhcfull:thetap:configs rhcfull:thpipr:configs"
# once per tag (marker file): a second pass would move this run's files onto the archived ones
for s in $INCL; do IFS=: read -r c o _ <<< "$s"; TAGS="${TAGS:-} ${c}_${o}"; done
for s in $P1;   do IFS=: read -r c o _ <<< "$s"; TAGS="$TAGS 1p_${c}_${o}"; done
for tag in $TAGS; do
  [ -e "$A/.archived_$tag" ] && { echo "already archived: $tag"; continue; }
  n=0; for f in "$E"/*_${tag}_t[0-9]*; do [ -e "$f" ] || continue; mv -n "$f" "$A/"; n=$((n+1)); done
  touch "$A/.archived_$tag"; echo "archived $n files of $tag"
done
LOGF=rebuild_jobs_2026-09-19.txt
RES="--cpus-per-task=2 --mem=5G --time=16:00:00"
sub() {  # family config observable conf throttle
  local script=slurm_ensemble_cfg.sbatch lab="incl"; [ "$1" = 1p ] && { script=slurm_ensemble_1p.sbatch; lab="1p"; }
  local J; J=$(sbatch --parsable $RES --array=0-99%$5 --export=ALL,REPO=$REPO,PROC=$PROC,CONF=$4,UNIV_THREADS=2,OBS=$3,CFG=$2 $script) \
    || { echo "FAILED to submit $lab $2 $3"; return 1; }
  echo "$J ensemble $lab $2 $3 (100 members, release without the Bragg-pion cut) $(date)" >> $LOGF; echo "submitted $J $lab $2 $3"
}
for pass in comb other; do
  for s in $INCL; do IFS=: read -r c o conf <<< "$s"; [ "$c" = comb ] && [ $pass = comb ] && sub incl "$c" "$o" "$conf" 12; [ "$c" != comb ] && [ $pass = other ] && sub incl "$c" "$o" "$conf" 7; done
  for s in $P1;   do IFS=: read -r c o conf <<< "$s"; [ "$c" = comb ] && [ $pass = comb ] && sub 1p "$c" "$o" "$conf" 12; [ "$c" != comb ] && [ $pass = other ] && sub 1p "$c" "$o" "$conf" 7; done
done
