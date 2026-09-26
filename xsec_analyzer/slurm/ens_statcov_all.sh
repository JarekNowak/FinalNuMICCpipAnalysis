#!/usr/bin/env bash
# ens_statcov_all.sh -- run ens_statcov.sh for the ensembles of the release in parallel (incremental; members
# already done are skipped). 2026-09-26: the four rebinned inclusive ensembles, the two one-bin totals of
# 2026-09-20 and the proton-tagged W_had (two regions).
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
L=../logs/ens; mkdir -p $L
for pair in "fhc5 pmu" "fhc5 costhmu" "rhcfull costhmu" "comb ppi3bin" "fhc5 total" "rhcfull total" "1p_fhc5 Whad"; do
  set -- $pair; nice -n 10 bash slurm/ens_statcov.sh $1 $2 100 > $L/statcov_$1_$2.log 2>&1 &
done
wait; echo "ALL_STATCOV_DONE $(date +%H:%M)"; tail -n 1 $L/statcov_*.log
