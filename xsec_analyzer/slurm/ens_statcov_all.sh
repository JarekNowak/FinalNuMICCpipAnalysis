#!/usr/bin/env bash
# ens_statcov_all.sh -- run ens_statcov.sh for the six 2026-09-19 ensembles in parallel (incremental).
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
L=../logs/ens; mkdir -p $L
for pair in "fhc5 pmu" "fhc5 costhmu" "rhcfull costhmu" "comb ppi2bin" "fhc5 total" "rhcfull total"; do
  set -- $pair; nice -n 10 bash slurm/ens_statcov.sh $1 $2 100 > $L/statcov_$1_$2.log 2>&1 &
done
wait; echo "ALL_STATCOV_DONE $(date +%H:%M)"; tail -n 1 $L/statcov_*.log
