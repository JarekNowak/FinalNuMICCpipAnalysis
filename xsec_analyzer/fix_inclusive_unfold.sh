#!/usr/bin/env bash
# fix_inclusive_unfold.sh — re-unfold the INCLUSIVE (CC1mu1piXp) RHC + combined
# differential cross sections after the per-run POT normalization fix. FHC is NOT
# re-run (one MC file per run -> immune). Uses XSEC_FORCE_REBUILD so build_universes()
# recomputes the POT-summed universes with the corrected multi-file-per-run scaling
# instead of loading the stale (buggy) cache. Covers all 5 observables in both the
# Wiener-SVD (default) and D'Agostini (_dag) variants. Serial to limit I/O contention.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
export XSEC_FORCE_REBUILD=1
PROC=/data/uboone/processed
declare -A TAG=( [rhcfull]=RHCFULL [comb]=COMB )
echo "==== FIX INCLUSIVE RHC+COMB UNFOLD (force-rebuild) START $(date) ===="
for cfg in rhcfull comb; do
  t=${TAG[$cfg]}
  for o in pmu ppi costhmu costhpi thmupi; do
    uni=$PROC/ccpi_${t}_${o}_univmake.root
    [ -s "$uni" ] || { echo "  !! $cfg $o univmake missing"; continue; }
    # Wiener-SVD
    xc=configs/ccpi_xsec_config_numi_${o}_${cfg}.txt
    if [ -f "$xc" ]; then
      nice bin/UnfolderNuMI "$xc" configs/ccpi_${o}_slice_config_opt.txt \
        $PROC/xsec_${t}_${o}.root > logs/incfix_${cfg}_${o}.log 2>&1 \
        && echo "  OK WSVD $t $o  sigma_int=$(grep 'SYSTDUMP] sigma_int' logs/incfix_${cfg}_${o}.log | tail -1 | awk '{print $3}')" \
        || echo "  !! WSVD $t $o FAILED"
    fi
    # D'Agostini cross-check (if a _dag config exists)
    xcd=configs/ccpi_xsec_config_numi_${o}_${cfg}_dag.txt
    if [ -f "$xcd" ]; then
      nice bin/UnfolderNuMI "$xcd" configs/ccpi_${o}_slice_config_opt.txt \
        $PROC/xsec_${t}_${o}_dag.root > logs/incfix_${cfg}_${o}_dag.log 2>&1 \
        && echo "  OK DAG  $t $o  sigma_int=$(grep 'SYSTDUMP] sigma_int' logs/incfix_${cfg}_${o}_dag.log | tail -1 | awk '{print $3}')" \
        || echo "  !! DAG  $t $o FAILED"
    fi
    echo "DONE_${cfg}_${o}" >> logs/fix_inclusive.status
  done
done
echo "ALL_DONE" >> logs/fix_inclusive.status
echo "==== FIX INCLUSIVE DONE $(date) ===="
