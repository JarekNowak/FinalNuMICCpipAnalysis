#!/usr/bin/env bash
# consolidated_wtki.sh — one autonomous pass to fix the proton-tagged W/TKI results:
#   - tightened proton PID (LLR<0.05, baked into the recompiled CC1mu1pi1p selection)
#   - coarser 3-bin observables (stabilises the low-statistics RHC/comb Wiener-SVD,
#     which oscillated negative at 4 bins)
# Steps: reprocess the w/ files (FHC + RHC) with the new selection, re-throw fake data,
# then univmake + unfold FHC / RHC / combined with the 3-bin ccpi1p configs.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"
echo "==== CONSOLIDATED W/TKI (LLR<0.05, 3-bin) START $(date) ===="

# clear the old (LLR<0.2 / 4-bin) univmakes + closure sidecars so everything regenerates
rm -f /data/uboone/processed/ccpi1p_*_univmake.root
rm -f /data/uboone/processed/closure_hists_xsec_ccpi1p_*.root
rm -f /data/uboone/processed/xsec_ccpi1p_*.root

# 1. reprocess w/ with the tightened selection (FHC first: it also does the shared
#    dirt + combined beamoff that RHC symlinks to; then RHC).
echo "---- [1] reprocess w/ FHC (LLR<0.05) ----"
bash reprocess_w_fhc.sh
echo "---- [1] reprocess w/ RHC (LLR<0.05) ----"
bash reprocess_w_rhc.sh

# 2. univmake + unfold with the 3-bin configs
echo "---- [2] FHC univmake+unfold (3-bin) ----"
bash w_batch_fhc.sh 2
echo "---- [2] RHC+comb univmake+unfold (3-bin) ----"
bash w_batch_rhc_comb.sh 2

echo "==== CONSOLIDATED W/TKI DONE $(date) ===="
echo "closure sigma_int:"
grep -h "OK " ../logs/w_*_unfold.log 2>/dev/null | tail -20
