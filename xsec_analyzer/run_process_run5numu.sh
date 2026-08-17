#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
in=/data/uboone/new_numi_flux/reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root
out=/data/uboone/processed/xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root
echo "[proc] Run5 numu overlay"; ProcessNTuples "$in" numuMC CC1mu1piXp "$out" && echo "Run5numu OK" || echo "Run5numu FAILED"
echo "########## PROCESS RUN5NUMU DONE ##########"
ls -la "$out" 2>/dev/null | awk '{print $5/1e9" GB"}'
