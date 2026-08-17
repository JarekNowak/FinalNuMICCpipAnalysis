#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
in=/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root
out=/data/uboone/processed/xsec-ana-Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root
echo "[proc] Run5_fhc"; ProcessNTuples "$in" numuMC CC1mu1piXp "$out" && echo "Run5 OK" || echo "Run5 FAILED"
echo "########## PROCESS RUN5 DONE ##########"
ls -la "$out" 2>/dev/null | awk '{print $5/1e9" GB"}'
