#!/usr/bin/env bash
cd "$(dirname "${BASH_SOURCE[0]}")"; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
L=../logs/sidebands_perrun; mkdir -p $L
for m in fhc rhc comb; do
  root.exe -l -b -q "macros/sideband_compare.C(\"$m\",\"fake\",\"/data/uboone/processed/sb/\",\"CC1mu1piXp\",false)" > $L/incl_$m.log 2>&1
  root.exe -l -b -q "macros/sideband_compare.C(\"$m\",\"fake\",\"/data/uboone/processed/w/\",\"CC1mu1pi1p\",false)" > $L/1p_notag_$m.log 2>&1
  root.exe -l -b -q "macros/sideband_compare.C(\"$m\",\"fake\",\"/data/uboone/processed/w/\",\"CC1mu1pi1p\",true)" > $L/1p_tag_$m.log 2>&1
  echo "$m done $(date +%H:%M)"
done; echo SIDEBANDS_DONE
