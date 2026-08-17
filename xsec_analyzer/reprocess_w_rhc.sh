#!/usr/bin/env bash
# reprocess_w_rhc.sh — reprocess the 10 RHC numuMC files with CC1mu1pi1p into w/,
# symlink the 4 per-run RHC EXT to the (already reprocessed) combined beamoff, and
# re-throw the 4 per-run RHC fake data. dirt + beamoff already reprocessed by the FHC
# step. This completes the inputs for BOTH the RHC and the combined W/TKI configs
# (comb reuses the FHC+RHC w/ files). detVar deferred.  nice'd, batch-safe.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh 2>/dev/null; set -u
OUT=/data/uboone/processed/w
RAW=/data/uboone/new_numi_flux
LOG=../logs/reprocess_w_rhc; mkdir -p "$LOG" "$OUT"
STATUS=$LOG/status.txt; : > "$STATUS"
MAXJOBS=3; NICE="nice -n 12"; SEL=CC1mu1pi1p

proc_one() {  # rawname
  local in="$RAW/$1.root" out="$OUT/xsec-ana-$1.root"
  [[ -f "$in" ]] || { echo "FAIL $1 (no raw: $in)" >> "$STATUS"; return; }
  $NICE ProcessNTuples "$in" numuMC $SEL "$out" > "$LOG/$1.log" 2>&1
  local rc=$?
  local ok=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");printf(\"%d\n\",(f&&!f->IsZombie()&&((TTree*)f->Get(\"stv_tree\"))&&((TTree*)f->Get(\"stv_tree\"))->GetBranch(\"${SEL}_W_pipr_reco\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$ok" == "1" ]]; then echo "OK   $1" >> "$STATUS"; else echo "FAIL $1 rc=$rc W=$ok" >> "$STATUS"; fi
}
export -f proc_one; export OUT RAW LOG STATUS NICE SEL

echo "==== W/TKI RHC REPROCESS START $(date) ====" | tee -a "$STATUS"
mc=( Run1_rhc_new_numi_flux_rhc_pandora_ntuple Run2_rhc_new_numi_flux_rhc_pandora_ntuple
     Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab
     Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad
     Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae
     Run4a_rhc_new_numi_flux_rhc_pandora_ntuple Run4b_rhc_new_numi_flux_rhc_pandora_ntuple
     Run4c_rhc_new_numi_flux_rhc_pandora_ntuple )
for n in "${mc[@]}"; do
  proc_one "$n" &
  while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
done
wait
# per-run RHC EXT symlinks (same combined beamoff already in w/)
for r in run1 run2 run3 run4; do ln -sf "$OUT/xsec-ana-beamoff_run1Andrun3.root" "$OUT/xsec-ana-beamoff_rhc_${r}.root"; done
echo "---- RHC MC reprocess done $(date) ----" | tee -a "$STATUS"

$NICE root.exe -l -b -q -e 'gROOT->LoadMacro("macros/throw_perrun_w.C"); throw_perrun_w_rhc(1);' > "$LOG/throw_rhc.log" 2>&1 \
  && echo "OK   fakedata RHC re-throw (w)" >> "$STATUS" || echo "FAIL fakedata RHC re-throw" >> "$STATUS"
echo "==== W/TKI RHC REPROCESS DONE $(date) ====" | tee -a "$STATUS"
sort "$STATUS"
