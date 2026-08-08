#!/usr/bin/env bash
# Reprocess the 14 numuMC overlay run files with the rebuilt ProcessNTuples so each
# processed file gains the per-cut RECO cut-flow histograms (h_cf_*). Overwrites the
# existing xsec-ana-*.root (resets summed_pot to native — fine, combined build is done).
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:$LD_LIBRARY_PATH"
source ./setup_xsec_analyzer.sh 2>/dev/null
RAW=/data/uboone/new_numi_flux
OUT=/data/uboone/processed
LOG=../logs/reprocess_cf; mkdir -p "$LOG"
STATUS=$LOG/status.txt; : > "$STATUS"
MAXJOBS=5

names=(
  Run1_fhc_new_numi_flux_fhc_pandora_ntuple
  Run2_fhc_new_numi_flux_fhc_pandora_ntuple
  Run4_fhc_new_numi_flux_fhc_pandora_ntuple
  reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc
  Run1_rhc_new_numi_flux_rhc_pandora_ntuple
  Run2_rhc_new_numi_flux_rhc_pandora_ntuple
  Run4a_rhc_new_numi_flux_rhc_pandora_ntuple
  Run4b_rhc_new_numi_flux_rhc_pandora_ntuple
  Run4c_rhc_new_numi_flux_rhc_pandora_ntuple
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae
)

proc_one() {
  local n="$1"
  local in="$RAW/$n.root" out="$OUT/xsec-ana-$n.root"
  [[ -f "$in" ]] || { echo "FAIL $n (no raw input)" >> "$STATUS"; return; }
  ProcessNTuples "$in" numuMC CC1mu1piXp "$out" > "$LOG/$n.log" 2>&1
  local rc=$?
  # verify h_cf histograms present
  local ok=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");printf(\"%d\n\",(f&&!f->IsZombie()&&f->Get(\"h_cf_pmu_cut9_final\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$ok" == "1" ]]; then echo "OK   $n" >> "$STATUS"; else echo "FAIL $n rc=$rc h_cf=$ok" >> "$STATUS"; fi
}
export -f proc_one; export RAW OUT LOG STATUS

for n in "${names[@]}"; do
  proc_one "$n" &
  while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
done
wait
echo "########## REPROCESS CUTFLOW DONE ##########" >> "$STATUS"
sort "$STATUS"
