#!/usr/bin/env bash
# Reprocess numuMC(14) + EXT(beamoff_run1Andrun3) + dirt(1) with the rebuilt
# instrumented ProcessNTuples so every processed file carries the diagnostic
# histograms: h_cutflow_tot/sig, h_cf_*, h_nm1_topo/oa_{sig,bkg}, h_fin_*_{sig,bkg},
# h_bkgcat. Binning-independent (stv_tree events unchanged -> univmakes stay valid).
# Does NOT touch the fake-data throw or detVar files.
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh 2>/dev/null; set -u
OUT=/data/uboone/processed
LOG=../logs/reprocess_diag; mkdir -p "$LOG"
STATUS=$LOG/status.txt; : > "$STATUS"
MAXJOBS=5

proc_one() {  # rawpath  outname  filetype
  local in="$1" out="$OUT/xsec-ana-$2.root" ftype="$3"
  [[ -f "$in" ]] || { echo "FAIL $2 (no raw input: $in)" >> "$STATUS"; return; }
  ProcessNTuples "$in" "$ftype" CC1mu1piXp "$out" > "$LOG/$2.log" 2>&1
  local rc=$?
  local ok=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");printf(\"%d\n\",(f&&!f->IsZombie()&&f->Get(\"h_cutflow_tot\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$ok" == "1" ]]; then echo "OK   $2" >> "$STATUS"; else echo "FAIL $2 rc=$rc cutflow=$ok" >> "$STATUS"; fi
}
export -f proc_one; export OUT LOG STATUS

RAW=/data/uboone/new_numi_flux
numu=(
  Run1_fhc_new_numi_flux_fhc_pandora_ntuple Run2_fhc_new_numi_flux_fhc_pandora_ntuple
  Run4_fhc_new_numi_flux_fhc_pandora_ntuple reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc
  Run1_rhc_new_numi_flux_rhc_pandora_ntuple Run2_rhc_new_numi_flux_rhc_pandora_ntuple
  Run4a_rhc_new_numi_flux_rhc_pandora_ntuple Run4b_rhc_new_numi_flux_rhc_pandora_ntuple
  Run4c_rhc_new_numi_flux_rhc_pandora_ntuple
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad
  Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae
)
for n in "${numu[@]}"; do
  proc_one "$RAW/$n.root" "$n" numuMC &
  while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
done
# EXT + dirt
proc_one /data/uboone/EXT/beamoff_run1Andrun3.root beamoff_run1Andrun3 extBNB &
proc_one /data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot dirtMC &
wait
echo "########## REPROCESS DIAGNOSTICS DONE ##########" >> "$STATUS"
sort "$STATUS"
