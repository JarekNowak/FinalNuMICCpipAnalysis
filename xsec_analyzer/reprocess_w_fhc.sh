#!/usr/bin/env bash
# reprocess_w_fhc.sh — reprocess the FHC files with the CC1mu1pi1p (proton-tagged
# W/TKI) selection into a separate w/ dir, named to match file_properties so a
# w-variant properties file only swaps the directory. 4 numuMC + dirt + combined
# beamoff (symlinked to the 4 per-run EXT names), then re-throw the 4 per-run fake
# data from the w/ MC. detVar deferred to a follow-up (dominant flux/xsec/reint
# systematics come from the MC weight branches, present here).  nice'd, batch-safe.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh 2>/dev/null; set -u
OUT=/data/uboone/processed/w
RAW=/data/uboone/new_numi_flux
LOG=../logs/reprocess_w; mkdir -p "$LOG" "$OUT"
STATUS=$LOG/status.txt; : > "$STATUS"
MAXJOBS=3; NICE="nice -n 12"; SEL=CC1mu1pi1p

proc_one() {  # rawpath outname filetype
  local in="$1" out="$OUT/xsec-ana-$2.root" ftype="$3"
  [[ -f "$in" ]] || { echo "FAIL $2 (no raw: $in)" >> "$STATUS"; return; }
  $NICE ProcessNTuples "$in" "$ftype" $SEL "$out" > "$LOG/$2.log" 2>&1
  local rc=$?
  local ok=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");printf(\"%d\n\",(f&&!f->IsZombie()&&((TTree*)f->Get(\"stv_tree\"))&&((TTree*)f->Get(\"stv_tree\"))->GetBranch(\"${SEL}_W_pipr_reco\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$ok" == "1" ]]; then echo "OK   $2" >> "$STATUS"; else echo "FAIL $2 rc=$rc W=$ok" >> "$STATUS"; fi
}
export -f proc_one; export OUT LOG STATUS NICE SEL

echo "==== W/TKI FHC REPROCESS START $(date) ====" | tee -a "$STATUS"
mc=( Run1_fhc_new_numi_flux_fhc_pandora_ntuple Run2_fhc_new_numi_flux_fhc_pandora_ntuple
     Run4_fhc_new_numi_flux_fhc_pandora_ntuple reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc )
for n in "${mc[@]}"; do
  proc_one "$RAW/$n.root" "$n" numuMC &
  while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
done
proc_one /data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot dirtMC &
while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
proc_one /data/uboone/EXT/beamoff_run1Andrun3.root beamoff_run1Andrun3 extBNB &
wait
# per-run EXT symlinks (same combined beamoff, matching file_properties names)
for r in run1 run2 run4 run5; do ln -sf "$OUT/xsec-ana-beamoff_run1Andrun3.root" "$OUT/xsec-ana-beamoff_fhc_${r}.root"; done
echo "---- MC+dirt+EXT reprocess done $(date) ----" | tee -a "$STATUS"

$NICE root.exe -l -b -q "macros/throw_perrun_w.C(1)" > "$LOG/throw_w.log" 2>&1 \
  && echo "OK   fakedata re-throw (w)" >> "$STATUS" || echo "FAIL fakedata re-throw" >> "$STATUS"
echo "==== W/TKI FHC REPROCESS DONE $(date) ====" | tee -a "$STATUS"
sort "$STATUS"
