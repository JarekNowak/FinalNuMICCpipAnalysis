#!/usr/bin/env bash
# reprocess_sidebands.sh — reprocess the sideband-relevant files with the rebuilt
# ProcessNTuples so every stv_tree carries the four background-control sideband flags
# (CC1mu1piXp_sb_cc0pi / sb_multipi / sb_pi0 / sb_cosmic). Writes to a SEPARATE
# directory (/data/uboone/processed/sb) so it NEVER overwrites the standard
# xsec-ana-*.root files the running comb univmake batch is still reading. Runs nice'd
# at low parallelism so the batch keeps CPU priority. After the MC lands, re-throws the
# per-run fake data from the sb/ MC (CloneTree carries sb_ through) so the fake "data"
# also has the sideband flags. EXT + dirt reprocessed too.  usage: ./reprocess_sidebands.sh
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh 2>/dev/null; set -u
OUT=/data/uboone/processed/sb
RAW=/data/uboone/new_numi_flux
LOG=../logs/reprocess_sb; mkdir -p "$LOG"
STATUS=$LOG/status.txt; : > "$STATUS"
MAXJOBS=3
NICE="nice -n 15"

proc_one() {  # rawpath  outname  filetype
  local in="$1" out="$OUT/xsec-ana-$2.root" ftype="$3"
  [[ -f "$in" ]] || { echo "FAIL $2 (no raw input: $in)" >> "$STATUS"; return; }
  $NICE ProcessNTuples "$in" "$ftype" CC1mu1piXp "$out" > "$LOG/$2.log" 2>&1
  local rc=$?
  local ok=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");printf(\"%d\n\",(f&&!f->IsZombie()&&((TTree*)f->Get(\"stv_tree\"))&&((TTree*)f->Get(\"stv_tree\"))->GetBranch(\"CC1mu1piXp_sb_cosmic\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$ok" == "1" ]]; then echo "OK   $2" >> "$STATUS"; else echo "FAIL $2 rc=$rc sb=$ok" >> "$STATUS"; fi
}
export -f proc_one; export OUT LOG STATUS NICE

echo "==== SIDEBAND REPROCESS START $(date) ====" | tee -a "$STATUS"

# --- MC (14): FHC 4 + RHC 10 ---
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
# --- EXT + dirt (share the pool) ---
proc_one /data/uboone/EXT/beamoff_run1Andrun3.root beamoff_run1Andrun3 extBNB &
while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
proc_one /data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot dirtMC &
wait
echo "---- MC+EXT+dirt reprocess done $(date) ----" | tee -a "$STATUS"

# --- re-throw per-run fake data FROM the sb/ MC (sb_ flags carried through CloneTree) ---
$NICE root.exe -l -b -q "macros/throw_perrun_sb.C(1)" > "$LOG/throw_sb.log" 2>&1 \
  && echo "OK   fakedata re-throw (sb)" >> "$STATUS" || echo "FAIL fakedata re-throw" >> "$STATUS"

echo "==== SIDEBAND REPROCESS DONE $(date) ====" | tee -a "$STATUS"
sort "$STATUS"
