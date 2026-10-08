#!/bin/bash
# SLURM array task: univmake of the two-pion one-bin total with its sideband reco bins (Phase 4, item 3, of
# report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md) for ONE configuration (fhc5, rhcfull, comb), from the production in
# /data/uboone/processed/mp2pi (scripts/mp2pi_setup.py writes the file lists). univmake is multi-threaded and uses
# ~13 of 16 cores, so each task takes a whole node.
#   sbatch --array=0-2 --export=ALL,REPO=$PWD slurm/slurm_univmake_mp2pi.sbatch
#SBATCH --job-name=univmp2pi
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.out
set -u
REPO="${REPO:?}"
CFGS=(fhc5 rhcfull comb); declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
CFG=${CFGS[$SLURM_ARRAY_TASK_ID]}
cd "$REPO" || exit 1
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$REPO"
OUT="/data/uboone/processed/mp2pi_${TAG[$CFG]}_total_univmake.root"
FPM="configs/file_properties_numi_${CFG}_mp2pi.txt"; BIN="configs/mp2pi_total_bin_config.txt"
for f in "$FPM" "$BIN"; do [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
echo "host=$(hostname) cfg=$CFG FPM=$FPM OUT=$OUT"
TMP="${OUT}.partial.$$"; t0=$(date +%s)
"$REPO/bin/univmake" "$FPM" "$BIN" "$TMP" "$FPM"; rc=$?; t1=$(date +%s)
if [ $rc -ne 0 ] || [ ! -s "$TMP" ]; then echo "FAILED $CFG rc=$rc after $(( (t1-t0)/60 ))min" >&2; rm -f "$TMP"; exit 1; fi
mv -f "$TMP" "$OUT"
echo "OK $CFG $(( $(stat -c%s "$OUT")/1000000 ))MB in $(( (t1-t0)/60 ))min"
