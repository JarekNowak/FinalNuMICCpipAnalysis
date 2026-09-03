#!/usr/bin/env bash
# post_rethrow.sh -- after rethrow_fakedata.sh: (1) rerun the sideband macros on the re-thrown
# sb/ pseudo-data; (2) submit the 51-extraction universe rebuild (re-thrown fake data are
# univmake inputs) into a clean $PROC/rebuild_extfix for later promotion by promote_extfix.sh.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
L=../logs/sidebands_perrun
./regen_sidebands.sh > ../logs/sidebands_perrun_driver3.log 2>&1
for m in fhc rhc; do root -l -b -q "macros/sb_overlap.C(\"$m\")" > $L/overlap3_$m.log 2>&1; root -l -b -q "macros/sb_plots.C(\"$m\")" > $L/plots3_$m.log 2>&1; done
for m in fhc rhc comb; do root.exe -l -b -q "macros/sideband_compare.C(\"$m\",\"fakestack\",\"/data/uboone/processed/sb/\",\"CC1mu1piXp\",false)" > $L/fakestack3_$m.log 2>&1; done
echo "SIDEBANDS3_DONE $(date)"
PROC=/data/uboone/processed; rm -rf $PROC/rebuild_extfix; mkdir -p $PROC/rebuild_extfix
cd slurm && sbatch --array=0-50 --partition=main --export="ALL,REPO=$(cd .. && pwd),PROC=$PROC,MANIFEST=$(pwd)/rebuild_manifest.list" --chdir=$(pwd) slurm_extfix.sbatch
echo "SUBMITTED $(date)"
