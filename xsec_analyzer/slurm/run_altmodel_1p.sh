#!/usr/bin/env bash
# run_altmodel_1p.sh -- throw the proton-tagged (w/) alternative-model fake data (same universes as
# the inclusive test: GENIE multisim u545 = altgenie, Delta->N pi angular unisim = altdelta), then
# submit the 12 FHC extractions of slurm/altmodel_1p_manifest.list (4 CPU / 8 GB, backfill).
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$PWD"
L=../logs/altmodel; mkdir -p $L
U=${ALT_UNIV:-545}
echo "throws start $(date +%H:%M) (GENIE universe $U)"
nice -n 10 root -l -b -q "macros/throw_altmodel_w.C(\"weight_All_UBGenie\",$U,\"altgenie\",1)" > $L/throw_w_altgenie.log 2>&1 || { echo "altgenie throw FAILED"; exit 1; }
nice -n 10 root -l -b -q 'macros/throw_altmodel_w.C("weight_Theta_Delta2Npi_UBGenie",0,"altdelta",1)' > $L/throw_w_altdelta.log 2>&1 || { echo "altdelta throw FAILED"; exit 1; }
grep "POTSCALE" $L/throw_w_altgenie.log $L/throw_w_altdelta.log
for f in $(grep -hv "^#" configs/file_properties_numi_fhc5_w_altgenie.txt configs/file_properties_numi_fhc5_w_altdelta.txt | awk '{print $1}' | grep fakedata); do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done
J=$(sbatch --parsable --cpus-per-task=4 --mem=8G --array=0-11 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/altmodel_1p_manifest.list slurm/slurm_extract.sbatch)
echo "$J altmodel_1p_manifest.list (proton-tagged alt-model closure, FHC, 6 obs x altgenie/altdelta) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J $(date +%H:%M)"
