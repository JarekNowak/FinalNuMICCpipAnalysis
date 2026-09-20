#!/usr/bin/env bash
# run_altmodel_fhc_rethrow.sh -- the inclusive FHC alternative-model fake data were thrown on 2026-09-18
# at the old Run-1 exposure (summed_pot 3.283e20) and the closure built self-consistently on it; the
# file lists now carry 2.192e20, so the files are re-thrown (same universes 545 / Delta unisim, same
# seeds) and the 8 FHC extractions of slurm/altmodel_manifest.list re-run at 4 CPU / 8 GB.
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$PWD"
L=../logs/altmodel; mkdir -p $L
for r in 1 2 4 5; do for t in altgenie altdelta; do f=/data/uboone/processed/xsec-ana-fakedata_${t}_fhc_run$r.root; [ -f $f ] && mv $f ${f%.root}_pot3283.root; done; done
echo "throws start $(date +%H:%M)"
nice -n 10 root -l -b -q 'macros/throw_altmodel_fhc.C("weight_All_UBGenie",545,"altgenie",1)' > $L/throw_fhc_altgenie_rethrow.log 2>&1 || { echo "altgenie throw FAILED"; exit 1; }
nice -n 10 root -l -b -q 'macros/throw_altmodel_fhc.C("weight_Theta_Delta2Npi_UBGenie",0,"altdelta",1)' > $L/throw_fhc_altdelta_rethrow.log 2>&1 || { echo "altdelta throw FAILED"; exit 1; }
grep -h "POTSCALE" $L/throw_fhc_altgenie_rethrow.log $L/throw_fhc_altdelta_rethrow.log
J=$(sbatch --parsable --cpus-per-task=4 --mem=8G --array=0-7 --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/altmodel_manifest.list slurm/slurm_extract.sbatch)
echo "$J altmodel_manifest.list (inclusive FHC alt-model closure re-thrown at 2.192e20) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J $(date +%H:%M)"
