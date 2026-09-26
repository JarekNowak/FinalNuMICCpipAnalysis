#!/usr/bin/env bash
# run_mcs_050.sh -- the data-side MCS muon-momentum scale (+-5% on the MCS-measured muons of the fake data,
# nominal response) for the six-bin p_mu adopted under the 0.50 migration criterion (2026-09-26), in FHC,
# RHC and combined. The other inclusive observables carry no reco-p_mu cut and stay exactly zero.
# Configs derive from the staged release configs (configs/adopt050); replaced outputs go to
# rebuild_alt/pre050_20260926/ and the old configs to configs/mcs/pre050_20260926/.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
RA=/data/uboone/processed/rebuild_alt; BK=$RA/pre050_20260926; CB=configs/mcs/pre050_20260926; mkdir -p "$BK" "$CB"
M=slurm/mcsdata_050_manifest.list; : > "$M"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB ); declare -A SFX=( [fhc5]="" [rhcfull]="_rhcfull" [comb]="_comb" )
for c in fhc5 rhcfull comb; do T=${TAG[$c]}
  for v in up05 dn05; do
    for f in "$RA"/ccpi_${T}_pmu_mcsdata${v}_univmake.root "$RA"/xsec_ccpi_${T}_pmu_mcsdata${v}.root "$RA"/closure_hists_xsec_ccpi_${T}_pmu_mcsdata${v}.root; do
      [ -f "$f" ] && mv "$f" "$BK/"; done
    out=configs/mcs/ccpi_xsec_config_numi_pmu${SFX[$c]}_mcsdata_${v}.txt
    [ -f "$out" ] && [ ! -f "$CB/$(basename "$out")" ] && cp -p "$out" "$CB/"
    sed -E "s#^UnivFile .*#UnivFile $RA/ccpi_${T}_pmu_mcsdata${v}_univmake.root#; s#^FPFile .*#FPFile configs/mcs/file_properties_numi_${c}_mcsdata_${v}.txt#" \
      configs/adopt050/ccpi_xsec_config_numi_pmu_${c}.txt > "$out"
    echo "$out configs/adopt050/ccpi_pmu_bin_config_opt.txt configs/adopt050/ccpi_pmu_slice_config_opt.txt" >> "$M"
  done
done
J=$(sbatch --parsable --cpus-per-task=4 --mem=12G --array=0-5 --export=ALL,REPO=$PWD,MANIFEST=$PWD/$M,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
echo "$J $(basename $M) (data-side MCS scale +-5% for the six-bin p_mu of the 0.50 criterion, FHC/RHC/COMB) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J"
