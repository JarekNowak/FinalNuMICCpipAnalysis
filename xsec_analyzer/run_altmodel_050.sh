#!/usr/bin/env bash
# run_altmodel_050.sh -- alternative-model closure (GENIE multisim universe 545 and the Delta->N pi angular
# unisim, FHC, fake data already thrown at the corrected exposure) for the binnings adopted under the 0.50
# migration criterion on 2026-09-26. Re-run: inclusive p_mu, cos theta_mu and p_pi (three regions); proton-
# tagged W_had (two regions) and delta phi_T (three bins); with PP_K set, also the new proton momentum p_p.
# The unchanged binnings (cos theta_pi, the two-bin delta p_T, delta alpha_T, p_n) keep their results.
# With RHC=1, the RHC extractions of the inclusive p_mu, cos theta_mu and p_pi instead.
# Configs are derived from the staged release configs (configs/adopt050), so this can run before promotion;
# the outputs these replace are moved to rebuild_alt/pre050_20260926/ first.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
RA=/data/uboone/processed/rebuild_alt; BK=$RA/pre050_20260926; CB=configs/alt/pre050_20260926
mkdir -p "$BK" "$CB"
M=slurm/altmodel_050$([ "${RHC:-0}" = 1 ] && echo _rhc)_manifest.list; : > "$M"   # RHC mode has its own list: running jobs read theirs
add() {   # pfx obs bin slice fpsuffix [cfg: fhc5 | rhcfull]
  local pfx=$1 obs=$2 bin=$3 slice=$4 fps=$5 cfg=${6:-fhc5}
  local T=FHC5 cs=""; [ "$cfg" = rhcfull ] && { T=RHCFULL; cs="_rhc"; }
  for v in altgenie altdelta; do
    for f in "$RA"/${pfx}_${T}_${obs}_${v}_univmake.root "$RA"/xsec_${pfx}_${T}_${obs}_${v}.root "$RA"/closure_hists_xsec_${pfx}_${T}_${obs}_${v}.root "$RA"/xsec_${pfx}_${T}_${obs}_${v}.log; do
      [ -f "$f" ] && mv "$f" "$BK/"; done
    local out=configs/alt/${pfx}_xsec_config_numi_${obs}${cs}_${v}.txt
    [ -f "$out" ] && [ ! -f "$CB/$(basename "$out")" ] && cp -p "$out" "$CB/"
    sed -E "s#^UnivFile .*#UnivFile $RA/${pfx}_${T}_${obs}_${v}_univmake.root#; s#^FPFile .*#FPFile configs/file_properties_numi_${cfg}${fps}_${v}.txt#" \
      configs/adopt050/${pfx}_xsec_config_numi_${obs}_${cfg}.txt > "$out"
    echo "$out configs/adopt050/$bin configs/adopt050/$slice" >> "$M"
  done
}
# ONLY_PP=1: just the proton momentum (run after PP_K is chosen and staged), leaving the others untouched
# RHC=1: the RHC rows of the inclusive table (p_mu, cos theta_mu, p_pi), leaving the FHC results untouched
if [ "${RHC:-0}" = 1 ]; then
  add ccpi pmu ccpi_pmu_bin_config_opt.txt ccpi_pmu_slice_config_opt.txt "" rhcfull
  add ccpi costhmu ccpi_costhmu_bin_config_opt.txt ccpi_costhmu_slice_config_opt.txt "" rhcfull
  add ccpi ppi3bin ccpi_ppi_bin_config_3bin.txt ccpi_ppi_slice_config_3bin.txt "" rhcfull
elif [ "${ONLY_PP:-0}" != 1 ]; then
  add ccpi pmu ccpi_pmu_bin_config_opt.txt ccpi_pmu_slice_config_opt.txt ""
  add ccpi costhmu ccpi_costhmu_bin_config_opt.txt ccpi_costhmu_slice_config_opt.txt ""
  add ccpi ppi3bin ccpi_ppi_bin_config_3bin.txt ccpi_ppi_slice_config_3bin.txt ""
  add ccpi1p Whad ccpi1p_Whad_bin_config.txt ccpi1p_Whad_slice_config.txt "_w"
  add ccpi1p dphit3bin ccpi1p_dphit_bin_config_3bin.txt ccpi1p_dphit_slice_config_3bin.txt "_w"
fi
[ -n "${PP_K:-}" ] && [ "${RHC:-0}" != 1 ] && add ccpi1p pp ccpi1p_pp_bin_config.txt ccpi1p_pp_slice_config.txt "_w"
N=$(wc -l < "$M")
J=$(sbatch --parsable --cpus-per-task=4 --mem=12G --array=0-$((N-1)) --export=ALL,REPO=$PWD,MANIFEST=$PWD/$M,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
echo "$J $(basename $M) (alternative-model closure, $([ "${RHC:-0}" = 1 ] && echo RHC || echo FHC), binnings of the 0.50 criterion; $N extractions) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J ($N extractions)"
