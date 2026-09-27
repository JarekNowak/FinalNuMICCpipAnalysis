#!/usr/bin/env bash
# run_mcs_1p.sh -- the data-side MCS muon-momentum scale (+-5% on the MCS-measured muons of the fake data,
# nominal response) for the proton-tagged observables that read the muon momentum: p_mu, W_had, delta p_T,
# delta alpha_T and p_n, in FHC, RHC and combined (2026-09-27). delta phi_T reads only the muon direction
# and W_pipr, p_p, p_pi and the angles no muon at all, so their term is exactly zero.
# Fake data: w/xsec-ana-fakedata_<beam>_run<i>_mcs{up05,dn05}.root from macros/mcs_scale_fakedata_1p.C,
# which recomputes the TKI and W_had from the scaled muon (s=1 reproduces the stored branches exactly).
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
RA=/data/uboone/processed/rebuild_alt; mkdir -p "$RA" configs/mcs
M=slurm/mcsdata_1p_manifest.list; : > "$M"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
declare -A BIN=( [pmu]=ccpi1p_pmu_bin_config.txt [Whad]=ccpi1p_Whad_bin_config.txt [dpt2bin]=ccpi1p_dpt_bin_config_2bin.txt
                 [dalphat2bin]=ccpi1p_dalphat_bin_config_2bin.txt [pn2bin]=ccpi1p_pn_bin_config_2bin.txt )
for c in fhc5 rhcfull comb; do T=${TAG[$c]}
  for v in up05 dn05; do
    fp=configs/mcs/file_properties_numi_${c}_w_mcsdata_${v}.txt
    sed -E "s#^(/data/uboone/processed/w/xsec-ana-fakedata_(fhc|rhc)_run[0-9])\.root#\1_mcs${v}.root#" \
      configs/file_properties_numi_${c}_w.txt > "$fp"
    [ "$(grep -c "_mcs${v}.root" "$fp")" -eq "$(grep -c '^/data/uboone/processed/w/xsec-ana-fakedata_' configs/file_properties_numi_${c}_w.txt)" ] \
      || { echo "fake-data substitution incomplete in $fp"; exit 1; }
    for o in pmu Whad dpt2bin dalphat2bin pn2bin; do
      out=configs/mcs/ccpi1p_xsec_config_numi_${o}_${c}_mcsdata_${v}.txt
      sed -E "s#^UnivFile .*#UnivFile $RA/ccpi1p_${T}_${o}_mcsdata${v}_univmake.root#; s#^FPFile .*#FPFile $fp#" \
        configs/ccpi1p_xsec_config_numi_${o}_${c}.txt > "$out"
      b=configs/${BIN[$o]}; s=${b/_bin_config/_slice_config}
      echo "$out $b $s" >> "$M"
    done
  done
done
for f in $(awk '{print $2, $3}' "$M"); do [ -f "$f" ] || { echo "missing $f"; exit 1; }; done
N=$(wc -l < "$M")
J=$(sbatch --parsable --cpus-per-task=4 --mem=12G --array=0-$((N-1)) --export=ALL,REPO=$PWD,MANIFEST=$PWD/$M,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
echo "$J $(basename $M) (data-side MCS scale +-5% for the proton-tagged p_mu, W_had, dpt2bin, dalphat2bin, pn2bin; FHC/RHC/COMB) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J ($N tasks)"
