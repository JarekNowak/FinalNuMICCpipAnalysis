#!/usr/bin/env bash
# run_altmodel_v1.sh -- alternative-model closure (prerequisite V1, 2026-09-27) for the results proposed for
# approval that were outside the earlier alternative-model round: the proton angles theta_p and theta_pi_p
# (FHC, like the other proton-tagged observables) and the two combined-only two-dimensional results,
# cos theta_pi x cos theta_mu (inclusive) and theta_p x delta p_T (proton-tagged). Pseudo-data: GENIE
# multisim universe 545 (altgenie) and the Delta->N pi angular unisim (altdelta), response nominal.
# PART=a: the angles and the inclusive 2D result; PART=b: the proton-tagged 2D result, which needs the RHC
# proton-tagged alternative pseudo-data (THROW_DIR=/data/uboone/processed/w/ macros/throw_altmodel_rhc.C).
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
RA=/data/uboone/processed/rebuild_alt; P=/data/uboone/processed
PART=${PART:-a}
M=slurm/altmodel_v1${PART}_manifest.list; : > "$M"
for v in altgenie altdelta; do
  # combined file lists with the alternative pseudo-data in place of the nominal throws
  sed -E "s#^($P/xsec-ana-fakedata)_(fhc|rhc)_(run[0-9])\.root#\1_${v}_\2_\3.root#" configs/file_properties_numi_comb.txt > configs/file_properties_numi_comb_${v}.txt
  sed -E "s#^($P/w/xsec-ana-fakedata)_(fhc|rhc)_(run[0-9])\.root#\1_${v}_\2_\3.root#" configs/file_properties_numi_comb_w.txt > configs/file_properties_numi_comb_w_${v}.txt
  for f in configs/file_properties_numi_comb_${v}.txt configs/file_properties_numi_comb_w_${v}.txt; do
    n=$(grep -c "fakedata_${v}_" "$f"); [ "$n" -eq 8 ] || { echo "substitution incomplete in $f ($n)"; exit 1; }
    for d in $(grep -o "^/[^ ]*fakedata_${v}_[^ ]*" "$f"); do [ -f "$d" ] || [ "$PART" = a -a "${f#*comb_w}" != "$f" ] || { echo "missing $d"; exit 1; }; done
  done
  if [ "$PART" = a ]; then
    for o in thetap thpipr; do
      out=configs/alt/ccpi1p_xsec_config_numi_${o}_${v}.txt
      sed -E "s#^UnivFile .*#UnivFile $RA/ccpi1p_FHC5_${o}_${v}_univmake.root#; s#^FPFile .*#FPFile configs/file_properties_numi_fhc5_w_${v}.txt#" configs/ccpi1p_xsec_config_numi_${o}_fhc5.txt > "$out"
      echo "$out configs/ccpi1p_${o}_bin_config.txt configs/ccpi1p_${o}_slice_config.txt" >> "$M"
    done
    out=configs/alt/ccpi_xsec_config_numi_costhpi_costhmu_2d_comb_${v}.txt
    sed -E "s#^UnivFile .*#UnivFile $RA/ccpi_COMB_costhpi_costhmu_${v}_univmake.root#; s#^FPFile .*#FPFile configs/file_properties_numi_comb_${v}.txt#" configs/d2/ccpi_xsec_config_numi_costhpi_costhmu_2d_comb.txt > "$out"
    echo "$out configs/d2/ccpi_costhpi_costhmu_bin_config_2d.txt configs/d2/ccpi_costhpi_costhmu_slice_config_2d.txt" >> "$M"
  else
    out=configs/alt/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb_${v}.txt
    sed -E "s#^UnivFile .*#UnivFile $RA/ccpi1p_COMB_thetap_dpt_${v}_univmake.root#; s#^FPFile .*#FPFile configs/file_properties_numi_comb_w_${v}.txt#" configs/d2/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb.txt > "$out"
    echo "$out configs/d2/ccpi1p_thetap_dpt_bin_config_2d.txt configs/d2/ccpi1p_thetap_dpt_slice_config_2d.txt" >> "$M"
  fi
done
N=$(wc -l < "$M")
J=$(sbatch --parsable --cpus-per-task=4 --mem=12G --array=0-$((N-1)) --export=ALL,REPO=$PWD,MANIFEST=$PWD/$M,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
echo "$J $(basename $M) (alternative-model closure V1 part $PART: $N extractions) $(date)" >> slurm/rebuild_jobs_2026-09-19.txt
echo "submitted $J ($N extractions)"
