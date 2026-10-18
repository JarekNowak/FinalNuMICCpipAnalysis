#!/bin/bash
# mp2pi_2f_mcsvar_submit.sh -- MCS-calibration term of the two-pion one-bin total with the two-fold event classifier
# (production /data/uboone/processed/mp2pi_2f), as for the calibrated single-classifier production: for each of the
# seven variant maps the overlays and dirt of the pseudo-data plan are processed with that map
# (slurm/slurm_braggfix_reprocess.sbatch, the job of the nominal two-fold production, into var/<variant>), thrown with
# the nominal plan (scripts/mp2pi_mcscal_var_throw.sh with MP2PI_DIR) and extracted against the nominal response
# (slurm/slurm_mp2pi_mcscal_var.sbatch with PROD=2f). Then: python3 scripts/mp2pi_mcscal_var_compare.py --prod 2f
#   bash scripts/mp2pi_2f_mcsvar_submit.sh
set -euo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$R"
N=/data/uboone/processed/mp2pi_2f; MAPS=/data/uboone/processed/mcscal/maps
VARIANTS=(shrink coarse fine old_production new_production FHC RHC)
THROWS=()
for V in "${VARIANTS[@]}"; do
  D=$N/var/$V; mkdir -p "$D/logs"
  awk '$1=="IN"{print $2}' $N/throw_plan.txt | sort -u > "$D/plan_in.txt"
  awk -F'|' -v d="$D/" 'NR==FNR{a[$1];next} ($4 in a){n=split($4,p,"/"); $4=d p[n]; print}' OFS='|' \
      "$D/plan_in.txt" slurm/mp2pi_2f_manifest.list > "$D/manifest.list"
  n=$(wc -l < "$D/manifest.list"); [ "$n" = "$(wc -l < "$D/plan_in.txt")" ] || { echo "manifest $V: $n lines"; exit 1; }
  P=$(sbatch --parsable --chdir="$D/logs" --job-name=mcsvar2f_$V --mem=6G --time=08:00:00 --array=0-$((n - 1))%5 \
      --export=ALL,REPO=$R,MANIFEST=$D/manifest.list,CC1MU1PIXP_MP_MCS_MAP=$MAPS/$V.txt slurm/slurm_braggfix_reprocess.sbatch)
  T=$(sbatch --parsable --chdir="$D/logs" --job-name=mcsthrow2f_$V --mem=4G --time=06:00:00 --dependency=afterok:$P \
      --export=ALL,MP2PI_DIR=$N --wrap="bash $R/scripts/mp2pi_mcscal_var_throw.sh $V")
  echo "$V: production $P ($n samples), throw $T"; THROWS+=("$T")
done
E=$(sbatch --parsable --chdir=$N/var --job-name=mcsext2f --array=0-20%7 --dependency=afterok:$(IFS=:; echo "${THROWS[*]}") \
    --export=ALL,REPO=$R,PROD=2f slurm/slurm_mp2pi_mcscal_var.sbatch)
echo "extraction $E"
