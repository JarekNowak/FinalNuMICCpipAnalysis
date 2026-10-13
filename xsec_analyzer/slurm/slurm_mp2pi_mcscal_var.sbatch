#!/bin/bash
#SBATCH --job-name=mp2pivar
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# SLURM array task: calibration uncertainty of the two-pion one-bin total for ONE configuration and ONE variant MCS
# map, as slurm/slurm_mcscal_var.sbatch for the single-pion results. Task i: configuration i / 7 of (fhc5 rhcfull
# comb), variant i % 7 of VARIANTS. The calibrated universe file (mp2pi_<CFG>_total_mcscal_fluxw_univmake.root) is
# copied; univmake runs on the pseudo-data of the variant alone (scripts/mp2pi_mcscal_var_throw.sh), their per-file
# directories are added to the copy (macros/merge_univ_dirs.C), and both extractions (without and with the sideband
# constraint) run with the pseudo-data of the variant against the response calibrated with the nominal map.
# Output: /data/uboone/processed/mp2pi_mcscal/var/<variant>/unfold/<cfg><mode>.log (SYSTDUMP).
#   sbatch --array=0-20 --dependency=afterok:<throws> --export=ALL,REPO=$PWD slurm/slurm_mp2pi_mcscal_var.sbatch
set -u
REPO="${REPO:?}"; N=/data/uboone/processed/mp2pi_mcscal
VARIANTS=(${MP2PI_VARIANTS:-shrink coarse fine old_production new_production FHC RHC}); CFGS=(fhc5 rhcfull comb)   # MP2PI_VARIANTS: test
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
CFG=${CFGS[$(( SLURM_ARRAY_TASK_ID / ${#VARIANTS[@]} ))]}; V=${VARIANTS[$(( SLURM_ARRAY_TASK_ID % ${#VARIANTS[@]} ))]}
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
echo "host=$(hostname) cfg=$CFG variant=$V"
U=/data/uboone/processed/mp2pi_${TAG[$CFG]}_total_mcscal_fluxw_univmake.root
FP=configs/file_properties_numi_${CFG}_mp2pi_mcscal_fluxw.txt; BIN=configs/mp2pi_total_bin_config.txt
D="$N/var/$V/unfold"; WK="$N/var/work"; mkdir -p "$D" "$WK"
UV="$WK/${CFG}_${V}_univmake.root"
cp -p "$U" "$UV.tmp" && mv "$UV.tmp" "$UV" || { echo "FAIL copy" >&2; exit 1; }
FPV="$WK/${CFG}_${V}_fp.txt"; sed "s|^$N/xsec-ana-fakedata_|$N/var/$V/xsec-ana-fakedata_|" "$FP" > "$FPV"
NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FP" | wc -l); NV=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPV" | grep -c "^$N/var/$V/")
[ "$NON" -gt 0 ] && [ "$NON" = "$NV" ] || { echo "FAIL file list: $NON onBNB, $NV variant" >&2; exit 1; }
LIST="$WK/${CFG}_${V}_fakedata.list"; awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPV" > "$LIST"
SMALL="$WK/${CFG}_${V}_small_univmake.root"; command rm -f -- "$SMALL"
UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$LIST" "$BIN" "$SMALL" "$FPV" > "$WK/${CFG}_${V}_univmake.log" 2>&1
echo "univmake rc=$? (the closing total is expected to fail)"
root -l -b -q "macros/merge_univ_dirs.C(\"$SMALL\",\"$UV\")" > "$WK/${CFG}_${V}_merge.log" 2>&1
grep -q "^$NON per-file directories copied under 1 root" "$WK/${CFG}_${V}_merge.log" || { echo "FAIL merge" >&2; cat "$WK/${CFG}_${V}_merge.log"; exit 1; }
W="$WK/${CFG}_${V}_unfold"; command rm -rf -- "$W"; mkdir -p "$W/unfold_output"; ln -s "$REPO/configs" "$W/configs"
ok=1
for mode in "" _constr; do
  XC=configs/mp2pi_xsec_config_numi_total_${CFG}${mode}_dvfrac_fluxw_mcscal.txt; XV="$WK/${CFG}${mode}_${V}_xsec.txt"
  sed -e "s|^UnivFile .*|UnivFile $UV|" -e "s|^FPFile .*|FPFile $FPV|" "$XC" > "$XV"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "$XV" configs/mp2pi_total_slice_config.txt "$D/xsec_mp2pi_${CFG}${mode}_${V}_total.root" ) \
    > "$D/${CFG}${mode}.log" 2>&1
  rc=$?; s=$(grep 'SYSTDUMP] sigma_int' "$D/${CFG}${mode}.log" | awk '{print $3}')
  echo "$CFG$mode rc=$rc sigma_int=$s"; [ "$rc" = 0 ] && [ -n "$s" ] || ok=0
done
command rm -rf -- "$W" "$UV" "$SMALL"
[ $ok = 1 ] && echo "DONE $CFG $V" || { echo "FAIL $CFG $V"; exit 1; }
