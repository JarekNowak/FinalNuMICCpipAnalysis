#!/bin/bash
#SBATCH --job-name=kinxsec
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=06:00:00
# SLURM job: the extractions of one variable of the two-pion differential measurement (configs of
# scripts/mp2pi_kin_configs.py, universes of slurm/slurm_univmake_mp2pi.sbatch with MP2PI_SUFFIX=_kin_fluxw and
# MP2PI_OUTTAG=_<variable>): FHC, RHC and combined, without and then with the sideband constraint, on the pseudo-data.
# Output: /data/uboone/processed/mp2pi_kin/unfold/<variable>/{xsec_<cfg><mode>.root,<cfg><mode>.log} and the closure
# sidecars.
#   sbatch --export=ALL,REPO=$PWD,VAR=pmu slurm/slurm_mp2pi_kin_extract.sbatch
set -u
REPO="${REPO:?}"; VAR="${VAR:?}"; PROD="${PROD:-}"; P=/data/uboone/processed/mp2pi_${PROD:-kin}/unfold/$VAR   # PROD=2f: the two-fold production (VAR=total: the one-bin total)
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
mkdir -p "$P"; ok=1
for cfg in fhc5 rhcfull comb; do
  W=$P/$cfg; mkdir -p "$W"; ln -sfn "$REPO/configs" "$W/configs"
  for mode in "" _constr; do
    ( cd "$W" && "$REPO/bin/UnfolderNuMI" configs/mp2pi_xsec_config_numi_kin${PROD}_${VAR}_${cfg}${mode}.txt \
        configs/mp2pi_kin_${VAR}_slice_config.txt "$P/xsec_${cfg}${mode}.root" ) > "$P/${cfg}${mode}.log" 2>&1
    rc=$?; echo "$VAR $cfg$mode rc=$rc total=$(grep 'SYSTDUMP] total ' "$P/${cfg}${mode}.log" | awk '{print $3}')"
    [ $rc = 0 ] || ok=0
  done
done
[ $ok = 1 ] && echo "DONE $VAR" || { echo "FAIL $VAR"; exit 1; }
