#!/bin/bash
# mp2pi_extract.sh -- the two-pion one-bin total through the framework (Phase 4, item 3, of
# report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md): for each configuration the extraction without and with the
# sideband constraint, from the universe files of slurm/slurm_univmake_mp2pi.sbatch. Unconstrained first (it builds
# the universe cache), then constrained; the configurations in parallel, each in its own work directory.
R="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P=/data/uboone/processed/mp2pi/unfold
export LD_LIBRARY_PATH="$R/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$R"
for cfg in fhc5 rhcfull comb; do
  (
    W=$P/$cfg; mkdir -p $W; ln -sfn $R/configs $W/configs
    for mode in "" _constr; do
      ( cd $W && nice -n 10 $R/bin/UnfolderNuMI configs/mp2pi_xsec_config_numi_total_${cfg}${mode}.txt \
          configs/mp2pi_total_slice_config.txt $P/xsec_mp2pi_${cfg}${mode}_total.root ) > $P/${cfg}${mode}.log 2>&1
      echo "$cfg${mode} rc=$?"
    done
  ) &
done
wait
for f in $P/*.log; do echo "== $(basename $f)"; grep "\[SYSTDUMP\] \(sigma_int\|total \|flux_total\|xsec_total\|detVar_total\|reint\|DataStats\|SimulationStats\|POT\|numTargets\)" $f; done
