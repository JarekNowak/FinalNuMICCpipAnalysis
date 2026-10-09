#!/bin/bash
# mp2pi_fluxw_univ.sh -- two-pion universe files and extractions with the detector-variation samples weighted to the
# new NuMI flux (scripts/detvar_fluxw.py, slurm/slurm_detvar_fluxw.sbatch). Instead of re-running univmake on all 46
# inputs, univmake runs on the 18 reweighted samples alone (they carry no universe weights, so only their unweighted
# universe is built) and macros/merge_univ_dirs.C adds their per-file directories to copies of the existing universe
# files, mp2pi_<CFG>_total_univmake.root -> mp2pi_<CFG>_total_fluxw_univmake.root. Then the extractions without and
# with the sideband constraint, DVfrac, into /data/uboone/processed/mp2pi/unfold/*_dvfrac_fluxw*.
#   bash scripts/mp2pi_fluxw_univ.sh > /data/uboone/processed/mp2pi/detvar_fluxw/logs/univ.log 2>&1
set -uo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$R"
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$R/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$R"
PROC=/data/uboone/processed; D=$PROC/mp2pi/detvar_fluxw; P=$PROC/mp2pi/unfold
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )

# 1. univmake on the 18 reweighted samples (its closing total over the full file list fails without the other
#    per-file directories; only the per-file directories are used)
grep -o "$D/xsec-ana-detvar_[^ ]*\.root" configs/file_properties_numi_comb_mp2pi_fluxw.txt > $D/univ_list.txt
n=$(wc -l < $D/univ_list.txt); [ "$n" -eq 18 ] || { echo "expected 18 reweighted samples, found $n"; exit 1; }
for f in $(cat $D/univ_list.txt); do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done
echo "==== univmake on $n samples $(date)"
UNIV_THREADS=${UNIV_THREADS:-3} nice -n 10 bin/univmake $D/univ_list.txt configs/mp2pi_total_bin_config.txt \
  $D/detvar_fluxw_univmake.root configs/file_properties_numi_comb_mp2pi_fluxw.txt > $D/logs/univmake.log 2>&1
echo "univmake rc=$? (the closing total is expected to fail)"
ndir=$(python3 -c "
import uproot; f=uproot.open('$D/detvar_fluxw_univmake.root'); r=f[f.keys(recursive=False)[0].split(';')[0]]
print(len({k.split(';')[0] for k in r.keys(recursive=False) if k.startswith('+')}))")
echo "per-file directories: $ndir"; [ "$ndir" -ge 18 ] || { echo "univmake output incomplete"; exit 1; }

# 2. merged universe files
for cfg in fhc5 rhcfull comb; do
  src=$PROC/mp2pi_${TAG[$cfg]}_total_univmake.root; dst=$PROC/mp2pi_${TAG[$cfg]}_total_fluxw_univmake.root
  cp -p "$src" "$dst.tmp" && root -l -b -q "macros/merge_univ_dirs.C(\"$D/detvar_fluxw_univmake.root\",\"$dst.tmp\")" \
    > $D/logs/merge_$cfg.log 2>&1 && grep -q "18 per-file directories copied under 1 root" $D/logs/merge_$cfg.log \
    && mv "$dst.tmp" "$dst" && echo "merged $cfg" || { echo "merge $cfg FAILED"; exit 1; }
done

# 3. extractions (unconstrained first: it builds the universe cache of the merged file)
for cfg in fhc5 rhcfull comb; do
  (
    W=$P/${cfg}_fluxw; mkdir -p $W; ln -sfn $R/configs $W/configs
    for mode in "" _constr; do
      ( cd $W && nice -n 10 $R/bin/UnfolderNuMI configs/mp2pi_xsec_config_numi_total_${cfg}${mode}_dvfrac_fluxw.txt \
          configs/mp2pi_total_slice_config.txt $P/xsec_mp2pi_${cfg}${mode}_dvfrac_fluxw_total.root ) \
          > $P/${cfg}${mode}_dvfrac_fluxw.log 2>&1
      echo "$cfg$mode rc=$?"
    done
  ) &
done
wait
for f in $P/*_dvfrac_fluxw.log; do echo "== $(basename $f)"; grep "\[SYSTDUMP\] \(sigma_int\|total \|flux_total\|xsec_total\|detVar_total\|DataStats\)" $f; done
echo "==== done $(date)"
