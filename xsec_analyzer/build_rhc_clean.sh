#!/usr/bin/env bash
# Clean RHC-full {Run1,2,3,4a,4b,4c} build: setpot -> throw -> univmake(5) -> unfold+flux+closure.
# detVar now = Run4 RHC native. No stale waits; all inputs ready. Disjoint files from FHC -> parallel-safe.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
echo "[setpot] RHC combined 1.546330e22 on 10 RHC files"
root.exe -l -b -q -e '{const char* rhc[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; for(auto r:rhc){TFile f(Form("/data/uboone/processed/xsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",r),"UPDATE");TParameter<float> s("summed_pot",1.546330e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();} for(auto s3:{"aa","ab","ac","ad","ae"}){TFile f(Form("/data/uboone/processed/xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",s3),"UPDATE");TParameter<float> s("summed_pot",1.546330e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();} printf("summed_pot set\n");}' 2>/dev/null | grep set
echo "[throw] RHC fake data (seed 1)..."; root.exe -l -b -q "macros/throw_cv_rhcfull.C(1)" 2>/dev/null | grep "RHCFULL seed"
echo "[build] RHC univmakes (5 observables)..."
pids=(); for o in pmu ppi costhmu costhpi thmupi; do
  FPM=configs/file_properties_numi_rhcfull.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt OUT=$PROC/ccpi_RHCFULL_${o}_univmake.root ./run_universe_maker.sh > "../logs/rhcfull_${o}_univ.log" 2>&1 & pids+=($!)
done
fail=0; i=0; for o in pmu ppi costhmu costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || { echo "  univmake[$o] FAIL"; fail=1; }; i=$((i+1)); done
echo "===== RHC-FULL {Run1,2,3,4a,4b,4c} FLUX + CLOSURE (detVar=Run4 RHC) ====="
for o in pmu ppi costhmu costhpi thmupi; do
  fb=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_rhcfull.txt "$o" 2>/dev/null | grep "flux RECO")
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_rhcfull.txt configs/ccpi_${o}_slice_config_opt.txt /data/uboone/processed/xsec_RHCFULL_${o}.root > "../logs/rhcfull_${o}_unfold.log" 2>&1
  cl=$(grep -iE "^truth: " "../logs/rhcfull_${o}_unfold.log" | head -1)
  echo "--- $o ---"; echo "  $fb"; echo "  closure $cl"
done
echo "===== RHC BUILD DONE (fail=$fail) ====="
