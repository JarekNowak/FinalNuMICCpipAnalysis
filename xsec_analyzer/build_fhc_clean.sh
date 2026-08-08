#!/usr/bin/env bash
# Clean FHC {Run1,2,4,5} build: setpot -> throw -> univmake(5) -> unfold+flux+closure.
# detVar now = Run4 FHC native (stand-in). No stale waits; all inputs ready.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
echo "[setpot] FHC combined 9.585254e21 on Run1,2,4 + Run5(reweightedPPFX)"
root.exe -l -b -q -e '{const char* fs[]={"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root"}; for(auto r:fs){TFile f(Form("/data/uboone/processed/%s",r),"UPDATE");TParameter<float> sp("summed_pot",9.585254e21f);sp.Write("summed_pot",TObject::kOverwrite);f.Close();} printf("summed_pot set\n");}' 2>/dev/null | grep set
echo "[throw] FHC fake data (seed 1)..."; root.exe -l -b -q "macros/throw_cv_fhc5.C(1)" 2>/dev/null | grep "FHC5 seed"
echo "[build] FHC univmakes (5 observables)..."
pids=(); for o in pmu ppi costhmu costhpi thmupi; do
  FPM=configs/file_properties_numi_fhc5.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt OUT=$PROC/ccpi_FHC5_${o}_univmake.root ./run_universe_maker.sh > "../logs/fhc5_${o}_univ.log" 2>&1 & pids+=($!)
done
fail=0; i=0; for o in pmu ppi costhmu costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || { echo "  univmake[$o] FAIL"; fail=1; }; i=$((i+1)); done
echo "===== FHC {Run1,2,4,5} FLUX + CLOSURE (detVar=Run4 FHC) ====="
for o in pmu ppi costhmu costhpi thmupi; do
  fb=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_fhc5.txt "$o" 2>/dev/null | grep "flux RECO")
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_fhc5.txt configs/ccpi_${o}_slice_config_opt.txt /data/uboone/processed/xsec_FHC5_${o}.root > "../logs/fhc5_${o}_unfold.log" 2>&1
  cl=$(grep -iE "^truth: " "../logs/fhc5_${o}_unfold.log" | head -1)
  dv=$(grep -iE "detVar_total" "../logs/fhc5_${o}_unfold.log" | head -1)
  echo "--- $o ---"; echo "  $fb"; echo "  closure $cl"
done
echo "===== FHC BUILD DONE (fail=$fail) ====="
