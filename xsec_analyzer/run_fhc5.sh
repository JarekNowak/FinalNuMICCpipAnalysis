#!/usr/bin/env bash
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed; mkdir -p unfold_output
echo "[wait] Run5 numu processing..."; for i in $(seq 1 240); do grep -q "PROCESS RUN5NUMU DONE" ../logs/process_run5numu.log 2>/dev/null && break; sleep 30; done
grep -E "Run5numu OK|Run5numu FAILED" ../logs/process_run5numu.log
echo "[wait] Phase3 to finish (frees shared FHC summed_pot)..."; for i in $(seq 1 400); do grep -q "COMB ALLDONE" ../logs/run_comb.log 2>/dev/null && break; sleep 30; done
echo "[setpot] combined 9.585254e21 in Run1,2,4 + Run5(reweightedPPFX)"
root.exe -l -b -q -e '{const char* fs[]={"xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root","xsec-ana-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root"}; for(auto r:fs){TFile f(Form("/data/uboone/processed/%s",r),"UPDATE");TParameter<float> sp("summed_pot",9.585254e21f);sp.Write("summed_pot",TObject::kOverwrite);f.Close();} printf("summed_pot set\n");}' 2>/dev/null | grep set
echo "[throw] FHC5 fake data..."; root.exe -l -b -q "macros/throw_cv_fhc5.C(1)" 2>/dev/null | grep "FHC5 seed"
echo "[build] FHC5 univmakes..."
pids=(); for o in pmu ppi costhmu costhpi thmupi; do
  FPM=configs/file_properties_numi_fhc5.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt OUT=$PROC/ccpi_FHC5_${o}_univmake.root ./run_universe_maker.sh > "../logs/fhc5_${o}_univ.log" 2>&1 & pids+=($!)
done
i=0; for o in pmu ppi costhmu costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || echo "  univmake[$o] FAIL"; i=$((i+1)); done
echo "===== FHC {Run1,2,4,5} FLUX BREAKDOWN + CLOSURE ====="
for o in pmu ppi costhmu costhpi thmupi; do
  fb=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_fhc5.txt "$o" 2>/dev/null | grep "flux RECO")
  cl=$(bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_fhc5.txt configs/ccpi_${o}_slice_config_opt.txt /tmp/fhc5_${o}.root 2>/dev/null | grep -iE "^truth: " | head -1)
  echo "--- $o ---"; echo "  $fb"; echo "  closure $cl"
done
echo "===== FHC5 ALLDONE ====="
