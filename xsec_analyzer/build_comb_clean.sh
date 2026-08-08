#!/usr/bin/env bash
# Combined FHC+RHC build with PER-MODE summed_pot (numuMC + detVar) so each mode
# scales to its own data exposure. Run ONLY after FHC/RHC builds are done (shares files).
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed

echo "[setpot] per-mode: FHC numuMC=2.157846e22, RHC numuMC=2.782194e22; detVar native*mode-factor"
root.exe -l -b -q -e '{
  // FHC numuMC
  const char* fF[]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
  for(auto r:fF){TFile f(Form("/data/uboone/processed/xsec-ana-%s.root",r),"UPDATE");TParameter<float> s("summed_pot",2.157846e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  // RHC numuMC
  const char* fR[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
  for(auto r:fR){TFile f(Form("/data/uboone/processed/xsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",r),"UPDATE");TParameter<float> s("summed_pot",2.782194e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  for(auto s3:{"aa","ab","ac","ad","ae"}){TFile f(Form("/data/uboone/processed/xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",s3),"UPDATE");TParameter<float> s("summed_pot",2.782194e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  printf("numuMC set\n");
}' 2>/dev/null | grep set
# detVar per-mode: summed_pot = native * (FHC 2.25121 | RHC 1.79922), read native snapshot
while read tag native; do
  [[ "$tag" =~ ^# || -z "$tag" ]] && continue
  mult=$([[ "$tag" == run4fhc* ]] && echo 2.25121 || echo 1.79922)
  root.exe -l -b -q -e "{double v=${native}*${mult}; TFile f(\"/data/uboone/processed/xsec-ana-detvar_${tag}.root\",\"UPDATE\"); TParameter<float> s(\"summed_pot\",(float)v); s.Write(\"summed_pot\",TObject::kOverwrite); f.Close();}" 2>/dev/null
done < configs/detvar_native_pot.txt
echo "  detVar per-mode summed_pot set"

echo "[throw] combined-full fake data (seed 1)..."; root.exe -l -b -q "macros/throw_cv_comb_full.C(1)" 2>/dev/null | grep "COMBFULL seed"
echo "[build] combined univmakes (5 observables)..."
pids=(); for o in pmu ppi costhmu costhpi thmupi; do
  FPM=configs/file_properties_numi_comb.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt OUT=$PROC/ccpi_COMB_${o}_univmake.root ./run_universe_maker.sh > "../logs/comb_${o}_univ.log" 2>&1 & pids+=($!)
done
fail=0; i=0; for o in pmu ppi costhmu costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || { echo "  univmake[$o] FAIL"; fail=1; }; i=$((i+1)); done
echo "===== COMBINED FHC+RHC FLUX + CLOSURE (detVar=Run4 FHC+RHC per-mode) ====="
for o in pmu ppi costhmu costhpi thmupi; do
  fb=$(/tmp/flux_breakdown configs/ccpi_xsec_config_numi_${o}_comb.txt "$o" 2>/dev/null | grep "flux RECO")
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_comb.txt configs/ccpi_${o}_slice_config_opt.txt /data/uboone/processed/xsec_COMB_${o}.root > "../logs/comb_${o}_unfold.log" 2>&1
  cl=$(grep -iE "^truth: " "../logs/comb_${o}_unfold.log" | head -1)
  echo "--- $o ---"; echo "  $fb"; echo "  closure $cl"
done
echo "===== COMB BUILD DONE (fail=$fail) ====="
