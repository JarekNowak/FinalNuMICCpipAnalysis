#!/usr/bin/env bash
# Re-run FHC univmake+unfold for costhpi+thmupi ONLY, at the coarsened binning.
# Reuses the existing fake-data throw + summed_pot (both binning-independent).
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
echo "[univmake] costhpi + thmupi (coarsened) in parallel..."
pids=()
for o in costhpi thmupi; do
  FPM=configs/file_properties_numi_fhc5.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT=$PROC/ccpi_FHC5_${o}_univmake.root ./run_universe_maker.sh > "../logs/fhc5_${o}_univ_coarse.log" 2>&1 & pids+=($!)
done
i=0; for o in costhpi thmupi; do wait "${pids[$i]}" && echo "  univmake[$o] OK" || echo "  univmake[$o] FAIL"; i=$((i+1)); done
echo "[unfold] costhpi + thmupi..."
for o in costhpi thmupi; do
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_fhc5.txt configs/ccpi_${o}_slice_config_opt.txt \
    $PROC/xsec_FHC5_${o}.root > "../logs/fhc5_${o}_unfold_coarse.log" 2>&1
  echo "  unfold[$o] done"
done
echo "===== COARSENED FHC integrated sigma (unfolded / truth) ====="
root.exe -l -b -q -e '{const char* obs[]={"costhpi","thmupi"};for(auto o:obs){TFile*f=TFile::Open(Form("/data/uboone/processed/closure_hists_xsec_FHC5_%s.root",o));if(!f||f->IsZombie()){printf("  %s NO FILE\n",o);continue;}TH1*u=(TH1*)f->Get("h_unfolded_nuwro");TH1*t=(TH1*)f->Get("h_fakedata_truth");printf("  %-8s nbins=%d unfolded=%.3f truth=%.3f\n",o,u->GetNbinsX(),u->Integral("width"),t->Integral("width"));f->Close();}}' 2>/dev/null | grep -E "costhpi|thmupi"
echo "===== closure chi2 ====="
for o in costhpi thmupi; do grep -iE "^truth: " "../logs/fhc5_${o}_unfold_coarse.log" | head -1 | sed "s/^/  $o: /"; done
echo "===== COARSEN TEST DONE ====="
