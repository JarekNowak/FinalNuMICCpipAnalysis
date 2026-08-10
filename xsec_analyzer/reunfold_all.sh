#!/usr/bin/env bash
# Re-unfold ALL configs x observables with the rebuilt binary (per-config Flux) so
# every result uses its correct flux normalisation. Reads existing univmakes (which
# already baked in the per-mode POT pooling) -> setpot-independent, fast. FHC flux
# unchanged (identical results); RHC (+5.7%) and combined (+3.1%) corrected.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
for cfg in fhc5 rhcfull comb; do
  t=${TAG[$cfg]}
  for o in pmu ppi costhmu costhpi thmupi; do
    bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_${cfg}.txt configs/ccpi_${o}_slice_config_opt.txt \
      $PROC/xsec_${t}_${o}.root > "../logs/reunfold_${cfg}_${o}.log" 2>&1 \
      && echo "  unfold[$cfg $o] ok" || echo "  unfold[$cfg $o] FAIL"
  done
done
echo "===== integrated sigma, all configs (Wiener-SVD, corrected per-config flux) ====="
root.exe -l -b -q -e '{
  const char* cf[3]={"FHC5","RHCFULL","COMB"}; const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  printf("  %-8s %8s %8s %8s %8s %8s   (truth integral, 1e-38 cm^2/Ar)\n","config",obs[0],obs[1],obs[2],obs[3],obs[4]);
  for(auto c:cf){ printf("  %-8s",c);
    for(auto o:obs){ TFile*f=TFile::Open(Form("/data/uboone/processed/closure_hists_xsec_%s_%s.root",c,o));
      double t=0; if(f&&!f->IsZombie()){TH1*h=(TH1*)f->Get("h_fakedata_truth"); if(h)t=h->Integral("width"); f->Close();}
      printf(" %8.3f",t);} printf("\n"); }
}' 2>/dev/null | grep -E "config|FHC5|RHCFULL|COMB"
echo "########## REUNFOLD ALL DONE ##########"
