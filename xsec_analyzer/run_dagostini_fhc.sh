#!/usr/bin/env bash
# D'Agostini (iterative Bayesian) unfold for all 5 FHC observables, on the SAME
# univmake response matrices + fake data as Wiener-SVD. Sweep 1/2/4 iterations to
# check robustness. Cross-check: does the per-observable integrated-sigma spread
# persist without the Wiener-SVD A_C smearing, and is it stable vs iteration count?
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
TMP=/tmp/claude-661400007/-home-t2k-nowak-MicroBooNE-working-xsec-analyzer/e25e697c-fac6-48db-a45d-09ecd317057a/scratchpad
mkdir -p "$TMP"
for it in 1 2 4; do
  for o in pmu ppi costhmu costhpi thmupi; do
    cfg="$TMP/dag_${o}_it${it}.txt"
    sed "s|^Unfold DAgostini.*|Unfold DAgostini iter ${it}|" configs/ccpi_xsec_config_numi_${o}_fhc5_dag.txt > "$cfg"
    bin/UnfolderNuMI "$cfg" configs/ccpi_${o}_slice_config_opt.txt \
      $PROC/xsec_FHC5dag${it}_${o}.root > "../logs/fhc5_${o}_dagostini_it${it}.log" 2>&1 \
      && echo "  dag[$o it$it] done" || echo "  dag[$o it$it] FAIL"
  done
done
echo "===== Wiener-SVD vs D'Agostini integrated sigma (10^-38 cm^2/Ar) ====="
root.exe -l -b -q -e '{
  const char* obs[]={"pmu","ppi","costhmu","costhpi","thmupi"};
  const char* PROC="/data/uboone/processed";
  auto ig=[&](const char* path,const char* hn)->double{TFile*f=TFile::Open(path);if(!f||f->IsZombie())return 0;TH1*h=(TH1*)f->Get(hn);double v=h?h->Integral("width"):0;f->Close();return v;};
  printf("  %-9s %8s %8s | %8s | %8s %8s %8s  (DAg unfolded vs iterations)\n","obs","WSVD_un","WSVD_tr","DAg_tr","DAg_i1","DAg_i2","DAg_i4");
  for(auto o:obs){
    double wu=ig(Form("%s/closure_hists_xsec_FHC5_%s.root",PROC,o),"h_unfolded_nuwro");
    double wt=ig(Form("%s/closure_hists_xsec_FHC5_%s.root",PROC,o),"h_fakedata_truth");
    double dt=ig(Form("%s/closure_hists_xsec_FHC5dag4_%s.root",PROC,o),"h_fakedata_truth");
    double d1=ig(Form("%s/closure_hists_xsec_FHC5dag1_%s.root",PROC,o),"h_unfolded_nuwro");
    double d2=ig(Form("%s/closure_hists_xsec_FHC5dag2_%s.root",PROC,o),"h_unfolded_nuwro");
    double d4=ig(Form("%s/closure_hists_xsec_FHC5dag4_%s.root",PROC,o),"h_unfolded_nuwro");
    printf("  %-9s %8.3f %8.3f | %8.3f | %8.3f %8.3f %8.3f\n",o,wu,wt,dt,d1,d2,d4);
  }}' 2>/dev/null | grep -E "obs|pmu|ppi|costh|thmupi"
echo "===== DAGOSTINI COMPARE DONE ====="
