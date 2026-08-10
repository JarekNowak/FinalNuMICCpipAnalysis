#!/usr/bin/env bash
# D'Agostini (iterative Bayesian) cross-check for ALL 3 configs x 5 observables,
# sweeping 1/2/4 iterations, on the same univmakes + fake data + per-config flux as
# Wiener-SVD. Confirms whether the per-observable integrated-sigma spread is the
# Wiener-SVD A_C smearing (D'Agostini truth should be flat, robust vs iterations).
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
TMP=/tmp/claude-661400007/-home-t2k-nowak-MicroBooNE-working-xsec-analyzer/e25e697c-fac6-48db-a45d-09ecd317057a/scratchpad
mkdir -p "$TMP"
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
for cfg in fhc5 rhcfull comb; do t=${TAG[$cfg]}
  for it in 1 2 4; do
    for o in pmu ppi costhmu costhpi thmupi; do
      c="$TMP/dag_${cfg}_${o}_it${it}.txt"
      sed "s|^Unfold DAgostini.*|Unfold DAgostini iter ${it}|" configs/ccpi_xsec_config_numi_${o}_${cfg}_dag.txt > "$c"
      bin/UnfolderNuMI "$c" configs/ccpi_${o}_slice_config_opt.txt \
        $PROC/xsec_${t}dag${it}_${o}.root > "../logs/dag_${cfg}_${o}_it${it}.log" 2>&1 || echo "  FAIL $cfg $o it$it"
    done
  done
  echo "  $cfg D'Agostini done"
done
echo "===== Wiener-SVD vs D'Agostini integrated sigma (truth, 1e-38 cm^2/Ar) ====="
root.exe -l -b -q -e '{
  const char* CF[3]={"FHC5","RHCFULL","COMB"}; const char* obs[5]={"pmu","ppi","costhmu","costhpi","thmupi"};
  auto ig=[&](const char* p,const char* h)->double{TFile*f=TFile::Open(p);if(!f||f->IsZombie())return 0;TH1*x=(TH1*)f->Get(h);double v=x?x->Integral("width"):0;f->Close();return v;};
  for(auto c:CF){ printf("\n[%s]  %-9s %8s | %8s %8s %8s %8s\n",c,"obs","WSVD_tr","DAg_tr","DAg_i1","DAg_i2","DAg_i4");
    for(auto o:obs){
      double wt=ig(Form("%s/closure_hists_xsec_%s_%s.root",PROC="/data/uboone/processed",c,o),"h_fakedata_truth");
      double dt=ig(Form("%s/closure_hists_xsec_%sdag4_%s.root",PROC,c,o),"h_fakedata_truth");
      double d1=ig(Form("%s/closure_hists_xsec_%sdag1_%s.root",PROC,c,o),"h_unfolded_nuwro");
      double d2=ig(Form("%s/closure_hists_xsec_%sdag2_%s.root",PROC,c,o),"h_unfolded_nuwro");
      double d4=ig(Form("%s/closure_hists_xsec_%sdag4_%s.root",PROC,c,o),"h_unfolded_nuwro");
      printf("       %-9s %8.3f | %8.3f %8.3f %8.3f %8.3f\n",o,wt,dt,d1,d2,d4);
    } }
}' 2>/dev/null
echo "########## DAGOSTINI ALL DONE ##########"
