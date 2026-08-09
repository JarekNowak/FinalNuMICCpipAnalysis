#!/usr/bin/env bash
# Extend the coarsened costhpi/thmupi binning to RHC + combined (targeted: only these
# two observables; pmu/ppi/costhmu univmakes are unchanged and kept). Runs AFTER the
# diagnostics reprocess, which reset numuMC summed_pot to native. detVar files are
# currently at comb per-mode values, so:
#   RHC  : reset run4rhc detVar -> native, setpot RHC numuMC -> 1.546330e22
#   comb : setpot detVar -> per-mode, setpot numuMC -> per-mode
# RHC and comb SHARE detVar files with different summed_pot -> must run sequentially.
# Fake-data throws are event-level/deterministic (seed 1) and NOT re-thrown.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
export XSEC_ANALYZER_DIR="$PWD"; PROC=/data/uboone/processed
SNAP=configs/detvar_native_pot.txt

setpot_detvar() {  # $1 = fhc_mult, $2 = rhc_mult (empty -> skip that mode; use "native" to write native as-is)
  local fm="$1" rm="$2"
  while read tag native; do
    [[ "$tag" =~ ^# || -z "$tag" ]] && continue
    local mult=""
    [[ "$tag" == run4fhc* ]] && mult="$fm"
    [[ "$tag" == run4rhc* ]] && mult="$rm"
    [[ -z "$mult" ]] && continue
    local val; if [[ "$mult" == "native" ]]; then val="$native"; else val=$(awk "BEGIN{printf \"%.6e\", $native*$mult}"); fi
    root.exe -l -b -q -e "{TFile f(\"$PROC/xsec-ana-detvar_${tag}.root\",\"UPDATE\");TParameter<float> s(\"summed_pot\",(float)${val});s.Write(\"summed_pot\",TObject::kOverwrite);f.Close();}" 2>/dev/null
  done < "$SNAP"
}

# ======================= RHC =======================
echo "[RHC] reset run4rhc detVar -> native"
setpot_detvar "" native
echo "[RHC] setpot numuMC -> 1.546330e22"
root.exe -l -b -q -e '{const char* rhc[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"}; for(auto r:rhc){TFile f(Form("/data/uboone/processed/xsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",r),"UPDATE");TParameter<float> s("summed_pot",1.546330e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();} for(auto s3:{"aa","ab","ac","ad","ae"}){TFile f(Form("/data/uboone/processed/xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",s3),"UPDATE");TParameter<float> s("summed_pot",1.546330e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();} printf("rhc numuMC set\n");}' 2>/dev/null | grep set
echo "[RHC] univmake costhpi + thmupi (coarsened)..."
pids=(); for o in costhpi thmupi; do
  FPM=configs/file_properties_numi_rhcfull.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT=$PROC/ccpi_RHCFULL_${o}_univmake.root ./run_universe_maker.sh > "../logs/rhcfull_${o}_univ_coarse.log" 2>&1 & pids+=($!)
done
i=0; for o in costhpi thmupi; do wait "${pids[$i]}" && echo "  RHC univmake[$o] OK" || echo "  RHC univmake[$o] FAIL"; i=$((i+1)); done
for o in costhpi thmupi; do
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_rhcfull.txt configs/ccpi_${o}_slice_config_opt.txt \
    $PROC/xsec_RHCFULL_${o}.root > "../logs/rhcfull_${o}_unfold_coarse.log" 2>&1; echo "  RHC unfold[$o] done"
done

# ======================= comb =======================
echo "[comb] setpot detVar per-mode (fhc x2.25121, rhc x1.79922)"
setpot_detvar 2.25121 1.79922
echo "[comb] setpot numuMC per-mode"
root.exe -l -b -q -e '{
  const char* fF[]={"Run1_fhc_new_numi_flux_fhc_pandora_ntuple","Run2_fhc_new_numi_flux_fhc_pandora_ntuple","Run4_fhc_new_numi_flux_fhc_pandora_ntuple","reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc"};
  for(auto r:fF){TFile f(Form("/data/uboone/processed/xsec-ana-%s.root",r),"UPDATE");TParameter<float> s("summed_pot",2.157846e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  const char* fR[]={"Run1_rhc","Run2_rhc","Run4a_rhc","Run4b_rhc","Run4c_rhc"};
  for(auto r:fR){TFile f(Form("/data/uboone/processed/xsec-ana-%s_new_numi_flux_rhc_pandora_ntuple.root",r),"UPDATE");TParameter<float> s("summed_pot",2.782194e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  for(auto s3:{"aa","ab","ac","ad","ae"}){TFile f(Form("/data/uboone/processed/xsec-ana-Run3_rhc_new_numi_flux_rhc_pandora_ntuple_%s.root",s3),"UPDATE");TParameter<float> s("summed_pot",2.782194e22f);s.Write("summed_pot",TObject::kOverwrite);f.Close();}
  printf("comb numuMC set\n");}' 2>/dev/null | grep set
echo "[comb] univmake costhpi + thmupi (coarsened)..."
pids=(); for o in costhpi thmupi; do
  FPM=configs/file_properties_numi_comb.txt BIN_CONFIG=configs/ccpi_${o}_bin_config_opt.txt \
    OUT=$PROC/ccpi_COMB_${o}_univmake.root ./run_universe_maker.sh > "../logs/comb_${o}_univ_coarse.log" 2>&1 & pids+=($!)
done
i=0; for o in costhpi thmupi; do wait "${pids[$i]}" && echo "  comb univmake[$o] OK" || echo "  comb univmake[$o] FAIL"; i=$((i+1)); done
for o in costhpi thmupi; do
  bin/UnfolderNuMI configs/ccpi_xsec_config_numi_${o}_comb.txt configs/ccpi_${o}_slice_config_opt.txt \
    $PROC/xsec_COMB_${o}.root > "../logs/comb_${o}_unfold_coarse.log" 2>&1; echo "  comb unfold[$o] done"
done

# ======================= report =======================
echo "===== RHC + comb coarsened integrated sigma (unfolded / truth) ====="
root.exe -l -b -q -e '{const char* cf[]={"RHCFULL","COMB"};const char* obs[]={"costhpi","thmupi"};for(auto c:cf)for(auto o:obs){TFile*f=TFile::Open(Form("/data/uboone/processed/closure_hists_xsec_%s_%s.root",c,o));if(!f||f->IsZombie()){printf("  %s %s NO FILE\n",c,o);continue;}TH1*u=(TH1*)f->Get("h_unfolded_nuwro");TH1*t=(TH1*)f->Get("h_fakedata_truth");printf("  %-8s %-8s nbins=%d unfolded=%.3f truth=%.3f\n",c,o,u->GetNbinsX(),u->Integral("width"),t->Integral("width"));f->Close();}}' 2>/dev/null | grep -E "RHCFULL|COMB"
for c in RHCFULL COMB; do for o in costhpi thmupi; do
  lg=$([ "$c" = RHCFULL ] && echo "../logs/rhcfull_${o}_unfold_coarse.log" || echo "../logs/comb_${o}_unfold_coarse.log")
  grep -iE "^truth: " "$lg" | head -1 | sed "s|^|  $c $o: |"
done; done
echo "########## EXTEND RHC+COMB DONE ##########"
