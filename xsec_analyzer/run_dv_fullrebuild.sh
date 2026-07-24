#!/usr/bin/env bash
# run_dv_fullrebuild.sh — CORRECT detector-variation-as-fake-data study for
# costhmu: for each variation, rebuild the univmake with that variation in the
# onBNB (fake-data) slot, then unfold. Unlike the histogram-swap shortcut, this
# rebuilds the measured data (onBNB_reco) through the real data path, so the
# unfolded result genuinely responds to each variation.
set -uo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; cd "$HERE"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$PWD/lib:${LD_LIBRARY_PATH:-}"
PROC=/data/uboone/processed; LOG=../logs; mkdir -p "$LOG"
# variation : label : POT
VARS=( "LYDown:7.053740e20" "LYRayleigh:7.623732e20" "Recombination:7.537606e20" "SCE:7.532321e20"
       "WireModThetaXZ:7.628598e20" "WireModThetaYZ:7.653130e20" "WireModX:7.600698e20" "WireModYZ:7.509717e20" )

echo "########## SETUP: file_properties + configs ##########"
for vp in "${VARS[@]}"; do
  label="${vp%%:*}"; pot="${vp#*:}"
  trig=$(python3 -c "print(round(7809962*$pot/3.283e20))")
  onbnb="/data/uboone/processed/xsec-ana-genie_dv${label}_fakedata_run1_fhc.root 1 onBNB $trig $pot"
  # file_properties: replace the onBNB line
  awk -v repl="$onbnb" '/1 onBNB / && /_fakedata_run1_fhc.root/ {print repl; next} {print}' \
    configs/file_properties_numi.txt > configs/file_properties_numi_dv${label}.txt
  # xsec config: point UnivFile at this variation's univmake
  sed "s#ccpi_Run1_costhmu_univmake.root#ccpi_Run1_costhmu_dv${label}_univmake.root#" \
    configs/ccpi_xsec_config_numi_costhmu.txt > configs/ccpi_xsec_config_numi_costhmu_dv${label}.txt
done

echo "########## STAGE A: univmake (parallel) ##########"
pids=()
for vp in "${VARS[@]}"; do
  label="${vp%%:*}"
  FPM=configs/file_properties_numi_dv${label}.txt \
  BIN_CONFIG=configs/ccpi_costhmu_bin_config.txt \
  OUT=$PROC/ccpi_Run1_costhmu_dv${label}_univmake.root \
  ./run_universe_maker.sh > "$LOG/dvrebuild_${label}_univ.log" 2>&1 &
  pids+=($!)
  echo "  launched univmake[$label] pid=${pids[-1]}"
done
fail=0
i=0
for vp in "${VARS[@]}"; do
  label="${vp%%:*}"
  if wait "${pids[$i]}"; then echo "  univmake[$label] OK"; else echo "  univmake[$label] FAILED"; fail=1; fi
  i=$((i+1))
done

echo "########## STAGE B: unfold + extract totals ##########"
CH=$LOG/dvrebuild_totals.txt; : > "$CH"
mkdir -p unfold_output/dv_rebuild
for vp in "${VARS[@]}"; do
  label="${vp%%:*}"
  univ=$PROC/ccpi_Run1_costhmu_dv${label}_univmake.root
  if [[ ! -s "$univ" ]]; then echo "  [skip] $label: no univ"; continue; fi
  out=$(UnfolderNuMI configs/ccpi_xsec_config_numi_costhmu_dv${label}.txt configs/ccpi_costhmu_slice_config.txt \
          unfold_output/ccpi_Run1_costhmu_dv${label}_xsec.root 2>&1)
  chi2=$(echo "$out" | grep -m1 "truth:" | sed 's/^[[:space:]]*//')
  cp -f unfold_output/plot_costhetamu_0.pdf unfold_output/dv_rebuild/${label}_costhmu.pdf 2>/dev/null
  tots=$(root -l -b -q -e "TFile f(\"unfold_output/closure_hists_ccpi_Run1_costhmu_dv${label}_xsec.root\"); TH1D*u=(TH1D*)f.Get(\"h_unfolded_nuwro\"); TH1D*t=(TH1D*)f.Get(\"h_fakedata_truth\"); TH1D*a=(TH1D*)f.Get(\"h_gen_truthAC\"); printf(\"unf=%.4f trueRaw=%.4f trueAC=%.4f\n\", u->Integral(\"width\"), t->Integral(\"width\"), a?a->Integral(\"width\"):-1);" 2>/dev/null | grep -m1 "unf=")
  echo "$label  $tots  |  $chi2" | tee -a "$CH"
done
echo "########## DV-REBUILD DONE (univmake fail=$fail) ##########"
