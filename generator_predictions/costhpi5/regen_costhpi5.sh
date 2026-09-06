#!/usr/bin/env bash
# regen_costhpi5.sh -- regenerate the four generator predictions with the forward-split cos(theta_pi)
# binning (adopted 2026-09-06) for FHC, RHC and combined; verify the unchanged observables against the
# live FTE files before installing.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/../xsec_analyzer/lib:${LD_LIBRARY_PATH:-}"
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
echo "==== COSTHPI5 REGEN START $(date) ===="
nice root.exe -l -b -q 'costhpi5/analyze_gst_c5.C("newg4/g_numu.gst.root",1.390129e-37,"costhpi5/genie_c5_numu.root")' 2>&1 | tail -2
nice root.exe -l -b -q 'costhpi5/analyze_gst_c5.C("newg4/g_numubar.gst.root",4.208102e-38,"costhpi5/genie_c5_numubar.root")' 2>&1 | tail -2
echo "GENIE done $(date)"
nice root.exe -l -b -q 'costhpi5/gibuu_cc1pi_c5.C("newg4/gibuu/run_numu/FinalEvents.dat",20,"costhpi5/gibuu_c5_numu.root")' 2>&1 | tail -2
nice root.exe -l -b -q 'costhpi5/gibuu_cc1pi_c5.C("newg4/gibuu/run_numubar/FinalEvents.dat",20,"costhpi5/gibuu_c5_numubar.root")' 2>&1 | tail -2
echo "GiBUU done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash costhpi5/reana_nuwro_c5.sh 2>&1 | tail -4
echo "NuWro done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash costhpi5/reana_neut_c5.sh 2>&1 | tail -4
echo "NEUT done $(date)"
for g in genie gibuu nuwro neut; do
  nice root.exe -l -b -q "newg4/combine_newg4.C(\"costhpi5/${g}_c5_numu.root\",\"costhpi5/${g}_c5_numubar.root\",\"costhpi5/${g}_c5_final.root\")" 2>&1 | grep integral
  nice root.exe -l -b -q "newg4/combine_rhc.C(\"costhpi5/${g}_c5_numu.root\",\"costhpi5/${g}_c5_numubar.root\",\"costhpi5/${g}_c5_rhc_final.root\")" 2>&1 | grep integral
  nice root.exe -l -b -q "newg4/make_fte.C(\"costhpi5/${g}_c5_final.root\",\"costhpi5/${g}_c5_fte.root\")" 2>&1 | grep -i 'warn\|theta' | head -2
  nice root.exe -l -b -q "newg4/make_fte.C(\"costhpi5/${g}_c5_rhc_final.root\",\"costhpi5/${g}_c5_rhc_fte.root\")" 2>&1 | grep -i 'warn' | head -2
  nice root.exe -l -b -q "newg4/combine_comb_fte.C(\"costhpi5/${g}_c5_fte.root\",\"costhpi5/${g}_c5_rhc_fte.root\",\"costhpi5/${g}_c5_comb_fte.root\")" 2>&1 | tail -1
done
echo "==== COSTHPI5 REGEN DONE $(date) ===="
