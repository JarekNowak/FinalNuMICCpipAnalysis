#!/usr/bin/env bash
# gen_wtki_regen.sh — regenerate the four generator proton-tagged W/TKI predictions at
# the new 3-bin binning (wtki_gen.h updated). GENIE + GiBUU run in the native ROOT;
# NuWro + NEUT reprocess their event vectors inside the SL7 container. Then combine
# numu+numubar over the FHC flux into the per-Ar FTE overlays.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/../xsec_analyzer/lib:${LD_LIBRARY_PATH:-}"
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
echo "==== GEN W/TKI REGEN (3-bin) START $(date) ===="

# --- GENIE (gst; absolute sigma_tot_cc per Ar: numu 1.390129e-37, numubar 4.208102e-38) ---
nice root.exe -l -b -q 'genie/analyze_gst_1p.C("newg4/g_numu.gst.root",1.390129e-37,"newg4/genie_1p_numu.root")' >/dev/null 2>&1
nice root.exe -l -b -q 'genie/analyze_gst_1p.C("newg4/g_numubar.gst.root",4.208102e-38,"newg4/genie_1p_numubar.root")' >/dev/null 2>&1
nice root.exe -l -b -q 'combine_wtki.C("newg4/genie_1p_numu.root","newg4/genie_1p_numubar.root","newg4/genie_wtki_fte.root")' 2>&1 | grep Wpipr
echo "GENIE done"

# --- GiBUU (FinalEvents text, 20 runs) ---
( cd newg4/gibuu
  nice root.exe -l -b -q 'gibuu_cc1pi_1p.C("run_numu/FinalEvents.dat",20,"gibuu_1p_numu.root")' >/dev/null 2>&1
  nice root.exe -l -b -q 'gibuu_cc1pi_1p.C("run_numubar/FinalEvents.dat",20,"gibuu_1p_numubar.root")' >/dev/null 2>&1 )
nice root.exe -l -b -q 'combine_wtki.C("newg4/gibuu/gibuu_1p_numu.root","newg4/gibuu/gibuu_1p_numubar.root","newg4/gibuu_wtki_fte.root")' 2>&1 | grep Wpipr
echo "GiBUU done"

# --- NuWro (container: recompiles nuwro_cc1pi1p against event1.so) ---
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash newg4/reana_nuwro_1p.sh >/dev/null 2>&1
nice root.exe -l -b -q 'combine_wtki.C("newg4/nuwro_1p_numu.root","newg4/nuwro_1p_numubar.root","newg4/nuwro_wtki_fte.root")' 2>&1 | grep Wpipr
echo "NuWro done"

# --- NEUT (container: ACLiC neut_cc1pi1p) ---
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash newg4/neut/reana_neut_1p.sh >/dev/null 2>&1
nice root.exe -l -b -q 'combine_wtki.C("newg4/neut/neut_1p_numu.root","newg4/neut/neut_1p_numubar.root","newg4/neut_wtki_fte.root")' 2>&1 | grep Wpipr
echo "NEUT done"
echo "==== GEN W/TKI REGEN DONE $(date) ===="
