#!/usr/bin/env bash
# regen_ext.sh -- generator predictions for the 2026-09-16 observables (theta_p, theta_pipr, 8 2D
# pairs), with FHC, RHC and COMB FTE files for the inclusive AND the proton-tagged sample.
# Readers are gen2d copies of the released ones (obs_ext.h fills added); outputs stay in gen2d/.
# Verification: every released histogram the copies also write must be bit-identical to the
# released per-flavour files, and the combiner is checked against the released FTE files.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/../xsec_analyzer/lib:${LD_LIBRARY_PATH:-}"
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
G=gen2d
echo "==== GEN EXT START $(date) ===="
nice root.exe -l -b -q "$G/analyze_gst_incl.C(\"newg4/g_numu.gst.root\",1.390129e-37,\"$G/genie_ext_numu.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/analyze_gst_incl.C(\"newg4/g_numubar.gst.root\",4.208102e-38,\"$G/genie_ext_numubar.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/analyze_gst_1p.C(\"newg4/g_numu.gst.root\",1.390129e-37,\"$G/genie_1p_ext_numu.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/analyze_gst_1p.C(\"newg4/g_numubar.gst.root\",4.208102e-38,\"$G/genie_1p_ext_numubar.root\")" 2>&1 | tail -1
echo "GENIE done $(date)"
nice root.exe -l -b -q "$G/gibuu_incl.C(\"newg4/gibuu/run_numu/FinalEvents.dat\",20,\"$G/gibuu_ext_numu.root\")" 2>&1 | tail -2
nice root.exe -l -b -q "$G/gibuu_incl.C(\"newg4/gibuu/run_numubar/FinalEvents.dat\",20,\"$G/gibuu_ext_numubar.root\")" 2>&1 | tail -2
nice root.exe -l -b -q "$G/gibuu_1p.C(\"newg4/gibuu/run_numu/FinalEvents.dat\",20,\"$G/gibuu_1p_ext_numu.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/gibuu_1p.C(\"newg4/gibuu/run_numubar/FinalEvents.dat\",20,\"$G/gibuu_1p_ext_numubar.root\")" 2>&1 | tail -1
echo "GiBUU done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash $G/reana_nuwro_ext.sh 2>&1 | tail -6
echo "NuWro done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash $G/reana_neut_ext.sh 2>&1 | tail -8
echo "NEUT done $(date)"

INC=pmu,ppi,costhmu,costhpi,thmupi,costhpi_costhmu,costhmu_pmu,costhpi_thmupi,costhpi_ppi
ONEP=Wpipr,Whad,dpt,dalphat,dphit,pn,thetap,thpipr,thetap_dpt,thpipr_dpt,thetap_Wpipr,thpipr_Wpipr
root.exe -l -b -q -e '.L gen2d/combine_ext.C+' >/dev/null 2>&1
for g in genie gibuu nuwro neut; do
  for s in ext:$INC 1p_ext:$ONEP; do
    t=${s%%:*}; o=${s#*:}
    root.exe -l -b -q -e ".L $G/combine_ext.C+" -e "combine_ext(\"$G/${g}_${t}_numu.root\",\"$G/${g}_${t}_numubar.root\",\"$G/${g}_${t}_fhc_fte.root\",\"fhc\",\"$o\"); combine_ext(\"$G/${g}_${t}_numu.root\",\"$G/${g}_${t}_numubar.root\",\"$G/${g}_${t}_rhc_fte.root\",\"rhc\",\"$o\"); combine_comb_ext(\"$G/${g}_${t}_fhc_fte.root\",\"$G/${g}_${t}_rhc_fte.root\",\"$G/${g}_${t}_comb_fte.root\",\"$o\")" 2>&1 | grep -E "missing|sigma_int" | grep -E "missing|_|fhc" | head -30
  done
done
echo "---- verification: released histograms reproduced by the gen2d readers ----"
root.exe -l -b -q -e '.L gen2d/verify_ext.C+' >/dev/null 2>&1
for g in genie gibuu nuwro neut; do
  for fl in numu numubar; do
    root.exe -l -b -q "$G/verify_ext.C+(\"$G/${g}_ext_${fl}.root\",\"costhpi5/${g}_c5_${fl}.root\",\"pmu,ppi,costhmu,costhpi,thmupi\")" 2>&1 | grep VERIFY | sed "s/^/$g $fl incl /"
  done
done
for g in genie gibuu nuwro neut; do
  root.exe -l -b -q "$G/verify_ext.C+(\"$G/${g}_1p_ext_fhc_fte.root\",\"newg4/${g}_wtki_2bin_fte.root\",\"dpt_fte,dalphat_fte,dphit_fte,pn_fte,Wpipr_fte,Whad_fte\")" 2>&1 | grep VERIFY | sed "s/^/$g 1p-fhc vs wtki_2bin /"
done
echo "==== GEN EXT DONE $(date) ===="
