#!/usr/bin/env bash
# regen_1pinc.sh (2026-09-27) -- generator predictions for the muon/pion observables measured on the PROTON-TAGGED
# sample (p_mu, p_pi two regions, cos theta_mu, cos theta_pi four bins, theta_mu_pi) at the proton-tagged binnings,
# FHC/RHC/COMB, all four generators. Readers: gen2d/*_1p_inc.* = the 1p readers + the obs_1pinc.h fills; outputs are
# NEW files gen2d/<gen>_1p_inc_{numu,numubar,fhc_fte,rhc_fte,comb_fte}.root (nothing existing is overwritten).
# Verification: (1) every histogram shared with the released gen2d/<gen>_1p_ext_<flavour>.root is bit-identical
# (same selection and normalisation); (2) edges vs the bin configs; (3) each new observable integrates to the
# proton-tagged total of the same generator and configuration.
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/../xsec_analyzer/lib:${LD_LIBRARY_PATH:-}"
IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
G=gen2d
echo "==== GEN 1P-INC START $(date) ===="
nice root.exe -l -b -q "$G/analyze_gst_1p_inc.C(\"newg4/g_numu.gst.root\",1.390129e-37,\"$G/genie_1p_inc_numu.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/analyze_gst_1p_inc.C(\"newg4/g_numubar.gst.root\",4.208102e-38,\"$G/genie_1p_inc_numubar.root\")" 2>&1 | tail -1
echo "GENIE done $(date)"
nice root.exe -l -b -q "$G/gibuu_1p_inc.C(\"newg4/gibuu/run_numu/FinalEvents.dat\",20,\"$G/gibuu_1p_inc_numu.root\")" 2>&1 | tail -1
nice root.exe -l -b -q "$G/gibuu_1p_inc.C(\"newg4/gibuu/run_numubar/FinalEvents.dat\",20,\"$G/gibuu_1p_inc_numubar.root\")" 2>&1 | tail -1
echo "GiBUU done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash $G/reana_nuwro_1pinc.sh 2>&1 | tail -6
echo "NuWro done $(date)"
apptainer exec -B /cvmfs -B /home -B /data "$IMG" bash $G/reana_neut_1pinc.sh 2>&1 | tail -8
echo "NEUT done $(date)"

OBS=pmu1p,ppi1p,costhmu1p,costhpi1p,thmupi1p
for g in genie gibuu nuwro neut; do
  root.exe -l -b -q -e ".L $G/combine_ext.C+" -e "combine_ext(\"$G/${g}_1p_inc_numu.root\",\"$G/${g}_1p_inc_numubar.root\",\"$G/${g}_1p_inc_fhc_fte.root\",\"fhc\",\"$OBS\"); combine_ext(\"$G/${g}_1p_inc_numu.root\",\"$G/${g}_1p_inc_numubar.root\",\"$G/${g}_1p_inc_rhc_fte.root\",\"rhc\",\"$OBS\"); combine_comb_ext(\"$G/${g}_1p_inc_fhc_fte.root\",\"$G/${g}_1p_inc_rhc_fte.root\",\"$G/${g}_1p_inc_comb_fte.root\",\"$OBS\")" 2>&1 | grep -E "missing|sigma_int" | sed "s/^/$g /"
done
echo "---- verification 1: histograms shared with the released 1p_ext per-flavour files are bit-identical ----"
SHARED=Wpipr,Whad,dpt,dalphat,dphit,pn,thetap,thpipr,thetap_dpt,thpipr_dpt,thetap_Wpipr,thpipr_Wpipr,Whad2,dphit3,thetap4m,thpipr4m,pp2,pp3,pp4,pp5
for g in genie gibuu nuwro neut; do
  for fl in numu numubar; do
    root.exe -l -b -q "$G/verify_ext.C+(\"$G/${g}_1p_inc_${fl}.root\",\"$G/${g}_1p_ext_${fl}.root\",\"$SHARED\")" 2>&1 | grep VERIFY | awk -v t="$g $fl" '{n++; if($0~/SAME/)s++; else print t, $0} END{printf "%s %s: %d/%d shared histograms SAME\n", t, "", s, n}'
  done
done
echo "---- verification 2: edges against the proton-tagged bin configs ----"
python3 $G/check_1pinc_bins.py
echo "---- verification 3: each observable integrates to the proton-tagged total (dpt_fte of the 1p_ext files) ----"
python3 - <<'PYEOF'
import uproot
G='gen2d'; worst=0.
for g in ('genie','gibuu','neut','nuwro'):
  for m in ('fhc','rhc','comb'):
    tot=uproot.open(f'{G}/{g}_1p_ext_{m}_fte.root')['dpt_fte'].values().sum(); f=uproot.open(f'{G}/{g}_1p_inc_{m}_fte.root')
    v={o:f[o+'_fte'].values().sum() for o in ('pmu1p','ppi1p','costhmu1p','costhpi1p','thmupi1p')}
    d=max(abs(x-tot)/tot for x in v.values()); worst=max(worst,d)
    print(f'  {g:6s} {m:4s} total {tot:.5f}  ' + '  '.join(f'{o} {x:.5f}' for o,x in v.items()) + f'  max rel dev {d:.1e}')
print(f'integral check: worst relative deviation from the proton-tagged total {worst:.2e}')
PYEOF
echo "==== GEN 1P-INC DONE $(date) ===="
