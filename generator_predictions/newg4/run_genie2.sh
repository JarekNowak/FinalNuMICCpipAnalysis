HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
cd "$HERE/newg4"; echo "[genie2] start"
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000
setup genie_phyopt v3_04_00 -q dkcharmtau
FLUX="$HERE/newg4/flux_newg4_clean.root"
# gevgen segfaults on shutdown AFTER writing all events, so ignore its exit code
# and check the output file instead.
[[ -s g_numu    ]] || gevgen -n 200000 -p  14 -t 1000180400 -e 0.06,10 -f "$FLUX,flux_numu"    --cross-sections "$GENIEXSECFILE" --tune G18_10a_02_11a --seed 1 -o g_numu    || true
[[ -s g_numubar ]] || gevgen -n 200000 -p -14 -t 1000180400 -e 0.06,10 -f "$FLUX,flux_numubar" --cross-sections "$GENIEXSECFILE" --tune G18_10a_02_11a --seed 2 -o g_numubar || true
[[ -s g_numu ]] && [[ -s g_numubar ]] || { echo "[fatal] missing ghep"; exit 1; }
[[ -s g_numu.gst.root ]]    || gntpc -i g_numu    -f gst -o g_numu.gst.root
[[ -s g_numubar.gst.root ]] || gntpc -i g_numubar -f gst -o g_numubar.gst.root
[[ -s ../genie/splines_numubar_ar.root ]] || gspl2root -f "$GENIEXSECFILE" -p -14 -t 1000180400 -o ../genie/splines_numubar_ar.root
# flux-avg tot_cc per Ar over newg4 flux, per species
root -l -b -q genie_fluxavg.C 2>&1 | grep -iE "flux-avg|numu"
SNU=$(awk '/flux_numu /{print $2}' genie_sigmatot.txt)
SNB=$(awk '/flux_numubar/{print $2}' genie_sigmatot.txt)
echo "sigma_tot numu=$SNU numubar=$SNB [1e-38 per Ar]"
root -l -b -q "../genie/analyze_gst.C(\"g_numu.gst.root\",${SNU}e-38,\"pred_g_numu.root\")"    2>&1 | grep CC1pi
root -l -b -q "../genie/analyze_gst.C(\"g_numubar.gst.root\",${SNB}e-38,\"pred_g_numubar.root\")" 2>&1 | grep CC1pi
root -l -b -q "combine_newg4.C(\"pred_g_numu.root\",\"pred_g_numubar.root\",\"genie_newg4_final.root\")" 2>&1 | grep integral
echo "[genie2] done"
