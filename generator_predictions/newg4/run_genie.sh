HERE=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
cd "$HERE/newg4"; echo "[genie-newg4] start"
set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000
setup genie_phyopt v3_04_00 -q dkcharmtau
[[ -z "$GENIEXSECFILE" ]] && { echo fatal; exit 1; }
set -e
FLUX="$HERE/newg4/flux_newg4_clean.root"
# generate numu + numubar on Ar-40 with cleaned newg4 flux, E>60 MeV
[[ -s g_numu ]]    || gevgen -n 200000 -p  14 -t 1000180400 -e 0.06,10 -f "$FLUX,flux_numu"    --cross-sections "$GENIEXSECFILE" --tune G18_10a_02_11a --seed 1 -o g_numu
[[ -s g_numubar ]] || gevgen -n 200000 -p -14 -t 1000180400 -e 0.06,10 -f "$FLUX,flux_numubar" --cross-sections "$GENIEXSECFILE" --tune G18_10a_02_11a --seed 2 -o g_numubar
gntpc -i g_numu    -f gst -o g_numu.gst.root
gntpc -i g_numubar -f gst -o g_numubar.gst.root
# reuse the existing numu spline (nu_mu_Ar40 tot_cc); need numubar spline too
[[ -s ../genie/splines_numubar_ar.root ]] || gspl2root -f "$GENIEXSECFILE" -p -14 -t 1000180400 -o ../genie/splines_numubar_ar.root
# flux-avg tot_cc over the cleaned newg4 flux, per species
root -l -b -q "$HERE/newg4/genie_fluxavg.C" 2>&1 | grep -iE "sigma|numu"
echo "[genie-newg4] gen+convolve done"
