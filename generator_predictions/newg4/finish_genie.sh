source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4
# sigma_tot in 1e-38 units from genie_sigmatot.txt -> convert to cm^2 (x1e-38)
SNU=$(awk '/flux_numu /{printf "%.6e", $2*1e-38}' genie_sigmatot.txt)
SNB=$(awk '/flux_numubar/{printf "%.6e", $2*1e-38}' genie_sigmatot.txt)
echo "sigma_tot(cm2): numu=$SNU numubar=$SNB"
root -l -b -q "../genie/analyze_gst.C(\"g_numu.gst.root\",${SNU},\"pred_g_numu.root\")"       2>&1 | grep -i "gst:"
root -l -b -q "../genie/analyze_gst.C(\"g_numubar.gst.root\",${SNB},\"pred_g_numubar.root\")" 2>&1 | grep -i "gst:"
root -l -b -q "combine_newg4.C(\"pred_g_numu.root\",\"pred_g_numubar.root\",\"genie_newg4_final.root\")" 2>&1 | grep -i integral
