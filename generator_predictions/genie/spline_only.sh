set +e
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh >/dev/null 2>&1
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup >/dev/null 2>&1
setup genie v3_04_00d -q c14:prof >/dev/null 2>&1
setup genie_xsec v3_04_00 -q G1810a0211a:k250:e1000 >/dev/null 2>&1
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/genie
echo "start gspl2root $(free -g | awk 'NR==2{print \"mem free \" $4 \"G\"}')"
gspl2root -f "$GENIEXSECFILE" -p 14 -t 1000180400 -o splines_numu_ar.root
echo "gspl2root exit=$? ; file: $(ls -la splines_numu_ar.root 2>/dev/null | awk '{print $5}')"
