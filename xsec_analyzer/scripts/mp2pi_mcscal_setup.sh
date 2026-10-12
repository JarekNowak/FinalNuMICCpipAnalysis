#!/bin/bash
# mp2pi_mcscal_setup.sh -- framework inputs of the two-pion one-bin total with the MCS-calibrated classifier inputs,
# after the production into /data/uboone/processed/mp2pi_mcscal (slurm/slurm_mp2pi.sbatch for the overlays and the
# dirt, slurm/slurm_braggfix_reprocess.sbatch with slurm/mp2pi_mcscal_detvar_manifest.list for the detector variations,
# both with CC1MU1PIXP_MP_MCS_MAP; beam-off outputs linked from /data/uboone/processed/mp2pi):
#   1. file lists configs/file_properties_numi_<cfg>_mp2pi_mcscal.txt, links and the pseudo-data plan (scripts/mp2pi_setup.py)
#   2. the plan must equal that of the uncalibrated production apart from the directory (same files, scales, seeds)
#   3. the per-run pseudo-data (macros/throw_mp2pi.C); same seeds and weights, so the same events as the uncalibrated
#      pseudo-data (entry counts compared)
#   4. lists with the flux-weighted detector variations, configs/file_properties_numi_<cfg>_mp2pi_mcscal_fluxw.txt
#      (detector variations in detvar_fluxw/, scripts/detvar_fluxw.py with MP2PI_DIR)
#   sbatch --dependency=afterok:<production> --wrap="bash scripts/mp2pi_mcscal_setup.sh"
set -euo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$R"
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$R/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$R"
D=/data/uboone/processed/mp2pi_mcscal; U=/data/uboone/processed/mp2pi

echo "==== 1. file lists and plan $(date)"
MP2PI_DIR=$D MP2PI_SUFFIX=_mcscal python3 scripts/mp2pi_setup.py

echo "==== 2. plan against the uncalibrated one $(date)"
diff <(sed "s#$D/#DIR/#g" $D/throw_plan.txt) <(sed "s#$U/#DIR/#g" $U/throw_plan.txt) && echo "plan identical apart from the directory"

echo "==== 3. pseudo-data $(date)"
root -l -b -q "macros/throw_mp2pi.C(\"$D/throw_plan.txt\")"
# (entry counts with uproot: root -e returns the value of its last expression as the exit code)
python3 - "$D" "$U" <<'PYEOF'
import sys, uproot
d, u = sys.argv[1], sys.argv[2]
outs = [l.split()[1] for l in open(d + '/throw_plan.txt') if l.startswith('OUT ')]
bad = 0
for f in outs:
    a, b = (uproot.open(f'{x}/{f}')['stv_tree'].num_entries for x in (d, u))
    print(f'  {f}: {a} entries' + (', as uncalibrated' if a == b else f' against {b} uncalibrated')); bad |= a != b
sys.exit('pseudo-data differ from the uncalibrated throw' if bad else 0)
PYEOF

echo "==== 4. lists with the flux-weighted detector variations $(date)"
for cfg in fhc5 rhcfull comb; do
  src=configs/file_properties_numi_${cfg}_mp2pi_mcscal.txt; dst=configs/file_properties_numi_${cfg}_mp2pi_mcscal_fluxw.txt
  { echo "# two-pion one-bin total ($cfg), MCS-calibrated classifier inputs, with the detector-variation samples weighted to the new NuMI flux (scripts/detvar_fluxw.py): $(basename $src) with the detVar paths in detvar_fluxw/"
    tail -n +2 $src | sed "s#$D/xsec-ana-detvar_#$D/detvar_fluxw/xsec-ana-detvar_#"; } > $dst
  echo "  $dst: $(grep -c "$D/detvar_fluxw/" $dst) flux-weighted detector-variation entries"
done
echo "==== SETUP DONE $(date)"
