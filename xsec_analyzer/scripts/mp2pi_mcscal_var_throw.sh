#!/bin/bash
# mp2pi_mcscal_var_throw.sh <variant> -- pseudo-data of the two-pion one-bin total with the classifier inputs
# calibrated by a variant MCS map (shrink, coarse, fine, old_production, new_production, FHC, RHC of
# scripts/mcs_cal_derive.py), for the calibration uncertainty. The overlays and the dirt processed with that map
# (slurm/slurm_mp2pi.sbatch on slurm/mp2pi_mc_manifest.list, OUT=/data/uboone/processed/mp2pi_mcscal/var/<variant>)
# are thrown with the plan of the nominal calibration (same scales and seeds, so the same events) into the same
# directory; the entry counts must equal those of the nominal pseudo-data. Then the processed samples of the variant
# are deleted (only the pseudo-data are used).
#   sbatch --dependency=afterok:<variant production> --wrap="bash scripts/mp2pi_mcscal_var_throw.sh <variant>"
# MP2PI_DIR: another nominal production (e.g. /data/uboone/processed/mp2pi_2f, scripts/mp2pi_2f_mcsvar_submit.sh)
set -euo pipefail
V="${1:?variant}"
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$R"
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$R/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$R"
N="${MP2PI_DIR:-/data/uboone/processed/mp2pi_mcscal}"; D=$N/var/$V
echo "==== variant $V: plan $(date)"
awk -v n="$N/" -v d="$D/" '$1=="IN"{ sub(n, d, $2) } { print }' $N/throw_plan.txt > $D/throw_plan.txt
for f in $(awk '$1=="IN"{print $2}' $D/throw_plan.txt); do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done
diff <(sed "s#$D/#DIR/#g" $D/throw_plan.txt) <(sed "s#$N/#DIR/#g" $N/throw_plan.txt) && echo "plan identical apart from the directory"
grep -q "MCS map for the classifier inputs: /data/uboone/processed/mcscal/maps/$V.txt" $D/logs/*.out \
  || { echo "the variant production does not record the map $V"; exit 1; }
echo "==== throw $(date)"
root -l -b -q "macros/throw_mp2pi.C(\"$D/throw_plan.txt\")" 2>&1 | grep -v "^Error in <TBranchElement::Fill>"
python3 - "$D" "$N" <<'PYEOF'
import sys, uproot
d, n = sys.argv[1], sys.argv[2]
outs = [l.split()[1] for l in open(d + '/throw_plan.txt') if l.startswith('OUT ')]
bad = 0
for f in outs:
    a, b = (uproot.open(f'{x}/{f}')['stv_tree'].num_entries for x in (d, n))
    print(f'  {f}: {a} entries' + (', as nominal' if a == b else f' against {b} nominal')); bad |= a != b
sys.exit('pseudo-data differ from the nominal throw' if bad else 0)
PYEOF
echo "==== delete the processed samples of the variant $(date)"
for f in $(awk '$1=="IN"{print $2}' $D/throw_plan.txt | sort -u); do command rm -f -- "$f"; done
echo "==== THROW DONE $V $(date)"
