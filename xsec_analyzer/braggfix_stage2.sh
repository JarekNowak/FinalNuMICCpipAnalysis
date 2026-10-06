#!/usr/bin/env bash
# braggfix_stage2.sh -- after the staging reprocess (slurm/slurm_braggfix_reprocess.sbatch, all units DONE):
# throw the fake data of both trees from the staging overlays with the release seeds (identical events,
# since entries and CV weights are unchanged), then run the staging-vs-live gate (slurm/braggfix_gate.py).
#   ./braggfix_stage2.sh > ../logs/braggfix/stage2.log 2>&1        (SKIP_THROWS=1: alignment and gate only)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
PROC=/data/uboone/processed; L=../logs/braggfix; mkdir -p "$L"
n=$(grep -l "^DONE .* rc=0 .*write_errors=0" slurm/braggfix_*_*.out 2>/dev/null | wc -l)
[ "$n" -eq "$(wc -l < slurm/braggfix_manifest.list)" ] || { echo "only $n units done"; exit 1; }
if [ "${SKIP_THROWS:-0}" != 1 ]; then
echo "== throws $(date)"
THROW_DIR=$PROC/incl_braggfix/ root.exe -l -b -q 'macros/throw_perrun_fhc.C(1)' > $L/throw_incl_fhc.log 2>&1 &
THROW_DIR=$PROC/incl_braggfix/ root.exe -l -b -q 'macros/throw_perrun_rhc.C(1)' > $L/throw_incl_rhc.log 2>&1 &
THROW_W_DIR=$PROC/w_braggfix/ root.exe -l -b -q 'macros/throw_perrun_w.C(1)' > $L/throw_w_fhc.log 2>&1 &
# loaded with .L: running the macro file would also execute its default FHC throw, writing the FHC files
# a second time (and concurrently with the line above)
THROW_W_DIR=$PROC/w_braggfix/ root.exe -l -b -q -e '.L macros/throw_perrun_w.C' -e 'throw_perrun_w_rhc(1)' > $L/throw_w_rhc.log 2>&1 &
wait
grep -h "POTSCALE" $L/throw_*.log
[ "$(grep -h POTSCALE $L/throw_*.log | wc -l)" -eq 16 ] || { echo "throws incomplete"; exit 1; }
fi
echo "== summed_pot alignment $(date)"
# precedent of reprocess_thpipr.sh (2026-09-17): a native summed_pot that differs from the live one only
# by float rounding (relative < 1e-5; the FHC Run-1 dirt snapshot, 2e-6) is set to the live value so the
# normalisation stays bit-identical; anything larger is left for the gate to fail
python3 - <<'PY'
import os, uproot, ROOT
P = '/data/uboone/processed'
for st, lv in ((P + '/incl_braggfix', P), (P + '/w_braggfix', P + '/w'), (P + '/ext_perrun_braggfix', P + '/ext_perrun')):
    for b in sorted(os.listdir(st)):
        f = os.path.join(st, b)
        if not b.endswith('.root') or os.path.islink(f) or b.startswith('xsec-ana-fakedata'): continue
        n = float(uproot.open(f)['summed_pot'].member('fVal')); o = float(uproot.open(os.path.join(lv, b))['summed_pot'].member('fVal'))
        if n != o and abs(n - o) < 1e-5 * abs(o):
            t = ROOT.TFile.Open(f, 'UPDATE'); ROOT.TParameter('float')('summed_pot', o).Write('summed_pot', ROOT.TObject.kOverwrite); t.Close()
            print(f'summed_pot of {f} set to the live value {o:.8g} (native {n:.8g}, rel. diff {abs(n - o) / o:.1e})')
PY
echo "== gate $(date)"
python3 slurm/braggfix_gate.py
