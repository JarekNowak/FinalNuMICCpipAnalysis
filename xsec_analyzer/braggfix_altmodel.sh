#!/usr/bin/env bash
# braggfix_altmodel.sh -- second pass after the release without the Bragg-pion cut (version 1.7): the
# alternative-model closure. The reweighted pseudo-data (GENIE multisim universe 545 "altgenie", Delta->N pi
# angular unisim "altdelta") are re-thrown from the promoted trees with the recorded arguments and seeds
# (inclusive FHC/RHC, proton-tagged FHC/RHC); the event counts must equal the previous throws (weights and
# seeds unchanged). The previous files and extractions move to pre_braggfix directories, then the 38
# extractions of slurm/braggfix_altmodel_manifest.list are submitted (slurm_extract.sbatch).
# Afterwards: report/tools/altmodel_eval.py {FHC5,RHCFULL} incl, FHC5 1p, V1.
#   ./braggfix_altmodel.sh > ../logs/braggfix/altmodel.log 2>&1
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
PROC=/data/uboone/processed; RA=$PROC/rebuild_alt; BK=$PROC/pre_braggfix/alt_fakedata; RBK=$RA/pre_braggfix/alt; L=../logs/braggfix
[ -e "$BK" ] || [ -e "$RBK" ] && { echo "REFUSING: backup directories exist"; exit 1; }
mkdir -p "$BK/incl" "$BK/w" "$RBK"
mv $PROC/xsec-ana-fakedata_alt{genie,delta}_{fhc,rhc}_run[0-9].root "$BK/incl/"; mv $PROC/w/xsec-ana-fakedata_alt{genie,delta}_{fhc,rhc}_run[0-9].root "$BK/w/"   # the *_pot3283 throws stay
echo "== throws $(date)"
G='"weight_All_UBGenie",545,"altgenie",1'; D='"weight_Theta_Delta2Npi_UBGenie",0,"altdelta",1'
for a in "$G" "$D"; do
  t=$(echo "$a" | cut -d, -f3 | tr -d '"')
  root -l -b -q "macros/throw_altmodel_fhc.C($a)" > $L/alt_throw_incl_fhc_$t.log 2>&1 &
  root -l -b -q "macros/throw_altmodel_rhc.C($a)" > $L/alt_throw_incl_rhc_$t.log 2>&1 &
  root -l -b -q "macros/throw_altmodel_w.C($a)" > $L/alt_throw_w_fhc_$t.log 2>&1 &
  THROW_DIR=$PROC/w/ root -l -b -q "macros/throw_altmodel_rhc.C($a)" > $L/alt_throw_w_rhc_$t.log 2>&1 &
done
wait
python3 - <<'PY'
import glob, os, uproot, sys
P = '/data/uboone/processed/'; B = P + 'pre_braggfix/alt_fakedata/'; bad = []; n = 0
for d, b in (('', 'incl/'), ('w/', 'w/')):
    for old in sorted(glob.glob(B + b + 'xsec-ana-fakedata_alt*.root')):
        new = P + d + os.path.basename(old); n += 1
        if not os.path.exists(new): bad.append(('missing', new)); continue
        a, o = uproot.open(new)['stv_tree'].num_entries, uproot.open(old)['stv_tree'].num_entries
        if a != o: bad.append((os.path.basename(new), a, o))
print(f'{n} re-thrown files, event counts equal to the previous throws: {not bad}', bad)
sys.exit(1 if bad or n != 32 else 0)
PY
echo "== backups of the previous extractions $(date)"
while read -r xc b s; do
  U=$(awk '$1=="UnivFile"{print $2}' "$xc"); X=$(basename "$U" _univmake.root)
  for f in "$U" "$RA/xsec_$X.root" "$RA/closure_hists_xsec_$X.root" "$RA/closure_hists_all_xsec_$X.root" "$RA/univ_$X.log" "$RA/unfold_$X.log" "$RA/xsec_$X.log" "$RA/plots_$X"; do
    [ -e "$f" ] && mv "$f" "$RBK/"; done
done < slurm/braggfix_altmodel_manifest.list
echo "== extractions $(date)"
N=$(grep -c '' slurm/braggfix_altmodel_manifest.list)
J=$(sbatch --parsable --chdir=$PWD/slurm --cpus-per-task=4 --mem=12G --array=0-$((N-1)) --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/braggfix_altmodel_manifest.list,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
echo "$J braggfix_altmodel_manifest.list ($N tasks; alternative-model closure, release without the Bragg-pion cut) $(date)" | tee -a slurm/rebuild_jobs_2026-09-19.txt
