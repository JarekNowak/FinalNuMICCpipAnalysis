#!/usr/bin/env bash
# braggfix_mcs.sh -- second pass after the release without the Bragg-pion cut (version 1.7): the data-side MCS
# momentum-scale term from the new pseudo-data. The +-5% pseudo-data (x1.05 "up05", x0.95 "dn05" on the
# MCS-measured muons) are regenerated from the promoted nominal fake data of both trees, after the s=1 closure
# check of the proton-tagged recomputation; the previous files and extractions move to pre_braggfix
# directories; then the 38 extractions of slurm/mcsdata_{050,1p,2d}_manifest.list are submitted
# (slurm_extract.sbatch). Afterwards: report/tools/mcs_eval.py all incl / all 1p and the table chain.
#   ./braggfix_mcs.sh > ../logs/braggfix/mcs.log 2>&1
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
PROC=/data/uboone/processed; RA=$PROC/rebuild_alt; BK=$PROC/pre_braggfix/mcs_fakedata; RBK=$RA/pre_braggfix
[ -e "$BK" ] || [ -e "$RBK" ] && { echo "REFUSING: backup directories exist"; exit 1; }
echo "== proton-tagged s=1 closure check $(date)"
CHK=$( { root -l -b -e 'gSystem->AddIncludePath("-I'"$PWD"'/include")' -q 'macros/mcs_scale_fakedata_1p.C+(1.00,"chk")'
         root -l -b -e 'gSystem->AddIncludePath("-I'"$PWD"'/include")' -q 'macros/mcs_scale_fakedata_1p.C+(1.00,"chk","rhc")'; } 2>&1 | grep "deviate" )
echo "$CHK"
[ "$(echo "$CHK" | grep -c ', 0 deviate')" -eq 8 ] || { echo "REFUSING: the s=1 recomputation does not reproduce the stored branches"; exit 1; }
echo "== backups $(date)"
mkdir -p "$BK/incl" "$BK/w" "$RBK"
mv $PROC/xsec-ana-fakedata_*_mcs{up05,dn05}.root "$BK/incl/"
mv $PROC/w/xsec-ana-fakedata_*_mcs{up05,dn05}.root "$BK/w/"
for M in slurm/mcsdata_050_manifest.list slurm/mcsdata_1p_manifest.list slurm/mcsdata_2d_manifest.list; do
  awk '{print $1}' "$M" | while read -r xc; do
    U=$(awk '$1=="UnivFile"{print $2}' "$xc"); B=$(basename "$U" _univmake.root)
    for f in "$U" "$RA/xsec_$B.root" "$RA/closure_hists_xsec_$B.root" "$RA/univ_$B.log" "$RA/unfold_$B.log" "$RA/plots_$B"; do
      [ -e "$f" ] && mv "$f" "$RBK/"; done
  done
done
echo "== +-5% pseudo-data $(date)"
for b in fhc rhc; do
  for v in "1.05,\"up05\"" "0.95,\"dn05\""; do
    root -l -b -q "macros/mcs_scale_fakedata.C($v,\"$b\")" 2>&1 | grep -v "^Processing" | tail -5
    root -l -b -e 'gSystem->AddIncludePath("-I'"$PWD"'/include")' -q "macros/mcs_scale_fakedata_1p.C+($v,\"$b\")" 2>&1 | grep -v "^Processing" | tail -5
  done
done
ls $PROC/xsec-ana-fakedata_*_mcs*.root $PROC/w/xsec-ana-fakedata_*_mcs*.root | wc -l
echo "== extractions $(date)"
for M in mcsdata_050 mcsdata_1p mcsdata_2d; do
  N=$(grep -c '' slurm/${M}_manifest.list)
  J=$(sbatch --parsable --chdir=$PWD/slurm --cpus-per-task=4 --mem=12G --array=0-$((N-1)) --export=ALL,REPO=$PWD,MANIFEST=$PWD/slurm/${M}_manifest.list,UNIV_THREADS=4 slurm/slurm_extract.sbatch)
  echo "$J ${M}_manifest.list ($N tasks; data-side MCS scale +-5%, release without the Bragg-pion cut) $(date)" | tee -a slurm/rebuild_jobs_2026-09-19.txt
done
