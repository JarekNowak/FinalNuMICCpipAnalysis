#!/usr/bin/env bash
# reprocess_thpipr.sh -- rebuild the proton-tagged tree with the pion-proton opening angle
# (CC1mu1pi1p_pi_pr_opening_angle_{reco,true}, CandidateProtonIndex; 2026-09-16) into STAGING
# directories, validate them bit-for-bit against the live files, and only then promote.
#
#   stage 1  the 41 processed/w units of ../logs/rerun_beta_inputs.txt -> processed/w_thpipr
#            (detVar part1+part2 merged first, entry count verified, as in rerun_beta.sh)
#   stage 2  the 7 per-run EXT files, BOTH selections -> processed/ext_perrun_thpipr
#            (shared with the inclusive tree, so the CC1mu1piXp branches are compared too)
#   stage 3  proton-tagged fake data re-thrown from the staging MC with the SAME seeds
#   stage 4  symlinks copied verbatim from processed/w (absolute targets resolve after promotion)
#   stage 5  macros/compare_reprocess.C on every regular file: all old branches bit-identical,
#            entry count and summed_pot unchanged, new angle == friend tree
#   stage 6  PROMOTE=1 only: rotate w -> w_pre_thpipr, ext_perrun -> ext_perrun_pre_thpipr and
#            move the staging trees into place. Refuses unless stage 5 is all OK.
#
# Idempotent: per-unit done marks in ../logs/thpipr/done.
#   usage: nohup setsid ./reprocess_thpipr.sh > ../logs/thpipr/driver.log 2>&1 &
#          PROMOTE=1 ./reprocess_thpipr.sh     (after reviewing stage 5)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u

PROC=/data/uboone/processed
LIVE_W=$PROC/w;          NEW_W=$PROC/w_thpipr
LIVE_E=$PROC/ext_perrun; NEW_E=$PROC/ext_perrun_thpipr
LOG=../logs/thpipr; DONEDIR=$LOG/done; STATUS=$LOG/status.txt
mkdir -p "$NEW_W" "$NEW_E" "$DONEDIR"; touch "$STATUS"
NICE="nice -n 15"; NPROC=${NPROC:-3}

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$STATUS"; }
key(){ printf '%s' "$1" | tr ' /' '__'; }
done_already(){ [ -f "$DONEDIR/$(key "$1")" ]; }
mark(){ : > "$DONEDIR/$(key "$1")"; }

proc_one(){  # rawpath outpath filetype selections tag
  local in="$1" out="$2" ft="$3" sel="$4" tag="$5" merged=""
  done_already "s1 $tag" && { echo "  [skip] $tag"; return 0; }
  if [[ "$in" == *"+"* ]]; then
    local p1="${in%%+*}" p2="${in##*+}"
    merged="$PROC/merge_tmp_thpipr_${tag}.root"; rm -f "$merged"
    $NICE hadd -f "$merged" "$p1" "$p2" > "$LOG/s1_${tag}_hadd.log" 2>&1 || { say "FAIL s1 $tag (hadd)"; rm -f "$merged"; return 1; }
    local nexp nget
    nexp=$(root.exe -l -b -q -e "Long64_t n=0;for(auto p:{\"$p1\",\"$p2\"}){auto f=TFile::Open(p);auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;if(t)n+=t->GetEntries();if(f)f->Close();}printf(\"%lld\\n\",n);" 2>/dev/null | tail -1)
    nget=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$merged\");auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;printf(\"%lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | tail -1)
    [[ "$nexp" == "$nget" ]] || { say "FAIL s1 $tag (merge truncated $nget/$nexp)"; rm -f "$merged"; return 1; }
    in="$merged"
  fi
  [[ -f "$in" ]] || { say "FAIL s1 $tag (no input $in)"; return 1; }
  $NICE ./bin/ProcessNTuples "$in" "$ft" "$sel" "$out" > "$LOG/s1_$tag.log" 2>&1
  local rc=$?
  [[ -n "$merged" ]] && rm -f "$merged"
  local chk=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%d\n\",(t&&t->GetEntries()>0&&t->GetBranch(\"CC1mu1pi1p_pi_pr_opening_angle_reco\"))?1:0);" 2>/dev/null | tail -1)
  if [[ $rc -eq 0 && "$chk" == "1" ]]; then mark "s1 $tag"; say "OK   s1 $tag"; else say "FAIL s1 $tag rc=$rc chk=$chk"; fi
}

if [ "${PROMOTE:-0}" != 1 ]; then
  say "===== STAGE 1+2: reprocessing (NPROC=$NPROC) ====="
  while IFS='|' read -r raw outdir ft sel tag; do
    [[ -z "${raw// }" || "${raw:0:1}" == "#" ]] && continue
    raw=$(echo $raw); outdir=$(echo $outdir); ft=$(echo $ft); sel=$(echo $sel); tag=$(echo $tag)
    [[ "$outdir" == "$LIVE_W" ]] || continue
    proc_one "$raw" "$NEW_W/xsec-ana-${tag#w_}.root" "$ft" "$sel" "$tag" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || sleep 5; done
  done < ../logs/rerun_beta_inputs.txt
  for s in neutrinoselection_filt_run1_beamoff neutrinoselection_filt_run3b_beamoff \
           numi_pelee_ntuple_beam_off_run4a_rhc_ana numi_pelee_ntuple_beam_off_run4b_rhc_ana \
           numi_pelee_ntuple_beam_off_run4c_fhc_ana numi_pelee_ntuple_beam_off_run4d_fhc_ana \
           numi_pelee_ntuple_beam_off_run5_fhc_ana; do
    proc_one "/data/uboone/EXT/$s.root" "$NEW_E/xsec-ana-$s.root" extBNB CC1mu1piXp,CC1mu1pi1p "ext_$s" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || sleep 5; done
  done
  wait

  say "===== STAGE 3: fake-data re-throw (seeds unchanged) ====="
  if ! done_already "s3 throw"; then
    export THROW_W_DIR="$NEW_W/"
    $NICE root.exe -l -b -q 'macros/throw_perrun_w.C(1)' > "$LOG/s3_fhc.log" 2>&1 && \
    $NICE root.exe -l -b -q macros/throw_perrun_w.C -e 'throw_perrun_w_rhc(1)' > "$LOG/s3_rhc.log" 2>&1 && \
      { mark "s3 throw"; say "OK   s3 throw"; } || say "FAIL s3 throw"
    unset THROW_W_DIR
  fi

  say "===== STAGE 4: symlinks ====="
  for l in "$LIVE_W"/*.root; do [ -L "$l" ] && cp -P "$l" "$NEW_W/"; done
  say "   $(find "$NEW_W" -maxdepth 1 -type l | wc -l) symlinks copied"

  say "===== STAGE 5: bit-for-bit comparison ====="
  cmp_one(){  # old new friend tag
    done_already "s5 $4" && return 0
    $NICE root.exe -l -b -q "macros/compare_reprocess.C+(\"$1\",\"$2\",\"$3\")" > "$LOG/s5_$4.log" 2>&1
    grep -q "^\[CMPREPRO\] .* OK" "$LOG/s5_$4.log" && { mark "s5 $4"; say "OK   s5 $4"; } || say "FAIL s5 $4"
  }
  root.exe -l -b -q -e '.L macros/compare_reprocess.C+' > /dev/null 2>&1   # compile once, not in parallel
  for f in "$NEW_W"/*.root; do
    [ -L "$f" ] && continue; b=$(basename "$f"); fr="$LIVE_W/friends_pa/${b%.root}.pa.root"; [ -f "$fr" ] || fr=""
    cmp_one "$LIVE_W/$b" "$f" "$fr" "w_$b" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || sleep 5; done
  done
  for f in "$NEW_E"/*.root; do
    b=$(basename "$f"); cmp_one "$LIVE_E/$b" "$f" "" "e_$b" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || sleep 5; done
  done
  wait
  nok=$(grep -l "^\[CMPREPRO\] .* OK" "$LOG"/s5_*.log 2>/dev/null | wc -l)
  nall=$(( $(find "$NEW_W" -maxdepth 1 -type f -name '*.root' | wc -l) + $(find "$NEW_E" -maxdepth 1 -type f -name '*.root' | wc -l) ))
  say "===== STAGE 5 done: $nok / $nall files identical ====="
  exit 0
fi

# ---------------------------------------------------------------- stage 6 (explicit)
nok=$(grep -l "^\[CMPREPRO\] .* OK" "$LOG"/s5_*.log 2>/dev/null | wc -l)
nall=$(( $(find "$NEW_W" -maxdepth 1 -type f -name '*.root' | wc -l) + $(find "$NEW_E" -maxdepth 1 -type f -name '*.root' | wc -l) ))
[ "$nall" -eq 56 ] && [ "$nok" -eq "$nall" ] || { say "!! REFUSING TO PROMOTE: $nok/$nall identical (expect 56/56)"; exit 1; }
[ -e "$PROC/w_pre_thpipr" ] || [ -e "$PROC/ext_perrun_pre_thpipr" ] && { say "!! REFUSING: a *_pre_thpipr rotation already exists"; exit 1; }
mv "$LIVE_W/friends_pa" "$NEW_W/friends_pa"
mv "$LIVE_W" "$PROC/w_pre_thpipr" && mv "$NEW_W" "$LIVE_W"
mv "$LIVE_E" "$PROC/ext_perrun_pre_thpipr" && mv "$NEW_E" "$LIVE_E"
say "===== PROMOTED: $LIVE_W and $LIVE_E now carry the opening angle; old trees kept as *_pre_thpipr ====="
