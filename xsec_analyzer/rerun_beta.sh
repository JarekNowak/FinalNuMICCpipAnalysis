#!/usr/bin/env bash
# rerun_beta.sh -- the full re-run for the beam-frame (beta) correction, start to
# finish, unattended.
#
#   stage 1  reprocess every input ntuple with the beta-corrected selections
#   stage 1b re-throw the per-run Poisson fake data from the new MC
#   stage 2  univmake for every (config, observable) whose observable changed
#   stage 3  UnfolderNuMI extraction + figure regeneration
#
# IDEMPOTENT AND RESUMABLE. Every stage skips work whose output already exists and
# validates; killing and re-running picks up where it left off. Progress is appended
# to ../logs/rerun_beta/status.txt, one line per unit of work.
#
# What is deliberately NOT rebuilt: the univmake outputs for p_mu, p_pi, theta_mupi
# and W_pipr. Those observables are frame-independent and were verified unchanged
# event-by-event over 300k events (zero differences on signal events; the 8941
# differences are all on non-signal entries, where the truth branches are filled
# opportunistically and carry no meaning). Rebuilding them would cost ~1h20m each for
# an identical result. Pass REBUILD_ALL=1 to rebuild them anyway.
#
#   usage:  nohup setsid ./rerun_beta.sh > ../logs/rerun_beta/driver.log 2>&1 &
#           STAGE=2 ./rerun_beta.sh        # start from a given stage
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:${LD_LIBRARY_PATH:-}"
set +u; source ./setup_xsec_analyzer.sh 2>/dev/null; set -u

PROC=/data/uboone/processed
LOG=../logs/rerun_beta; mkdir -p "$LOG"
STATUS=$LOG/status.txt; touch "$STATUS"
NICE="nice -n 15"
NPROC=${NPROC:-3}          # parallel ProcessNTuples; /data is NFS, >3 does not help
NUNIV=${NUNIV:-2}          # parallel univmake, as in perrun_batch2.sh
START_STAGE=${STAGE:-1}
REBUILD_ALL=${REBUILD_ALL:-0}

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$STATUS"; }

# Completion marks are PER-UNIT FILES, not appends to a shared log.
#
# They used to be `echo "DONE $tag" >> $STATUS`. $STATUS lives on NFS, where a
# concurrent append from two parallel units finishing seconds apart is NOT atomic
# and one of the two lines is silently lost. The unit had really completed, but the
# driver had no record of it and rebuilt it from scratch on the next restart.
# Measured on the 2026-08-22 run: 36 "OK s2" lines against only 11 surviving DONE
# marks, 11 of 24 units built twice, roughly 8 hours of wasted compute.
#
# Creating a distinct file per unit is atomic on NFS, so marks cannot collide.
DONEDIR="$LOG/done"; mkdir -p "$DONEDIR"
_donekey(){ printf '%s' "$1" | tr ' /' '__'; }
done_already(){ [ -f "$DONEDIR/$(_donekey "$1")" ]; }
mark(){ : > "$DONEDIR/$(_donekey "$1")"; echo "DONE $1" >> "$STATUS"; }

# ---------------------------------------------------------------- stage 1
# Reprocess one ntuple and validate that the tree and a beta-dependent branch exist.
proc_one(){  # rawpath outpath filetype selection tag
  local in="$1" out="$2" ft="$3" sel="$4" tag="$5"
  done_already "s1 $tag" && { echo "  [skip] $tag"; return 0; }
  # detVar inputs arrive as "part1+part2" and must be merged first. The merge is kept
  # in a scratch file so the two parts are never modified.
  local merged=""
  if [[ "$in" == *"+"* ]]; then
    local p1="${in%%+*}" p2="${in##*+}"
    [[ -f "$p1" && -f "$p2" ]] || { say "FAIL s1 $tag (missing part)"; return 1; }
    # Scratch name keyed on the TAG, so the inclusive and proton-tagged units of the
    # same detVar never share it. They previously did, and the second unit both read a
    # half-written file and had it deleted underneath it by the first.
    merged="$PROC/merge_tmp_${tag}.root"
    rm -f "$merged"
    # NO -k here. "-k" means "skip corrupt or unreadable files and keep going", which
    # over NFS turns a transient read failure into a SILENTLY TRUNCATED merge: the
    # first attempt at this produced a detVar merge holding 96,900 of 923,357 events,
    # exit status 0, nothing in the log. Fail loudly instead, and verify the merged
    # entry count against the sum of the parts before processing anything.
    $NICE hadd -f "$merged" "$p1" "$p2" > "$LOG/s1_${tag}_hadd.log" 2>&1 \
      || { say "FAIL s1 $tag (hadd rc=$?)"; rm -f "$merged"; return 1; }
    local nexp nget
    nexp=$(root.exe -l -b -q -e "Long64_t n=0;for(auto p:{\"$p1\",\"$p2\"}){auto f=TFile::Open(p);auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;if(t)n+=t->GetEntries();if(f)f->Close();}printf(\"%lld\\n\",n);" 2>/dev/null | tail -1)
    nget=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$merged\");auto t=f?(TTree*)f->Get(\"nuselection/NeutrinoSelectionFilter\"):nullptr;printf(\"%lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | tail -1)
    if [[ "$nexp" != "$nget" ]]; then
      say "FAIL s1 $tag (merge truncated: got ${nget:-?} of ${nexp:-?} entries)"
      rm -f "$merged"; return 1
    fi
    echo "  [merge ok] $tag $nget entries"
    in="$merged"
  fi
  [[ -f "$in" ]] || { say "FAIL s1 $tag (no input: $in)"; return 1; }
  $NICE ./bin/ProcessNTuples "$in" "$ft" "$sel" "$out" > "$LOG/s1_$tag.log" 2>&1
  local rc=$?
  local chk=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");auto t=f&&!f->IsZombie()?(TTree*)f->Get(\"stv_tree\"):nullptr;printf(\"%d\n\",(t&&t->GetEntries()>0&&t->GetBranch(\"${sel}_candidate_muon_costh_reco\"))?1:0);" 2>/dev/null | tail -1)
  [[ -n "$merged" ]] && rm -f "$merged"
  if [[ $rc -eq 0 && "$chk" == "1" ]]; then mark "s1 $tag"; say "OK   s1 $tag"; else say "FAIL s1 $tag rc=$rc chk=$chk"; fi
}

stage1(){
  say "===== STAGE 1: reprocessing (NPROC=$NPROC) ====="
  # Each line: rawpath | outdir | filetype | selection | tag
  # The inclusive (CC1mu1piXp) set writes to processed/, the proton-tagged
  # (CC1mu1pi1p) set to processed/w/, matching the two file_properties families.
  while IFS='|' read -r raw outdir ft sel tag; do
    [[ -z "${raw// }" || "${raw:0:1}" == "#" ]] && continue
    raw=$(echo $raw); outdir=$(echo $outdir); ft=$(echo $ft); sel=$(echo $sel); tag=$(echo $tag)
    # The output name ALWAYS comes from the tag, never from the input basename. Two
    # cases previously broke: merged "part1+part2" inputs, whose basename is the part2
    # filename, and single-file detVar inputs whose output name differs from the input
    # name entirely (e.g. prod_numi_..._detvar_cv_run4_rhc_ana.root must be written as
    # detvar_run4rhc_CV.root, which is what file_properties expects). Deriving the name
    # from the tag covers both; the tag carries a "w_" prefix for the proton-tagged
    # family, stripped here.
    local oname="${tag#w_}.root"
    proc_one "$raw" "$outdir/xsec-ana-$oname" "$ft" "$sel" "$tag" &
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || sleep 5; done
  done < "../logs/rerun_beta_inputs.txt"
  wait
  # per-run EXT symlinks, both families
  for fam in "$PROC" "$PROC/w"; do
    [[ -f "$fam/xsec-ana-beamoff_run1Andrun3.root" ]] || continue
    for r in run1 run2 run4 run5; do ln -sf "$fam/xsec-ana-beamoff_run1Andrun3.root" "$fam/xsec-ana-beamoff_fhc_${r}.root"; done
    for r in run1 run2 run3 run4; do ln -sf "$fam/xsec-ana-beamoff_run1Andrun3.root" "$fam/xsec-ana-beamoff_rhc_${r}.root"; done
  done
  say "===== STAGE 1 done ====="
}

# ---------------------------------------------------------------- stage 1b
stage1b(){
  say "===== STAGE 1b: fake-data re-throw ====="
  # Four separate throws, not two. The integer argument to these macros is a SEED,
  # not a mode: calling throw_perrun_w.C(0) then (1) simply re-threw the same four
  # FHC proton-tagged files twice and left the inclusive family and all of RHC
  # untouched, still carrying the pre-beta kinematics.
  #   macros/throw_perrun_fhc.C   -> inclusive FHC        (processed/)
  #   macros/throw_perrun_rhc.C   -> inclusive RHC        (processed/)
  #   macros/throw_perrun_w.C     -> proton-tagged FHC    (processed/w/)
  #   throw_perrun_w_rhc()        -> proton-tagged RHC    (processed/w/), same file
  local calls=( "throw_perrun_fhc.C" "throw_perrun_rhc.C" "throw_perrun_w.C"
                "throw_perrun_w.C+throw_perrun_w_rhc" )
  local names=( "incl_fhc" "incl_rhc" "w_fhc" "w_rhc" )
  local i
  for i in "${!calls[@]}"; do
    local tag="throw_${names[$i]}"
    done_already "s1b $tag" && { echo "  [skip] $tag"; continue; }
    local spec="${calls[$i]}" cmd
    if [[ "$spec" == *"+"* ]]; then
      # load the file, then invoke the named function inside it
      cmd="macros/${spec%%+*}+ -e ${spec##*+}(1)"
      $NICE root.exe -l -b -q "macros/${spec%%+*}" -e "${spec##*+}(1)" > "$LOG/s1b_$tag.log" 2>&1
    else
      $NICE root.exe -l -b -q "macros/$spec(1)" > "$LOG/s1b_$tag.log" 2>&1
    fi
    if [ $? -eq 0 ]; then mark "s1b $tag"; say "OK   s1b $tag"; else say "FAIL s1b $tag"; fi
  done
}


# ---------------------------------------------------------------- stage gates
# A stage must never consume the output of an incomplete or stale earlier stage.
# The first pass of this re-run did exactly that: 20 detVar units failed in stage 1,
# stage 2 went ahead and built universes from the pre-beta detVar files still sitting
# on disk from 2026-08-09, and the detector systematic in all 47 extractions was
# silently invalid. Nothing in the pipeline noticed, because every individual univmake
# and unfold succeeded. Hence two gates: completeness and freshness.

# Newest mtime among the sources that determine the branch VALUES. Any processed file
# older than this was produced by different physics logic.
code_stamp(){
  local newest=0 t
  for f in src/selections/*.cxx include/XSecAnalyzer/Selections/*.hh \
           include/XSecAnalyzer/NuMIBeamFrame.hh include/XSecAnalyzer/Branches.hh \
           include/XSecAnalyzer/AnalysisEvent.hh src/utils/STVTools.cxx; do
    [ -f "$f" ] || continue
    t=$(stat -c %Y "$f"); (( t > newest )) && newest=$t
  done
  echo "$newest"
}

# Gate 1: every unit in the manifest carries a DONE mark.
gate_complete(){
  local missing=0 tag
  while IFS='|' read -r raw outdir ft sel tag; do
    [[ -z "${raw// }" || "${raw:0:1}" == "#" ]] && continue
    tag=$(echo $tag)
    if ! done_already "s1 $tag"; then
      [ "$missing" -lt 20 ] && say "   not done: s1 $tag"
      missing=$((missing+1))
    fi
  done < ../logs/rerun_beta_inputs.txt
  if (( missing )); then say "!! GATE FAILED: $missing stage-1 unit(s) incomplete"; return 1; fi
  say "GATE ok: all stage-1 units complete"
  return 0
}

# Gate 2: every processed file that stage 2 will read is newer than the code stamp.
# This is the gate that would have caught the detVar problem: those files existed and
# opened cleanly, they were simply built by ten-day-old code.
gate_fresh(){
  local stamp=$(code_stamp) stale=0 shown=0 p r t
  say "   code stamp: $(date -d @$stamp '+%Y-%m-%d %H:%M')"
  for fp in configs/file_properties_numi_fhc5.txt configs/file_properties_numi_rhcfull.txt \
            configs/file_properties_numi_comb.txt configs/file_properties_numi_fhc5_w.txt \
            configs/file_properties_numi_rhcfull_w.txt configs/file_properties_numi_comb_w.txt; do
    [ -f "$fp" ] || continue
    while read -r p rest; do
      [[ -z "$p" || "${p:0:1}" == "#" ]] && continue
      [[ "$p" == /data/uboone/* ]] || continue
      r=$(readlink -f "$p")
      if [ ! -f "$r" ]; then
        [ "$shown" -lt 15 ] && { say "   MISSING  $(basename $p)"; shown=$((shown+1)); }
        stale=$((stale+1)); continue
      fi
      t=$(stat -c %Y "$r")
      if (( t < stamp )); then
        [ "$shown" -lt 15 ] && { say "   STALE    $(basename $r)  ($(date -d @$t '+%m-%d %H:%M'))"; shown=$((shown+1)); }
        stale=$((stale+1))
      fi
    done < "$fp"
  done
  if (( stale )); then say "!! GATE FAILED: $stale processed file(s) missing or older than the code"; return 1; fi
  say "GATE ok: every processed file stage 2 reads is newer than the code"
  return 0
}

gate_before_stage2(){
  say "===== GATE: checking stage 1 before stage 2 ====="
  local bad=0
  gate_complete || bad=1
  gate_fresh    || bad=1
  if (( bad )); then
    say "!! REFUSING TO RUN STAGE 2. Completed work is kept in $STATUS;"
    say "!! fix the failures above and re-run this script to resume."
    return 1
  fi
  return 0
}

# ---------------------------------------------------------------- stage 2
univ_one(){  # fpm bincfg out tag
  local fpm="$1" bc="$2" out="$3" tag="$4"
  done_already "s2 $tag" && { echo "  [skip] $tag"; return 0; }
  local sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -gt 50000000 ] && [ "$REBUILD_ALL" != "1" ]; then mark "s2 $tag"; say "OK   s2 $tag (existing $((sz/1000000))MB)"; return 0; fi
  say "  start s2 $tag"
  # $NICE was defined but never applied here, so univmake ran at nice 5 and
  # competed with other users of this shared box. UNIV_THREADS (read by
  # univmake.C) caps each job's ROOT implicit-MT pool; without it a single job
  # sizes its pool to the whole machine and extra parallelism just thrashes.
  FPM="$fpm" BIN_CONFIG="$bc" OUT="$out" $NICE ./run_universe_maker.sh > "$LOG/s2_$tag.log" 2>&1
  sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -gt 50000000 ]; then mark "s2 $tag"; say "OK   s2 $tag ($((sz/1000000))MB)"; else say "FAIL s2 $tag (size $sz)"; fi
}
# setpot_dirt() REMOVED 2026-08-22.
#
# It used to overwrite the summed_pot parameter INSIDE the shared dirt file with a
# different value per configuration, so that one run1 FHC dirt sample could stand in
# for the whole exposure of FHC, RHC and COMB alike. The inflations were x2.70 /
# x18.31 / x6.07, whose RHC/FHC ratio of 0.147 disagrees with the true dirt POT ratio
# of 0.567, and COMB came out BELOW both of its own inputs (0.749 against FHC 0.885
# and RHC 0.800) which a combined exposure cannot do.
#
# Mode- and run-matched dirt now exists and is processed with its native POT, so the
# framework scales each run itself exactly as it does for numuMC. Runs with no
# dedicated sample use a symlinked stand-in (FHC run2 <- run1; RHC runs 1,2 <- run3b),
# which is the same pattern already used for EXT.

stage2(){
  say "===== STAGE 2: univmake (NUNIV=$NUNIV) ====="
  # --- inclusive: only the observables whose definition changed ---
  # p_mu, p_pi and theta_mupi are frame-independent and verified unchanged, so their
  # univmake outputs stay valid (REBUILD_ALL=1 overrides).
  # theta_pi is deliberately NOT built. It was introduced (a32f1b4) only as the
  # symmetric partner of theta_mu, which exists because A_C crushed the sharply
  # forward-peaked cos(theta_mu); that same commit records that theta_pi "is
  # essentially unchanged (never sharply peaked)". It is not unfolded by stage 3 and
  # is not referenced in the note, so building it costs ~4.7 h for nothing. The
  # configs are kept on disk should it ever be wanted.
  local INCL_OBS="costhmu costhpi thetamu"
  [ "$REBUILD_ALL" == "1" ] && INCL_OBS="pmu ppi2bin costhmu costhpi thmupi thetamu"
  for spec in "fhc5 FHC5" "rhcfull RHCFULL" "comb COMB"; do
    set -- $spec; local cfg=$1 TAG=$2
    for o in $INCL_OBS; do
      local bc=configs/ccpi_${o}_bin_config_opt.txt
      [ -f "$bc" ] || bc=configs/ccpi_${o}_bin_config.txt
      [ -f "$bc" ] || { say "SKIP s2 ${TAG}_${o} (no bin config)"; continue; }
      univ_one "configs/file_properties_numi_${cfg}.txt" "$bc" \
               "$PROC/ccpi_${TAG}_${o}_univmake.root" "${TAG}_${o}" &
      while [ "$(jobs -rp | wc -l)" -ge "$NUNIV" ]; do wait -n 2>/dev/null || sleep 10; done
    done
    wait
  done
  # --- proton-tagged: the four TKI observables at the NEW two-bin binning ---
  for spec in "fhc5 FHC5" "rhcfull RHCFULL" "comb COMB"; do
    set -- $spec; local cfg=$1 TAG=$2
    for o in dpt dalphat dphit pn; do
      univ_one "configs/file_properties_numi_${cfg}_w.txt" "configs/ccpi1p_${o}_bin_config_2bin.txt" \
               "$PROC/ccpi1p_${TAG}_${o}2bin_univmake.root" "${TAG}_${o}2bin" &
      while [ "$(jobs -rp | wc -l)" -ge "$NUNIV" ]; do wait -n 2>/dev/null || sleep 10; done
    done
    wait
    # Proton-tagged angular/kinematic observables.
    #
    # This list MUST cover every observable stage 3 unfolds for this family, which is
    # "Wpipr Whad costhmu costhpi pmu ppi thmupi". It previously held only costhmu and
    # costhpi, so stage 3 unfolded FIVE observables whose univmake outputs stage 2 had
    # never built -- stale files from an older run. That was invisible until the dirt
    # change, when stage 3 failed with
    #   Missing TDirectoryFile ...xsec-ana-dirt_fhc_run1.root
    # because the stale outputs predate the per-run dirt entries. 15 units failed
    # (5 observables x 3 configs).
    for o in costhmu costhpi Wpipr Whad pmu ppi thmupi; do
      univ_one "configs/file_properties_numi_${cfg}_w.txt" "configs/ccpi1p_${o}_bin_config.txt" \
               "$PROC/ccpi1p_${TAG}_${o}_univmake.root" "1p_${TAG}_${o}" &
      while [ "$(jobs -rp | wc -l)" -ge "$NUNIV" ]; do wait -n 2>/dev/null || sleep 10; done
    done
    wait
  done
  say "===== STAGE 2 done ====="
}

# ---------------------------------------------------------------- stage 3
# ---- stage-2/stage-3 consistency guard -------------------------------------
# Stage 3 must never unfold a univmake file that stage 2 did not successfully
# build.  Three separate bugs (ppi2bin naming, the 50 MB size gate, and the
# proton-tagged observable loop) let stage 3 consume stale or missing inputs and
# only surfaced hours into a run.  The mark name is derived from the UnivFile in
# the xsec config -- the one string both stages already agree on -- so there is
# no third copy of the observable list to drift.
# ROOT creates the output at univmake START, so file existence is NOT proof of a
# finished build; only the mark, written after a clean exit, is.
s2_mark_for(){ local b; b=$(basename "$1" _univmake.root)
  case "$b" in
    ccpi1p_*2bin) echo "s2_${b#ccpi1p_}" ;;
    ccpi1p_*)     echo "s2_1p_${b#ccpi1p_}" ;;
    *)            echo "s2_${b#ccpi_}" ;;
  esac; }
unit_ready(){   # xseccfg -> 0 ready, 1 not; prints reason on failure
  local xc="$1" uf fp mk
  uf=$(awk '$1=="UnivFile"{print $2}' "$xc"); fp=$(awk '$1=="FPFile"{print $2}' "$xc")
  [ -n "$uf" ] || { echo "no UnivFile line in $(basename "$xc")"; return 1; }
  # s2_mark_for already returns the full marker basename (s2_<TAG>_<obs>), which is
  # exactly what mark()/_donekey() produce for "s2 <tag>" -- do NOT re-prefix it.
  mk="$DONEDIR/$(s2_mark_for "$uf")"
  [ -f "$mk" ] || { echo "stage-2 never completed ($(basename "$uf"))"; return 1; }
  [ -f "$uf" ] || { echo "marked but univmake missing ($(basename "$uf"))"; return 1; }
  if [ -n "$fp" ] && [ -f "$fp" ] && [ "$fp" -nt "$uf" ]; then
    echo "STALE: $(basename "$fp") newer than $(basename "$uf")"; return 1; fi
  return 0
}

unfold_one(){  # xseccfg slicecfg out tag
  local xc="$1" sc="$2" out="$3" tag="$4"
  done_already "s3 $tag" && { echo "  [skip] $tag"; return 0; }
  [ -f "$xc" ] && [ -f "$sc" ] || { say "SKIP s3 $tag (missing config)"; return 0; }
  local _why; if ! _why=$(unit_ready "$xc"); then
    say "SKIP s3 $tag (stage-2 not ready: $_why)"; return 0; fi
  ./bin/UnfolderNuMI "$xc" "$sc" "$out" > "$LOG/s3_$tag.log" 2>&1
  local s=$(grep -oE "SYSTDUMP\] sigma_int [0-9.eE+-]+" "$LOG/s3_$tag.log" | awk '{print $3}')
  if [ -n "$s" ]; then mark "s3 $tag"; say "OK   s3 $tag sigma_int=$s"; else say "FAIL s3 $tag"; fi
  rm -f "$out"
}

stage3(){
  say "===== STAGE 3: extraction ====="
  # Consistency summary before anything unfolds.  Never fatal: the per-unit
  # unit_ready() guard in unfold_one is what actually protects each unit, and
  # check_stages.sh exits 1 (problems) / 2 (stage 2 still running) by design,
  # so "|| true" keeps a future "set -e" from turning this into an outage.
  if [ -x ./check_stages.sh ]; then
    local _cs; _cs=$(./check_stages.sh 2>&1) || true
    printf '%s\n' "$_cs" >> "$STATUS"
    say "stage2/3 consistency: $(printf '%s\n' "$_cs" | grep -E '^  ready ' || echo 'unavailable')"
  fi
  for spec in "fhc5 FHC5" "rhcfull RHCFULL" "comb COMB"; do
    set -- $spec; local cfg=$1 TAG=$2
    for o in pmu costhmu costhpi thmupi thetamu; do
      unfold_one "configs/ccpi_xsec_config_numi_${o}_${cfg}.txt" "configs/ccpi_${o}_slice_config_opt.txt" \
                 "$PROC/xsec_${TAG}_${o}.root" "${cfg}_${o}"
    done
    unfold_one "configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt" "configs/ccpi_ppi_slice_config_2bin.txt" \
               "$PROC/xsec_${TAG}_ppi2bin.root" "${cfg}_ppi2bin"
    for o in dpt dalphat dphit pn; do
      unfold_one "configs/ccpi1p_xsec_config_numi_${o}2bin_${cfg}.txt" "configs/ccpi1p_${o}_slice_config_2bin.txt" \
                 "$PROC/xsec1p_${TAG}_${o}2bin.root" "1p_${cfg}_${o}2bin"
    done
    for o in Wpipr Whad costhmu costhpi pmu ppi thmupi; do
      unfold_one "configs/ccpi1p_xsec_config_numi_${o}_${cfg}.txt" "configs/ccpi1p_${o}_slice_config.txt" \
                 "$PROC/xsec1p_${TAG}_${o}.root" "1p_${cfg}_${o}"
    done
  done
  say "===== STAGE 3 done ====="
  say "  regenerating figures"
  ./regen_figs.sh      > "$LOG/regen_figs.log" 2>&1      && say "OK   regen_figs"      || say "FAIL regen_figs"
  ./reco_spectra_all.sh > "$LOG/reco_spectra.log" 2>&1   && say "OK   reco_spectra"    || say "FAIL reco_spectra"
}

say "######## RERUN BETA START (from stage $START_STAGE) ########"
[ "$START_STAGE" -le 1 ] && { stage1; stage1b; }
if [ "$START_STAGE" -le 2 ]; then
  if gate_before_stage2; then stage2; else
    say "######## RERUN BETA HALTED AT THE GATE ########"; exit 2
  fi
fi
# Stage 3 reads the univmake outputs, so it inherits the same gate by construction:
# stage 2 cannot have run without passing it.
[ "$START_STAGE" -le 3 ] && stage3
say "######## RERUN BETA COMPLETE ########"
grep -c "^DONE" "$STATUS" | xargs echo "  units completed:"
grep "FAIL" "$STATUS" | sort -u
