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
NICE="nice -n 12"
NPROC=${NPROC:-3}          # parallel ProcessNTuples; /data is NFS, >3 does not help
NUNIV=${NUNIV:-2}          # parallel univmake, as in perrun_batch2.sh
START_STAGE=${STAGE:-1}
REBUILD_ALL=${REBUILD_ALL:-0}

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$STATUS"; }
done_already(){ grep -qxF "DONE $1" "$STATUS"; }
mark(){ echo "DONE $1" >> "$STATUS"; }

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
    merged="$PROC/merge_tmp_$(basename ${out%.root}).root"
    if [[ ! -s "$merged" ]]; then
      $NICE hadd -f -k "$merged" "$p1" "$p2" > "$LOG/s1_${tag}_hadd.log" 2>&1 \
        || { say "FAIL s1 $tag (hadd)"; rm -f "$merged"; return 1; }
    fi
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
    proc_one "$raw" "$outdir/xsec-ana-$(basename $raw)" "$ft" "$sel" "$tag" &
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
  for m in 0 1; do
    local tag="throw_$m"
    done_already "s1b $tag" && { echo "  [skip] $tag"; continue; }
    $NICE root.exe -l -b -q "macros/throw_perrun_w.C($m)" > "$LOG/s1b_$tag.log" 2>&1 \
      && { mark "s1b $tag"; say "OK   s1b $tag"; } || say "FAIL s1b $tag"
  done
}

# ---------------------------------------------------------------- stage 2
univ_one(){  # fpm bincfg out tag
  local fpm="$1" bc="$2" out="$3" tag="$4"
  done_already "s2 $tag" && { echo "  [skip] $tag"; return 0; }
  local sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -gt 50000000 ] && [ "$REBUILD_ALL" != "1" ]; then mark "s2 $tag"; say "OK   s2 $tag (existing $((sz/1000000))MB)"; return 0; fi
  say "  start s2 $tag"
  FPM="$fpm" BIN_CONFIG="$bc" OUT="$out" ./run_universe_maker.sh > "$LOG/s2_$tag.log" 2>&1
  sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  if [ "$sz" -gt 50000000 ]; then mark "s2 $tag"; say "OK   s2 $tag ($((sz/1000000))MB)"; else say "FAIL s2 $tag (size $sz)"; fi
}
setpot_dirt(){  # summed_pot  [dirtfile]
  local d=${2:-$PROC/xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root}
  root.exe -l -b -q -e "TFile*f=TFile::Open(\"$d\",\"update\");TParameter<float> p(\"summed_pot\",(float)$1);p.Write(\"summed_pot\",TObject::kOverwrite);f->Close();" >/dev/null 2>&1
}

stage2(){
  say "===== STAGE 2: univmake (NUNIV=$NUNIV) ====="
  # --- inclusive: only the observables whose definition changed ---
  # p_mu, p_pi and theta_mupi are frame-independent and verified unchanged, so their
  # univmake outputs stay valid (REBUILD_ALL=1 overrides).
  local INCL_OBS="costhmu costhpi thetamu thetapi"
  [ "$REBUILD_ALL" == "1" ] && INCL_OBS="pmu ppi2bin costhmu costhpi thmupi thetamu thetapi"
  for spec in "fhc5 FHC5 6.2046e20" "rhcfull RHCFULL 9.1429e19" "comb COMB 2.757e20"; do
    set -- $spec; local cfg=$1 TAG=$2 dirt=$3
    setpot_dirt "$dirt"
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
  for spec in "fhc5 FHC5 6.2046e20" "rhcfull RHCFULL 9.1429e19" "comb COMB 2.757e20"; do
    set -- $spec; local cfg=$1 TAG=$2 dirt=$3
    setpot_dirt "$dirt"
    for o in dpt dalphat dphit pn; do
      univ_one "configs/file_properties_numi_${cfg}_w.txt" "configs/ccpi1p_${o}_bin_config_2bin.txt" \
               "$PROC/ccpi1p_${TAG}_${o}2bin_univmake.root" "${TAG}_${o}2bin" &
      while [ "$(jobs -rp | wc -l)" -ge "$NUNIV" ]; do wait -n 2>/dev/null || sleep 10; done
    done
    wait
    # proton-tagged angular observables also moved
    for o in costhmu costhpi; do
      univ_one "configs/file_properties_numi_${cfg}_w.txt" "configs/ccpi1p_${o}_bin_config.txt" \
               "$PROC/ccpi1p_${TAG}_${o}_univmake.root" "1p_${TAG}_${o}" &
      while [ "$(jobs -rp | wc -l)" -ge "$NUNIV" ]; do wait -n 2>/dev/null || sleep 10; done
    done
    wait
  done
  say "===== STAGE 2 done ====="
}

# ---------------------------------------------------------------- stage 3
unfold_one(){  # xseccfg slicecfg out tag
  local xc="$1" sc="$2" out="$3" tag="$4"
  done_already "s3 $tag" && { echo "  [skip] $tag"; return 0; }
  [ -f "$xc" ] && [ -f "$sc" ] || { say "SKIP s3 $tag (missing config)"; return 0; }
  ./bin/UnfolderNuMI "$xc" "$sc" "$out" > "$LOG/s3_$tag.log" 2>&1
  local s=$(grep -oE "SYSTDUMP\] sigma_int [0-9.eE+-]+" "$LOG/s3_$tag.log" | awk '{print $3}')
  if [ -n "$s" ]; then mark "s3 $tag"; say "OK   s3 $tag sigma_int=$s"; else say "FAIL s3 $tag"; fi
  rm -f "$out"
}

stage3(){
  say "===== STAGE 3: extraction ====="
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
[ "$START_STAGE" -le 2 ] && stage2
[ "$START_STAGE" -le 3 ] && stage3
say "######## RERUN BETA COMPLETE ########"
grep -c "^DONE" "$STATUS" | xargs echo "  units completed:"
grep "FAIL" "$STATUS" | sort -u
