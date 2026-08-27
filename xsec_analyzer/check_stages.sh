#!/usr/bin/env bash
# check_stages.sh -- verify every unit stage 3 intends to unfold has a stage-2 output
# that is BUILT (marked), FRESH, and (with --deep) COMPLETE.  Run before stage 3.
#
#   ./check_stages.sh          fast: mark + existence + freshness
#   ./check_stages.sh --deep   also opens each univmake file and counts TDirectoryFiles
#
# Three bugs of this exact shape have hit this chain, each surfacing hours into a run:
#   1. ppi2bin   -- stage 3 unfolded a file stage 2 never built (config named
#                   ccpi_ppi_bin_config_2bin.txt vs observable "ppi2bin"); silently
#                   used a stale file rather than failing.
#   2. 50 MB gate -- valid outputs marked FAIL on a byte-count heuristic.
#   3. proton-tagged Wpipr/Whad/pmu/ppi/thmupi -- stage 2 built 2 of the 7 observables
#                   stage 3 unfolds; the other 5 read files predating the per-run dirt.
#
# Design note: the unit list below mirrors STAGE 3's loops (the consumer), and the
# stage-2 mark name is DERIVED from the UnivFile path in each xsec config rather than
# from stage 2's loops.  Both stages already agree on the UnivFile, so there is no
# third copy of the observable list to drift out of sync -- which is what caused all
# three bugs above.
#
# Why the mark and not just the file: ROOT creates the output at univmake START and
# fills it incrementally, so an in-progress or crashed job leaves a plausible-looking
# partial file on disk.  mark() is written only after univmake exits successfully.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
DEEP=0; FILTER=""
while [ $# -gt 0 ]; do
  case "$1" in
    --deep) DEEP=1 ;;
    --only) shift; FILTER="${1:-}" ;;   # substring match on the unit short name
    *) echo "unknown arg: $1" >&2; exit 3 ;;
  esac; shift
done
DONEDIR=../logs/rerun_beta/done
bad=0; ok=0; building=0

# mark name from the univmake path -- mirrors stage 2's tag construction:
#   ccpi_<TAG>_<o>_univmake.root       -> s2_<TAG>_<o>
#   ccpi1p_<TAG>_<o>2bin_univmake.root -> s2_<TAG>_<o>2bin
#   ccpi1p_<TAG>_<o>_univmake.root     -> s2_1p_<TAG>_<o>
markof(){ local b; b=$(basename "$1" _univmake.root)
  case "$b" in
    ccpi1p_*2bin) echo "s2_${b#ccpi1p_}" ;;
    ccpi1p_*)     echo "s2_1p_${b#ccpi1p_}" ;;
    *)            echo "s2_${b#ccpi_}" ;;
  esac; }

units=()
for cfg in fhc5 rhcfull comb; do
  for o in pmu costhmu costhpi thmupi thetamu; do units+=("configs/ccpi_xsec_config_numi_${o}_${cfg}.txt"); done
  units+=("configs/ccpi_xsec_config_numi_ppi2bin_${cfg}.txt")
  for o in dpt dalphat dphit pn;                      do units+=("configs/ccpi1p_xsec_config_numi_${o}2bin_${cfg}.txt"); done
  for o in Wpipr Whad costhmu costhpi pmu ppi thmupi; do units+=("configs/ccpi1p_xsec_config_numi_${o}_${cfg}.txt"); done
done

printf "  %-40s %s\n" "unit" "status"
printf "  %-40s %s\n" "----------------------------------------" "------"
for xc in "${units[@]}"; do
  short=$(basename "$xc" .txt | sed 's/_xsec_config_numi//')
  if [ -n "$FILTER" ]; then case "$short" in *"$FILTER"*) ;; *) continue ;; esac; fi
  if [ ! -f "$xc" ]; then printf "  %-40s NO XSEC CONFIG\n" "$short"; bad=$((bad+1)); continue; fi
  uf=$(awk '$1=="UnivFile"{print $2}' "$xc"); fp=$(awk '$1=="FPFile"{print $2}' "$xc")
  if [ -z "$uf" ]; then printf "  %-40s no UnivFile line\n" "$short"; bad=$((bad+1)); continue; fi
  mk="$DONEDIR/$(markof "$uf")"
  if [ ! -f "$mk" ]; then
    if [ -f "$uf" ]; then
      printf "  %-40s UNMARKED (partial/in-progress: %s)\n" "$short" "$(du -h "$uf" 2>/dev/null|cut -f1)"
      building=$((building+1))
    else
      printf "  %-40s NOT BUILT (no mark, no file)\n" "$short"; bad=$((bad+1))
    fi; continue
  fi
  if [ ! -f "$uf" ]; then printf "  %-40s MARKED BUT FILE GONE: %s\n" "$short" "$(basename "$uf")"; bad=$((bad+1)); continue; fi
  if [ -n "$fp" ] && [ -f "$fp" ] && [ "$fp" -nt "$uf" ]; then
    printf "  %-40s STALE (file_properties newer than univmake)\n" "$short"; bad=$((bad+1)); continue; fi
  if [ "$DEEP" -eq 1 ]; then
    # skip comment lines: headers mention "xsec-ana-*.root", and the real beam-on
    # data line is deliberately commented out for the blind analysis -- counting
    # either would fake an INCOMPLETE verdict.
    nexp=$(grep -vE '^[[:space:]]*#' "$fp" 2>/dev/null | grep -cE '\.root' || echo 0)
    ngot=$(root -l -b -q -e "auto f=TFile::Open(\"$uf\");
      if(!f||f->IsZombie()){printf(\"-1\n\");}else{
        TString dn=((TKey*)f->GetListOfKeys()->At(0))->GetName();
        auto d=(TDirectory*)f->Get(dn); int n=0;
        if(d) for(auto k : *d->GetListOfKeys()) if(TString(((TKey*)k)->GetClassName()).Contains(\"TDirectory\")) n++;
        printf(\"%d\n\",n);}" 2>/dev/null | tail -1)
    if [ "${ngot:-0}" -lt "$nexp" ]; then
      printf "  %-40s INCOMPLETE (%s dirs, expected >= %s)\n" "$short" "$ngot" "$nexp"; bad=$((bad+1)); continue; fi
  fi
  ok=$((ok+1))
done
echo
echo "  ready $ok | in-progress $building | problems $bad   (of ${#units[@]} units)"
[ $building -gt 0 ] && echo "  note: in-progress units are stage 2 still running -- not errors, but stage 3 must wait"
if [ $bad -gt 0 ]; then echo "  -> stage 3 would FAIL or use bad inputs for the units above"; exit 1; fi
[ $building -gt 0 ] && exit 2
echo "  -> stage 2 covers everything stage 3 unfolds"
