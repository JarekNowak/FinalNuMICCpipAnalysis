#!/bin/bash
# detVar processing campaign: chain <stem>*.root from the MERGED detvars folder
# (handles part1/part2/... splits) -> processed detVar file. Does NOT touch /temp.
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer
export XSEC_ANALYZER_DIR="$PWD"
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$PWD/lib:$LD_LIBRARY_PATH"
source ./setup_xsec_analyzer.sh 2>/dev/null

MAP=configs/detvar_campaign_map.txt
DVDIR=/data/uboone/detvars
LOGDIR=../logs/detvar; mkdir -p "$LOGDIR"
OUTDIR=/data/uboone/processed
STATUS=$LOGDIR/campaign_status.txt
: > "$STATUS"
MAXJOBS=6

process_one() {
  local stem="$1" type="$2" tag="$3"
  local out="$OUTDIR/xsec-ana-detvar_${tag}.root"
  local log="$LOGDIR/${tag}.log"
  ProcessNTuples "$DVDIR/${stem}*.root" "$type" CC1mu1piXp "$out" > "$log" 2>&1
  local rc=$?
  local chk=$(root.exe -l -b -q -e "auto f=TFile::Open(\"$out\");if(!f||f->IsZombie()){printf(\"BAD\n\");}else{auto t=(TTree*)f->Get(\"stv_tree\");printf(\"ev=%lld\n\",t?t->GetEntries():-1);}" 2>/dev/null | grep -E "ev=|BAD")
  if [[ $rc -eq 0 && "$chk" != "BAD" && "$chk" != "ev=-1" && "$chk" != "ev=0" ]]; then
    echo "OK   $tag  $type  $chk" >> "$STATUS"
  else
    echo "FAIL $tag  $type  rc=$rc  $chk" >> "$STATUS"
  fi
}
export -f process_one
export DVDIR OUTDIR LOGDIR STATUS

njobs=0
while read stem type tag; do
  [[ "$stem" =~ ^# || -z "$stem" ]] && continue
  process_one "$stem" "$type" "$tag" &
  njobs=$((njobs+1))
  while (( $(jobs -r | wc -l) >= MAXJOBS )); do wait -n 2>/dev/null || sleep 2; done
done < "$MAP"
wait
echo "########## DETVAR CAMPAIGN DONE ##########" >> "$STATUS"
sort "$STATUS"
