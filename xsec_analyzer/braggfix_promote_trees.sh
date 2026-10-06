#!/usr/bin/env bash
# braggfix_promote_trees.sh -- install the staging trees reprocessed without the Bragg-pion cut. Refuses
# unless ../logs/braggfix/gate.tsv has every file OK (101: 47 inclusive, 47 proton-tagged, 7 beam-off).
# Live regular files are moved to processed/pre_braggfix/{incl,w,ext_perrun}/ and the staging files moved
# into their place; symlinks (beam-off and dirt aliases) are not touched, their targets are the
# replaced files. Every move is logged in ../logs/braggfix/promote_trees.log (reverse it to undo).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PROC=/data/uboone/processed; L=../logs/braggfix; G=$L/gate.tsv; BK=$PROC/pre_braggfix
[ -f "$G" ] || { echo "no gate"; exit 1; }
nok=$(grep -c "^OK" "$G"); nall=$(grep -c "" "$G")
[ "$nall" -eq 101 ] && [ "$nok" -eq "$nall" ] || { echo "REFUSING: gate $nok/$nall OK (expect 101/101)"; exit 1; }
[ -e "$BK" ] && { echo "REFUSING: $BK exists"; exit 1; }
mkdir -p "$BK/incl" "$BK/w" "$BK/ext_perrun"
for pair in "incl_braggfix:.:incl" "w_braggfix:w:w" "ext_perrun_braggfix:ext_perrun:ext_perrun"; do
  IFS=: read -r st lv bk <<< "$pair"
  for f in "$PROC/$st"/*.root; do
    b=$(basename "$f"); live="$PROC/$lv/$b"
    [ -L "$live" ] && { echo "REFUSING: live $live is a symlink"; exit 1; }
    if [ -f "$live" ]; then mv "$live" "$BK/$bk/$b"; echo "mv $live $BK/$bk/$b" >> $L/promote_trees.log; fi
    mv "$f" "$live"; echo "mv $f $live" >> $L/promote_trees.log
  done
done
echo "promoted $(grep -c '' $L/promote_trees.log) moves; old files in $BK"
