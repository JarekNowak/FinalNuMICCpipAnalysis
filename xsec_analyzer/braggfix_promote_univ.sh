#!/usr/bin/env bash
# braggfix_promote_univ.sh -- install the 65 universe files rebuilt without the Bragg-pion cut
# (slurm/slurm_braggfix_univmake.sbatch -> processed/rebuild_braggfix). Refuses unless every task of the
# array reported identical bin specifications. Live files are moved to
# processed/univmake_backup_pre_braggfix/ (same relative path), the rebuilt ones moved into place; every
# move is logged in ../logs/braggfix/promote_univ.log.
#   JOB=<univmake array id> ./braggfix_promote_univ.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PROC=/data/uboone/processed; RB=$PROC/rebuild_braggfix; BK=$PROC/univmake_backup_pre_braggfix; L=../logs/braggfix
JOB="${JOB:?set JOB to the univmake array id}"
M=slurm/braggfix_univ_manifest.list; n=$(grep -c '' $M)
nok=$(grep -l "^bin specs identical" slurm/bfunivmake_${JOB}_*.out 2>/dev/null | xargs -r grep -l "^OK " | wc -l)
[ "$nok" -eq "$n" ] || { echo "REFUSING: $nok/$n univmake tasks OK"; exit 1; }
[ -e "$BK" ] && { echo "REFUSING: $BK exists"; exit 1; }
while IFS='|' read -r live fpm bins; do
  new="$RB/$(basename "$live")"; [ -s "$new" ] || { echo "missing $new"; exit 1; }
  rel="${live#$PROC/}"; mkdir -p "$BK/$(dirname "$rel")"
  mv "$live" "$BK/$rel"; echo "mv $live $BK/$rel" >> $L/promote_univ.log
  mv "$new" "$live";     echo "mv $new $live" >> $L/promote_univ.log
done < $M
echo "promoted $n universe files; old ones in $BK"
