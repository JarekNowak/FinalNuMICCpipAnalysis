#!/bin/bash
L=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/promote_run1fix.log
until grep -q "RUN1FIX PROMOTE+REGEN DONE" $L; do sleep 300; done
echo "promotion finished $(date +%H:%M)"; tail -8 $L
until [ "$(grep -h '^OK' /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/rebuild_3352400_*.out 2>/dev/null | wc -l)" -ge 4 ] || ! squeue -u $USER -h -o "%i" | grep -q 3352400; do sleep 300; done
echo "totals job 3352400: OK $(grep -h '^OK' /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/rebuild_3352400_*.out 2>/dev/null | wc -l)/4 $(date +%H:%M)"
