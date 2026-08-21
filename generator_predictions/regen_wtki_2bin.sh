#!/usr/bin/env bash
# regen_wtki_2bin.sh -- rebuild the four generator W/TKI predictions at the TWO-bin
# analysis binning (wtki_gen.h edges, see configs/ccpi1p_TKI_binning.README) and
# publish them under the *_wtki_2bin_fte.root names the 2bin xsec configs expect.
# The six-bin versions are preserved as *_wtki_6bin_fte.root and restored at the end,
# so *_wtki_fte.root keeps meaning "the binning the six-bin configs use".
set -uo pipefail
cd /home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions
./gen_wtki_regen.sh > ../logs/gen_wtki_2bin.log 2>&1
rc=$?
for g in genie gibuu neut nuwro; do
  if [ -f newg4/${g}_wtki_fte.root ]; then
    cp newg4/${g}_wtki_fte.root newg4/${g}_wtki_2bin_fte.root
    echo "published ${g}_wtki_2bin_fte.root"
  else
    echo "MISSING ${g}_wtki_fte.root after regen"
  fi
done
# restore the six-bin files under the plain name
for g in genie gibuu neut nuwro; do
  [ -f newg4/${g}_wtki_6bin_fte.root ] && cp newg4/${g}_wtki_6bin_fte.root newg4/${g}_wtki_fte.root
done
echo "regen rc=$rc done $(date '+%H:%M:%S')"
