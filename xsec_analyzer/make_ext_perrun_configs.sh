#!/usr/bin/env bash
# make_ext_perrun_configs.sh -- stage 2 of the EXT fix: per-RUN-PERIOD beam-off files.
# Beam-off data are cosmic-only, so the horn label on a sample is irrelevant (analyser,
# 2026-09-02): every beam-on run, FHC or RHC, uses the beam-off files taken in that run
# period, pooled within the run by the framework (ext gates summed per run, each file
# scaled by bnb_gates[run]/sum). Gate counts from the analyser's table:
#   run1  4582248.27 (also stands in for run 2, which has no processed beam-off; a
#         stand-in is scaled by ITS OWN gates)
#   run3b 18512842.15 + 14136286.5 = 32649128.65
#   run4  4a+4b = 17398802.925 - 8060024.70 = 9338778.225 (split evenly over 4a/4b),
#         4c 8060024.70, 4d 17432345.70   (all four pooled under run 4)
#   run5  19256341.475
# The same staged file (both selections) serves processed/ and processed/w/; a distinct
# symlink name per (file, run id) keeps the framework's per-file maps and the univmake
# subdirectories unique when a run appears twice (FHC run 4 and RHC run 14 in COMB).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/configs"
ST=/data/uboone/processed/ext_perrun
declare -A SRC=( [run1]=neutrinoselection_filt_run1_beamoff [run3b]=neutrinoselection_filt_run3b_beamoff
  [run4a]=numi_pelee_ntuple_beam_off_run4a_rhc_ana [run4b]=numi_pelee_ntuple_beam_off_run4b_rhc_ana
  [run4c]=numi_pelee_ntuple_beam_off_run4c_fhc_ana [run4d]=numi_pelee_ntuple_beam_off_run4d_fhc_ana
  [run5]=numi_pelee_ntuple_beam_off_run5_fhc_ana )
declare -A GATE=( [run1]=4582248.27 [run3b]=32649128.65 [run4a]=4669389.1125 [run4b]=4669389.1125 [run4c]=8060024.70 [run4d]=17432345.70 [run5]=19256341.475 )
for k in "${!SRC[@]}"; do [ -s "$ST/xsec-ana-${SRC[$k]}.root" ] || { echo "missing staged $k"; exit 1; }; done
# files per run id
files_for_run() { case "$1" in 1|11) echo run1;; 2|12) echo run1;; 3|13) echo run3b;; 4|14) echo "run4a run4b run4c run4d";; 5) echo run5;; esac; }
mklines() {  # dir runid
  local dir=$1 rid=$2 line
  for k in $(files_for_run "$rid"); do
    local ln="$dir/xsec-ana-ext_${k}_run${rid}.root"
    ln -sfn "$ST/xsec-ana-${SRC[$k]}.root" "$ln"
    echo "$ln $rid extBNB ${GATE[$k]} 0"
  done
}
hdr='# PER-RUN beam-off (EXT), 2026-09-02: cosmic-only data, matched by run PERIOD irrespective of horn
#   label; files pooled within a run (framework sums ext gates per run). Gate counts from the
#   analyser: run1 4582248.27 (stands in for run 2, scaled by its own gates); run3b pre+post
#   32649128.65; run4 4a+4b 9338778.225 (split evenly) + 4c 8060024.70 + 4d 17432345.70;
#   run5 19256341.475. Symlink names are unique per (file, run id).'
gen() {  # src dst dir runids...
  local src=$1 dst=$2 dir=$3; shift 3
  { local done=0
    while IFS= read -r line; do
      if [[ "$line" =~ [[:space:]]extBNB[[:space:]] && ! "$line" =~ ^# ]]; then
        if [ $done -eq 0 ]; then echo "$hdr"; for rid in "$@"; do mklines "$dir" "$rid"; done; done=1; fi
        continue
      fi
      echo "$line"
    done < "$src"; } > "$dst"
  echo "wrote $dst ($(grep -c extBNB "$dst") extBNB lines)"
}
P=/data/uboone/processed; W=/data/uboone/processed/w
gen file_properties_numi_fhc5.txt      file_properties_numi_fhc5_perrun.txt      $P 1 2 4 5
gen file_properties_numi_rhcfull.txt   file_properties_numi_rhcfull_perrun.txt   $P 1 2 3 4
gen file_properties_numi_comb.txt      file_properties_numi_comb_perrun.txt      $P 1 2 4 5 11 12 13 14
gen file_properties_numi_fhc5_w.txt    file_properties_numi_fhc5_w_perrun.txt    $W 1 2 4 5
gen file_properties_numi_rhcfull_w.txt file_properties_numi_rhcfull_w_perrun.txt $W 1 2 3 4
gen file_properties_numi_comb_w.txt    file_properties_numi_comb_w_perrun.txt    $W 1 2 4 5 11 12 13 14
