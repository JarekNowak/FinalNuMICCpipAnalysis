#!/usr/bin/env bash
# Runs the pi0-definition comparison on the sb_pi0/ reprocess: yields/purity/overlap/shapes
# (sb_pi0_final_yields.C), transfer factors and their multisim uncertainty (sb_transfer_pf /
# sb_transfer_unc_pf), the detVar term and its MC-stat floor (sb_detvar_pf / sb_detvar_split_pf),
# and the SAME macros with the frozen flag on the same files (_p0) as the consistency check
# against the sb/ numbers in the note. Logs: logs/sidebands_perrun/pi0final/.
cd "$(dirname "${BASH_SOURCE[0]}")/../.."; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
L=../logs/sidebands_perrun/pi0final; mkdir -p $L
for m in fhc rhc; do
  root.exe -l -b -q "macros/pi0final/sb_pi0_final_yields.C(\"$m\")" > $L/yields_$m.log 2>&1
  for v in pf p0; do for mac in sb_transfer sb_transfer_unc sb_detvar sb_detvar_split; do
    root.exe -l -b -q "macros/pi0final/${mac}_$v.C(\"$m\")" > $L/${mac}_${v}_$m.log 2>&1
  done; done
  echo "$m done $(date +%H:%M)"
done; echo PI0FINAL_DONE
