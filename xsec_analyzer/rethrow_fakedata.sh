#!/usr/bin/env bash
# rethrow_fakedata.sh -- re-throw every per-run fake-data set with the framework safe-weight
# rule (non-finite/negative/>30 CV weights -> 1, not skipped). Seeds unchanged.
cd "$(dirname "${BASH_SOURCE[0]}")"; set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
L=../logs/rethrow; mkdir -p $L
root -l -b -q 'macros/throw_perrun_fhc.C(1)' > $L/fhc.log 2>&1 && echo "OK fhc" || echo "FAIL fhc"
root -l -b -q 'macros/throw_perrun_rhc.C(1)' > $L/rhc.log 2>&1 && echo "OK rhc" || echo "FAIL rhc"
root -l -b -q 'macros/throw_perrun_w.C(1)'   > $L/w_fhc.log 2>&1 && echo "OK w_fhc" || echo "FAIL w_fhc"
root -l -b -q -e '.L macros/throw_perrun_w.C' -e 'throw_perrun_w_rhc(1)' > $L/w_rhc.log 2>&1 && echo "OK w_rhc" || echo "FAIL w_rhc"
root -l -b -q 'macros/throw_perrun_sb.C(1)'  > $L/sb.log 2>&1 && echo "OK sb" || echo "FAIL sb"
echo RETHROW_DONE
