#!/usr/bin/env bash
# promote_050.sh -- install the binnings adopted under the 0.50 migration criterion (user decision 2026-09-26,
# after the binning review of the technical supplement, sec:binreview) and regenerate the release.
#   1. check: the 27 universe files built by slurm/slurm_adopt050.sbatch (SLURM 3433061) and the gen2d
#      prediction histograms the new cross-section configs read;
#   2. back up every live config and universe file that is replaced, and archive the release products of
#      the retired extractions (inclusive two-bin p_pi -> ppi3bin, theta_mu dropped, proton-tagged
#      two-bin delta phi_T -> dphit3bin), so no tool can pick up a stale file;
#   3. install configs/adopt050/* and the universe files;
#   4. run regen_release.sh (re-unfold the 48 extractions, covariances, release export, figures, counting,
#      cut-flows) and re-unfold the proton-angle study extractions (theta_p, theta_pip, outside the 48).
# The report tables (closure_tables.C, rebuild_tables.py, ...) are run afterwards, once the documents
# carry the new table structure.
# PP_K (required): bins of the new proton-momentum observable p_p, chosen from its review candidates; its
# configs are staged by scripts/adopt_050_configs.py $PP_K and its universe files are the candidates' own
# (the unfolder reads the first directory of a universe file, so the file name is all that changes).
set -uo pipefail
PP_K=${PP_K:?set PP_K to the chosen number of p_p bins}
REPO=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer; cd "$REPO"
PROC=/data/uboone/processed; RB=$PROC/rebuild_adopt050; STAMP=pre050_20260926
GP=/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d
REL=../report/data_release
T3="FHC5 RHCFULL COMB"

# ---- 1. checks -----------------------------------------------------------------------------------------
miss=0
while read -r PFX CFG OBS BIN; do
  case $CFG in fhc5) T=FHC5;; rhcfull) T=RHCFULL;; comb) T=COMB;; esac
  [ -s "$RB/${PFX}_${T}_${OBS}_univmake.root" ] || { echo "MISSING $RB/${PFX}_${T}_${OBS}_univmake.root"; miss=1; }
done < slurm/adopt050_manifest.list
for T in $T3; do [ -s "$PROC/rebuild_binreview/ccpi1p_${T}_ppK${PP_K}_univmake.root" ] || { echo "MISSING p_p universe file $T"; miss=1; }; done
for c in fhc5 rhcfull comb; do [ -f "configs/adopt050/ccpi1p_xsec_config_numi_pp_${c}.txt" ] || { echo "MISSING staged p_p config $c (run scripts/adopt_050_configs.py $PP_K)"; miss=1; }; done
PP_K=$PP_K python3 - <<'PYEOF' || miss=1
import uproot, sys
G='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d'; bad=0
for g in ('genie','gibuu','neut','nuwro'):
    for m in ('fhc','rhc','comb'):
        for t,hs in (('ext',('pmu6','ppi3','costhmu8')),('1p_ext',('Whad2','dphit3','thetap4m','thpipr4m','pp'+__import__('os').environ['PP_K']))):
            f=uproot.open(f'{G}/{g}_{t}_{m}_fte.root')
            for h in hs:
                if h+'_fte' not in f: print('MISSING prediction', g, t, m, h); bad=1
sys.exit(bad)
PYEOF
[ $miss -eq 0 ] || { echo "==== PROMOTE 050 ABORTED: inputs missing ===="; exit 1; }
echo "inputs complete $(date)"

# ---- 2. backups and retirements ------------------------------------------------------------------------
CB=configs/backup_$STAMP; UB=$PROC/univmake_backup_$STAMP; RR=$PROC/release_retired_$STAMP; DR=$REL/retired_$STAMP
mkdir -p "$CB" "$UB" "$RR" "$DR/cov"
for f in configs/adopt050/*; do b=$(basename "$f"); [ -f "configs/$b" ] && [ ! -f "$CB/$b" ] && cp -p "configs/$b" "$CB/$b"; done
echo "configs backed up to $CB ($(ls "$CB" | wc -l) files)"
while read -r PFX CFG OBS BIN; do
  case $CFG in fhc5) T=FHC5;; rhcfull) T=RHCFULL;; comb) T=COMB;; esac
  u="$PROC/${PFX}_${T}_${OBS}_univmake.root"; [ -f "$u" ] && mv "$u" "$UB/"
done < slurm/adopt050_manifest.list
for T in $T3; do
  for f in ccpi_${T}_ppi2bin_univmake.root ccpi_${T}_thetamu_univmake.root ccpi1p_${T}_dphit2bin_univmake.root; do
    [ -f "$PROC/$f" ] && mv "$PROC/$f" "$UB/"; done
  for k in xsec_${T}_ppi2bin xsec_${T}_thetamu xsec_ccpi1p_${T}_dphit2bin; do
    for f in "$PROC/$k.root" "$PROC/closure_hists_$k.root"; do [ -f "$f" ] && mv "$f" "$RR/"; done; done
  for k in incl_${T}_ppi2bin incl_${T}_thetamu 1p_${T}_dphit2bin; do
    for f in "$REL/curves_$k.tsv" "$REL/A_C_$k.tsv"; do [ -f "$f" ] && mv "$f" "$DR/"; done
    [ -d "$REL/cov/$k" ] && mv "$REL/cov/$k" "$DR/cov/"; done
done
echo "retired extractions archived: $RR, $DR; replaced universe files in $UB"

# ---- 3. install ----------------------------------------------------------------------------------------
cp configs/adopt050/* configs/
mv "$RB"/*_univmake.root "$PROC/"
for T in $T3; do cp -p "$PROC/rebuild_binreview/ccpi1p_${T}_ppK${PP_K}_univmake.root" "$PROC/ccpi1p_${T}_pp_univmake.root"; done
echo "installed $(ls configs/adopt050 | wc -l) configs and the universe files $(date)"

# ---- 4. regenerate -------------------------------------------------------------------------------------
bash regen_release.sh > ../logs/fdfix/regen_release_050.log 2>&1
tail -n 14 ../logs/fdfix/regen_release_050.log
set +u; source setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="/usr/lib64/flexiblas:$(root-config --libdir):$REPO/lib:${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
for obs in thetap thpipr; do for cfg in fhc5 rhcfull comb; do
  case $cfg in fhc5) T=FHC5;; rhcfull) T=RHCFULL;; comb) T=COMB;; esac
  W=$PROC/unfold_scratch/pa_${T}_${obs}; rm -rf "$W"; mkdir -p "$W/unfold_output"
  ln -s "$REPO/configs" "$W/configs"; ln -s "$REPO/bin" "$W/bin"; ln -s "$REPO/lib" "$W/lib"
  ( cd "$W" && "$REPO/bin/UnfolderNuMI" "configs/ccpi1p_xsec_config_numi_${obs}_${cfg}.txt" "configs/ccpi1p_${obs}_slice_config.txt" \
      "$PROC/xsec_ccpi1p_${T}_${obs}.root" ) > "$PROC/unfold_ccpi1p_${T}_${obs}.log" 2>&1 && echo "  proton angle $T $obs OK" || echo "  proton angle $T $obs FAILED"
  rm -rf "$W"
done; done
echo "==== PROMOTE 050 DONE $(date) ===="
