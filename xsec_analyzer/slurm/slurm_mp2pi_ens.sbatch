#!/bin/bash
#SBATCH --job-name=mp2piens
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=03:00:00
# SLURM array task: one member of a pseudo-experiment ensemble of the two-pion one-bin total with the MCS-calibrated
# classifier inputs. Task i = member i of ensemble ENS (nominal; alt545: GENIE multisim universe 545; altdelta: the
# Delta->N pi angular unisim). The member's eight per-run pseudo-data files are thrown from the sources of
# macros/mp2pi_ens_source.C (macros/mp2pi_ens_throw.C, one directory per member); univmake runs on them alone and their
# per-file directories are added to a copy of each configuration's calibrated universe file (macros/merge_univ_dirs.C);
# the extractions without and with the sideband constraint run with the member as pseudo-data. Output per member and
# extraction: <ENS>/m<NNN>/<cfg><mode>.log (SYSTDUMP) and the closure sidecar; the thrown files are deleted.
#   sbatch --array=0-99 --export=ALL,REPO=$PWD,ENS=nominal slurm/slurm_mp2pi_ens.sbatch
# KINVARS="pmu cmu ..." (PROD=2f): the differential variables instead of the one-bin total, from the same thrown members:
# per variable univmake with configs/mp2pi_kin_<v>_bin_config.txt, merged into a copy of mp2pi_<CFG>_total_2f_fluxw_<v>,
# extractions mp2pi_xsec_config_numi_kin2f_<v>_<cfg>[_constr] for the configurations KINCFGS (default comb); output
# <ENS>_kin/m<NNN>/<v>/: logs, closure sidecars and the unfolded vector with its total and data-statistics covariance.
set -u
REPO="${REPO:?}"; ENS="${ENS:?}"; PROD="${PROD:-mcscal}"; N=/data/uboone/processed/mp2pi_$PROD; E=$N/ens; SRC="${SRC:-$E/source}"   # SRC: other sources
# PROD: the production (mcscal, the single-model classifier; 2f, the two-fold classifier, whose one-bin total runs through
# the kin machinery: bin config mp2pi_kin_total_bin_config.txt, universes *_2f_fluxw_total, configs mp2pi_xsec_config_numi_kin2f_total_*)
if [ "$PROD" = 2f ]; then
  FPTAG=_2f_fluxw; UNIVTAG=_2f_fluxw_total; XCFG(){ echo "configs/mp2pi_xsec_config_numi_kin2f_total_${1}${2:-}.txt"; }
else
  FPTAG=_mcscal_fluxw; UNIVTAG=_mcscal_fluxw; XCFG(){ echo "configs/mp2pi_xsec_config_numi_total_${1}${2:-}_dvfrac_fluxw_mcscal${BTAG:-}.txt"; }
fi
# nominal, alt545, altdelta: the first ensembles (consecutive seeds, throw mode 1). nominal2: nominal with hashed seeds.
# Stress models outside the covariance (TTreeFormula weights): stress_band, the score-band signal +20% against the signal
# region; stress_mult, three-or-more-pion background +30% and single-pion background -20%.
declare -A ALTW=( [nominal]="" [alt545]=ens_w545 [altdelta]=ens_wdelta [nominal2]=""
  [stress_band]="=1+0.2*(CC1mu2piBDT_MC_Signal&&CC1mu2piBDT_evt_region==3)"
  [stress_mult]="=1+0.3*((!CC1mu2piBDT_MC_Signal)&&CC1mu2piBDT_mc_topology==4)-0.2*((!CC1mu2piBDT_MC_Signal)&&CC1mu2piBDT_mc_topology==2)"
  [stress_plead]="=1+0.5*CC1mu2piBDT_MC_Signal*sin(CC1mu2piBDT_mc_pi_mom[0])"
  [stress_psoft]="=1+0.5*CC1mu2piBDT_MC_Signal*cos(3.14159265*TMath::Min(Min\$(CC1mu2piBDT_mc_pi_mom),0.5)/0.5)" )
declare -A SEED=( [nominal]=20261016 [alt545]=20261116 [altdelta]=20261216 [nominal2]=20261017 [stress_band]=20261117 [stress_mult]=20261217
  [stress_plead]=20261317 [stress_psoft]=20261417 )
declare -A MODE=( [nominal]=1 [alt545]=1 [altdelta]=1 [nominal2]=0 [stress_band]=0 [stress_mult]=0 [stress_plead]=0 [stress_psoft]=0 )
# signal-shape stress models (signal only, outside the covariance): stress_plead, 1 + 0.5 sin(leading true p_pi / GeV);
# stress_psoft, 1 + 0.5 cos(pi min(softest true p_pi, 0.5 GeV) / 0.5 GeV), which moves signal between the regions.
# BINV=noband: the same members (same seeds) extracted with the reduced-transfer configuration without the score band
# (configs/mp2pi_total_bin_config_noband.txt, universes *_noband, xsec configs *_noband); output <ENS>_noband/.
declare -A TAG=( [fhc5]=FHC5 [rhcfull]=RHCFULL [comb]=COMB )
BINV="${BINV:-}"; BTAG="${BINV:+_$BINV}"; BINC="configs/mp2pi_total_bin_config${BTAG}.txt"; [ "$PROD" = 2f ] && BINC=configs/mp2pi_kin_total_bin_config.txt
KINVARS="${KINVARS:-}"; KTAG="${KINVARS:+_kin}"
[ -z "$KINVARS" ] || [ "$PROD" = 2f ] || { echo "KINVARS needs PROD=2f"; exit 1; }
M=$SLURM_ARRAY_TASK_ID; MM=$(printf "m%03d" $M); D=$E/$ENS$BTAG$KTAG/$MM; WK=$E/work/${ENS}${BTAG}${KTAG}_$MM
cd "$REPO" || exit 1
set +u; source ./setup_xsec_analyzer.sh >/dev/null 2>&1; set -u
export LD_LIBRARY_PATH="$REPO/lib:/usr/lib64/flexiblas:$(root-config --libdir):${LD_LIBRARY_PATH:-}"; export XSEC_ANALYZER_DIR="$REPO"
echo "host=$(hostname) ensemble=$ENS member=$M"
command rm -rf -- "$WK"; mkdir -p "$D" "$WK/unfold_output"; ln -s "$REPO/configs" "$WK/configs"
# 1. throw
root -l -b -q "macros/mp2pi_ens_throw.C(\"$SRC\",\"$WK/fd\",$M,${SEED[$ENS]},\"${ALTW[$ENS]}\",${MODE[$ENS]})" > "$D/throw.log" 2>&1
NF=$(ls "$WK"/fd/xsec-ana-fakedata_*.root 2>/dev/null | wc -l); [ "$NF" = 8 ] || { echo "FAIL throw: $NF files"; exit 1; }
# 2. univmake on the member's pseudo-data (the COMB list carries all eight)
FPC=configs/file_properties_numi_comb_mp2pi${FPTAG}.txt
FPM="$WK/comb_fp.txt"; sed "s|^$N/xsec-ana-fakedata_|$WK/fd/xsec-ana-fakedata_|" "$FPC" > "$FPM"
awk '$1!~/^#/ && $3=="onBNB"{print $1}' "$FPM" > "$WK/fakedata.list"
[ "$(grep -c "^$WK/fd/" "$WK/fakedata.list")" = 8 ] || { echo "FAIL file list"; exit 1; }
ok=1
if [ -n "$KINVARS" ]; then
  for v in $KINVARS; do
    mkdir -p "$D/$v"
    UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$WK/fakedata.list" configs/mp2pi_kin_${v}_bin_config.txt \
      "$WK/small_$v.root" "$FPM" > "$WK/univmake_$v.log" 2>&1
    for cfg in ${KINCFGS:-comb}; do
      FP=configs/file_properties_numi_${cfg}_mp2pi${FPTAG}.txt; FPV="$WK/${cfg}_fp.txt"
      sed "s|^$N/xsec-ana-fakedata_|$WK/fd/xsec-ana-fakedata_|" "$FP" > "$FPV"
      UV="$WK/${cfg}_${v}_univmake.root"; cp -p /data/uboone/processed/mp2pi_${TAG[$cfg]}_total${FPTAG}_${v}_univmake.root "$UV"
      root -l -b -q "macros/merge_univ_dirs.C(\"$WK/small_$v.root\",\"$UV\")" > "$WK/${cfg}_${v}_merge.log" 2>&1
      grep -q "^8 per-file directories copied under 1 root" "$WK/${cfg}_${v}_merge.log" || { echo "FAIL merge $cfg $v"; ok=0; continue; }
      for mode in "" _constr; do
        XV="$WK/${cfg}${mode}_${v}_xsec.txt"
        sed -e "s|^UnivFile .*|UnivFile $UV|" -e "s|^FPFile .*|FPFile $FPV|" configs/mp2pi_xsec_config_numi_kin2f_${v}_${cfg}${mode}.txt > "$XV"
        command rm -f -- "$WK"/unfold_output/*
        ( cd "$WK" && "$REPO/bin/UnfolderNuMI" "$XV" configs/mp2pi_kin_${v}_slice_config.txt "$D/$v/xsec_${cfg}${mode}.root" ) \
          > "$D/$v/${cfg}${mode}.log" 2>&1
        rc=$?
        for t in vec_table_unfolded_signal mat_table_cov_total mat_table_cov_DataStats; do
          cp -p "$WK/unfold_output/$t.txt" "$D/$v/${cfg}${mode}_$t.txt" 2>/dev/null || rc=9
        done
        command rm -f -- "$D/$v/xsec_${cfg}${mode}.root"
        c=$(grep -m1 '^truth: χ²' "$D/$v/${cfg}${mode}.log")
        [ "$rc" = 0 ] && [ -n "$c" ] || ok=0
        echo "$v $cfg$mode rc=$rc $c"
      done
      command rm -f -- "$UV"
    done
  done
fi
[ -z "$KINVARS" ] && UNIV_THREADS=$SLURM_CPUS_PER_TASK bin/univmake "$WK/fakedata.list" "$BINC" "$WK/small_univmake.root" "$FPM" \
  > "$WK/univmake.log" 2>&1
for cfg in $([ -z "$KINVARS" ] && echo fhc5 rhcfull comb); do
  FP=configs/file_properties_numi_${cfg}_mp2pi${FPTAG}.txt; FPV="$WK/${cfg}_fp.txt"
  sed "s|^$N/xsec-ana-fakedata_|$WK/fd/xsec-ana-fakedata_|" "$FP" > "$FPV"
  NON=$(awk '$1!~/^#/ && $3=="onBNB"' "$FPV" | grep -c "^$WK/fd/")
  UV="$WK/${cfg}_univmake.root"; cp -p /data/uboone/processed/mp2pi_${TAG[$cfg]}_total${UNIVTAG}${BTAG}_univmake.root "$UV"
  root -l -b -q "macros/merge_univ_dirs.C(\"$WK/small_univmake.root\",\"$UV\")" > "$WK/${cfg}_merge.log" 2>&1
  grep -q "^8 per-file directories copied under 1 root" "$WK/${cfg}_merge.log" || { echo "FAIL merge $cfg"; ok=0; continue; }
  for mode in "" _constr; do
    XV="$WK/${cfg}${mode}_xsec.txt"
    sed -e "s|^UnivFile .*|UnivFile $UV|" -e "s|^FPFile .*|FPFile $FPV|" "$(XCFG $cfg $mode)" > "$XV"
    ( cd "$WK" && "$REPO/bin/UnfolderNuMI" "$XV" configs/mp2pi_total_slice_config.txt "$D/xsec_${cfg}${mode}.root" ) > "$D/${cfg}${mode}.log" 2>&1
    rc=$?; s=$(grep 'SYSTDUMP] sigma_int' "$D/${cfg}${mode}.log" | awk '{print $3}')
    [ "$rc" = 0 ] && [ -n "$s" ] || ok=0
    echo "$cfg$mode rc=$rc sigma_int=$s onBNB=$NON"
  done
  command rm -f -- "$UV"
done
grep -h "events$" "$D/throw.log" | tail -8
command rm -rf -- "$WK"
[ $ok = 1 ] && echo "DONE $ENS $M" || { echo "FAIL $ENS $M"; exit 1; }
