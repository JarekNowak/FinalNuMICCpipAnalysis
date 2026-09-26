#!/usr/bin/env python3
"""adopt_050_configs.py -- stage the binnings selected under the 0.50 migration criterion (user decision
2026-09-26, after the binning review of the technical supplement, sec:binreview).

Writes, into configs/adopt050/ and under their final live names, the bin, slice and cross-section
configs of every observable whose binning changes. Nothing live is touched here: promote_050.sh backs
up the files these replace and installs them together with the rebuilt universe files.

  inclusive      p_mu 6 bins (same names), p_pi 3 bins (new names *_3bin, ppi3bin), cos theta_mu 8 bins
  proton-tagged  p_mu and cos theta_mu at the inclusive edges, W_had 2 bins, delta phi_T 3 bins
                 (new names *_3bin, dphit3bin), theta_p and theta_pip 4 bins at their maximin edges
  unchanged      cos theta_pi (5 bins, user decision), theta_mupi, the two-bin delta p_T, delta alpha_T
                 and p_n, the proton-tagged p_pi (the 3-bin edges fail 0.50 there), the withdrawn W_pipr
  dropped        theta_mu (user decision: cos theta_mu carries the muon angle)

Bin and slice configs are the review candidates (configs/binreview/, written by make_bin_configs.py with
the conventions of the released configs) with the release name on the first line; the proton-tagged
p_mu and cos theta_mu are the inclusive ones with the selection prefix swapped, exactly as the released
files are. Generator predictions come from gen2d (obs_ext.h): pmu6, ppi3, costhmu8 in
<gen>_ext_{fhc,rhc,comb}_fte.root, Whad2, dphit3, thetap4m, thpipr4m in <gen>_1p_ext_*_fte.root.

    python3 scripts/adopt_050_configs.py [PP_K]

PP_K (added 2026-09-26, user decision): the number of bins chosen for the new proton-tagged observable p_p,
the proton momentum, from its review candidates ppK<PP_K>; the configs are staged as ccpi1p_pp_* with
generator predictions pp<PP_K>_fte in all three configurations.
"""
import os, re, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = os.path.join(HERE, 'configs'); BR = os.path.join(C, 'binreview'); OUT = os.path.join(C, 'adopt050')
GP = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d'
CFGS = ('fhc5', 'rhcfull', 'comb'); TAG = {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}
MODE = {'fhc5': 'fhc', 'rhcfull': 'rhc', 'comb': 'comb'}
GEN = {'GENIE': 'genie', 'GiBUU': 'gibuu', 'NEUT': 'neut', 'NuWro': 'nuwro'}


def read(p): return open(p).read().split('\n')


def put(name, lines):
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, name), 'w').write('\n'.join(lines))
    print('staged', name)


def renamed(src, first, swap=False):
    L = read(os.path.join(BR, src))
    if swap: L = [l.replace('CC1mu1piXp', 'CC1mu1pi1p') for l in L]
    L[0] = first
    return L


def xsec(src, univ=None, hist=None, sample='incl', cfg='fhc5'):
    """copy a live cross-section config, pointing UnivFile at a new universe file name and the four
    generator predictions at the gen2d file of this configuration with histogram HIST"""
    L = read(os.path.join(C, src)); out = []
    for l in L:
        if univ and l.startswith('UnivFile '):
            l = 'UnivFile /data/uboone/processed/' + univ
        m = re.match(r'Prediction (\w+) "([^"]+)" file (\S+) (\S+)', l)
        if m and hist and m.group(1) in GEN:
            f = f"{GP}/{GEN[m.group(1)]}_{'1p_ext' if sample == '1p' else 'ext'}_{MODE[cfg]}_fte.root"
            l = f'Prediction {m.group(1)} "{m.group(2)}" file {f} {hist}'
        out.append(l)
    return out


# ---- bin and slice configs ---------------------------------------------------------------------
put('ccpi_pmu_bin_config_opt.txt', renamed('ccpi_pmuK6_bin_config_br.txt', 'ccpi_CC1mu1piXp_pmu_1D'))
put('ccpi_pmu_slice_config_opt.txt', read(os.path.join(BR, 'ccpi_pmuK6_slice_config_br.txt')))
put('ccpi_ppi_bin_config_3bin.txt', renamed('ccpi_ppiK3_bin_config_br.txt', 'ccpi_CC1mu1piXp_ppi_1D'))
put('ccpi_ppi_slice_config_3bin.txt', read(os.path.join(BR, 'ccpi_ppiK3_slice_config_br.txt')))
put('ccpi_costhmu_bin_config_opt.txt', renamed('ccpi_costhmuK8_bin_config_br.txt', 'ccpi_CC1mu1piXp_costhmu_1D'))
put('ccpi_costhmu_slice_config_opt.txt', read(os.path.join(BR, 'ccpi_costhmuK8_slice_config_br.txt')))
put('ccpi1p_pmu_bin_config.txt', renamed('ccpi_pmuK6_bin_config_br.txt', 'ccpi_CC1mu1pi1p_pmu_1D', swap=True))
put('ccpi1p_pmu_slice_config.txt', read(os.path.join(BR, 'ccpi_pmuK6_slice_config_br.txt')))
put('ccpi1p_costhmu_bin_config.txt', renamed('ccpi_costhmuK8_bin_config_br.txt', 'ccpi_CC1mu1pi1p_costhmu_1D', swap=True))
put('ccpi1p_costhmu_slice_config.txt', read(os.path.join(BR, 'ccpi_costhmuK8_slice_config_br.txt')))
put('ccpi1p_Whad_bin_config.txt', renamed('ccpi1p_WhadK2_bin_config_br.txt', 'ccpi_CC1mu1pi1p_Whad_1D'))
put('ccpi1p_Whad_slice_config.txt', read(os.path.join(BR, 'ccpi1p_WhadK2_slice_config_br.txt')))
put('ccpi1p_dphit_bin_config_3bin.txt', renamed('ccpi1p_dphitK3_bin_config_br.txt', 'ccpi_CC1mu1pi1p_dphit_1D'))
put('ccpi1p_dphit_slice_config_3bin.txt', read(os.path.join(BR, 'ccpi1p_dphitK3_slice_config_br.txt')))
put('ccpi1p_thetap_bin_config.txt', renamed('ccpi1p_thetapK4_bin_config_br.txt', 'ccpi_CC1mu1pi1p_thetap_1D'))
put('ccpi1p_thetap_slice_config.txt', read(os.path.join(BR, 'ccpi1p_thetapK4_slice_config_br.txt')))
put('ccpi1p_thpipr_bin_config.txt', renamed('ccpi1p_thpiprK4_bin_config_br.txt', 'ccpi_CC1mu1pi1p_thpipr_1D'))
put('ccpi1p_thpipr_slice_config.txt', read(os.path.join(BR, 'ccpi1p_thpiprK4_slice_config_br.txt')))

# ---- cross-section configs ---------------------------------------------------------------------
for c in CFGS:
    T = TAG[c]
    put(f'ccpi_xsec_config_numi_pmu_{c}.txt', xsec(f'ccpi_xsec_config_numi_pmu_{c}.txt', hist='pmu6_fte', cfg=c))
    put(f'ccpi_xsec_config_numi_costhmu_{c}.txt', xsec(f'ccpi_xsec_config_numi_costhmu_{c}.txt', hist='costhmu8_fte', cfg=c))
    put(f'ccpi_xsec_config_numi_ppi3bin_{c}.txt',
        xsec(f'ccpi_xsec_config_numi_ppi2bin_{c}.txt', univ=f'ccpi_{T}_ppi3bin_univmake.root', hist='ppi3_fte', cfg=c))
    # proton-tagged: generator curves exist for FHC W/TKI and for theta_p / theta_pip in all three
    put(f'ccpi1p_xsec_config_numi_Whad_{c}.txt', xsec(f'ccpi1p_xsec_config_numi_Whad_{c}.txt', hist='Whad2_fte', sample='1p', cfg=c))
    put(f'ccpi1p_xsec_config_numi_dphit3bin_{c}.txt',
        xsec(f'ccpi1p_xsec_config_numi_dphit2bin_{c}.txt', univ=f'ccpi1p_{T}_dphit3bin_univmake.root', hist='dphit3_fte', sample='1p', cfg=c))
    put(f'ccpi1p_xsec_config_numi_thetap_{c}.txt', xsec(f'ccpi1p_xsec_config_numi_thetap_{c}.txt', hist='thetap4m_fte', sample='1p', cfg=c))
    put(f'ccpi1p_xsec_config_numi_thpipr_{c}.txt', xsec(f'ccpi1p_xsec_config_numi_thpipr_{c}.txt', hist='thpipr4m_fte', sample='1p', cfg=c))

# ---- proton momentum (new observable, 2026-09-26) ----------------------------------------------
if len(sys.argv) > 1:
    K = int(sys.argv[1])
    put('ccpi1p_pp_bin_config.txt', renamed(f'ccpi1p_ppK{K}_bin_config_br.txt', 'ccpi_CC1mu1pi1p_pp_1D'))
    put('ccpi1p_pp_slice_config.txt', read(os.path.join(BR, f'ccpi1p_ppK{K}_slice_config_br.txt')))
    for c in CFGS:
        L = read(os.path.join(BR, f'ccpi1p_xsec_config_numi_ppK{K}_br_{c}.txt'))
        L = ['# CC1mu1pi1p (proton-tagged) proton momentum p_p (leading proton, 0.3 GeV/c threshold), ' + c + '.',
             f'# Binning: maximin {K} bins under the 0.50 migration criterion (binning review, supplement sec:binreview).'] + \
            [(f'UnivFile /data/uboone/processed/ccpi1p_{TAG[c]}_pp_univmake.root' if l.startswith('UnivFile ') else l) for l in L if l.strip()]
        f = f"{GP}/{{}}_1p_ext_{MODE[c]}_fte.root"
        L += [f'Prediction {n} "{n}" file {f.format(GEN[n])} pp{K}_fte' for n in ('GENIE', 'GiBUU', 'NEUT', 'NuWro')]
        put(f'ccpi1p_xsec_config_numi_pp_{c}.txt', L)
