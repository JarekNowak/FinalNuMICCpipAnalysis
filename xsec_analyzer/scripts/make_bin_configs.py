#!/usr/bin/env python3
"""make_bin_configs.py -- write UniverseMaker bin configs, slice configs and study xsec configs
for 1D and 2D binning candidates (0.50-criterion study, 2026-09-16).

Conventions copied from the hand-written candidate configs (e.g. ccpi_costhpi6_*_opt.txt):
  * true bins  : "0 0 \"<S>_MC_Signal && <lo <= t < hi>\"" ; last bin open when open_top
  * reco bins  : "0 0 \"<S>_Selected && <lo <= r < hi>\""  ; first bin has NO lower edge and
                 the last bin none when open_top (as in the released configs)
  * one background bin "1 -1 \"!<S>_MC_Signal\""
2D: analysis bins are flattened Y slice by Y slice, so each slice is a contiguous block (the
per-slice A_C restriction in UnfolderNuMI needs that). The slice config has one slice per Y bin,
active variable X, other variable Y fixed to the slice range (UnfolderNuMI divides by the Y width
through other_vars), plus the "bin number" slice over everything.
Study xsec configs carry only the uBTune and FakeData predictions (generator files follow the
released binning only), like the 2026-09-05 candidates.

usage: scripts/make_bin_configs.py [module]   (writes every candidate in <module>.CANDS; default c50_candidates)
"""
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(HERE, 'configs', 'c50')  # default study subdirectory, kept out of norm_manifest.py's release glob
FLUX = {'fhc5': '6.81159e-10', 'rhcfull': '6.446460e-10', 'comb': '6.60865e-10'}
SYST = {'ccpi': 'configs/ccpi_systcalc_numi.conf', 'ccpi1p': 'configs/ccpi1p_systcalc_numi.conf'}
FPM = {('ccpi', 'fhc5'): 'configs/file_properties_numi_fhc5.txt',
       ('ccpi', 'rhcfull'): 'configs/file_properties_numi_rhcfull.txt',
       ('ccpi1p', 'fhc5'): 'configs/file_properties_numi_fhc5_w.txt',
       ('ccpi1p', 'rhcfull'): 'configs/file_properties_numi_rhcfull_w.txt',
       ('ccpi', 'comb'): 'configs/file_properties_numi_comb.txt',
       ('ccpi1p', 'comb'): 'configs/file_properties_numi_comb_w.txt'}
TAG = {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}
GENROOT = '/home/t2k/nowak/MICRO_GEN/gen2d'.replace('MICRO_GEN', 'MicroBooNE/working_xsec_analyzer/generator_predictions')
GENS = [('GENIE', 'GENIE', 'genie'), ('GiBUU', 'GiBUU', 'gibuu'), ('NEUT', 'NEUT', 'neut'), ('NuWro', 'NuWro', 'nuwro')]
GMODE = {'fhc5': 'fhc', 'rhcfull': 'rhc', 'comb': 'comb'}
PPR = ('sqrt(pow(sqrt(pow({S}_candidate_pion_mom_reco,2)+0.011164)-0.10566+0.13957,2)-0.019480)')

# variable: (true expr, reco expr, root label, units) ; {S} is the selection prefix
VARS = {
    'pmu':     ('{S}_candidate_muon_mom_true', '{S}_candidate_muon_mom_reco', 'p_{#mu}', ' (GeV/c)'),
    'costhmu': ('{S}_candidate_muon_costh_true', '{S}_candidate_muon_costh_reco', 'cos#theta_{#mu}', ''),
    'thetamu': ('TMath::ACos({S}_candidate_muon_costh_true)', 'TMath::ACos({S}_candidate_muon_costh_reco)', '#theta_{#mu}', ' (rad)'),
    'costhpi': ('{S}_candidate_pion_costh_true', '{S}_candidate_pion_costh_reco', 'cos#theta_{#pi}', ''),
    'thmupi':  ('{S}_true_mu_pi_opening_angle', '{S}_mu_pi_opening_angle', '#theta_{#mu#pi}', ' (rad)'),
    'ppi':     ('{S}_candidate_pion_mom_true', PPR, 'p_{#pi}', ' (GeV/c)'),
    'dpt':     ('{S}_deltaPt_true', '{S}_deltaPt_reco', '#deltap_{T}', ' (GeV/c)'),
    'dalphat': ('{S}_deltaAlphaT_true', '{S}_deltaAlphaT_reco', '#delta#alpha_{T}', ' (deg)'),
    'dphit':   ('{S}_deltaPhiT_true', '{S}_deltaPhiT_reco', '#delta#phi_{T}', ' (deg)'),
    'pn':      ('{S}_pn_true', '{S}_pn_reco', 'p_{n}', ' (GeV/c)'),
    'Wpipr':   ('{S}_W_pipr_true', '{S}_W_pipr_reco', 'W_{#pip}', ' (GeV/c^{2})'),
    'Whad':    ('{S}_W_had_true', '{S}_W_had_reco', 'W_{had}', ' (GeV/c^{2})'),
    'thetap':  ('TMath::ACos({S}_proton_costh_true)', 'TMath::ACos({S}_proton_costh_reco)', '#theta_{p}', ' (rad)'),
    'costhp':  ('{S}_proton_costh_true', '{S}_proton_costh_reco', 'cos#theta_{p}', ''),
    'thpipr':  ('{S}_pi_pr_opening_angle_true', '{S}_pi_pr_opening_angle_reco', '#theta_{#pip}', ' (rad)'),
}


def rng(expr, lo, hi, lo_open, hi_open):
    parts = []
    if not lo_open: parts.append(f'{expr} >= {lo:.3f}')
    if not hi_open: parts.append(f'{expr} < {hi:.3f}')
    return ' && '.join(parts)


def write(name, pfx, xvar, xedges, xopen, yvar=None, yedges=None, yopen=False, note='',
          xlow_open=False, ylow_open=False, subdir='c50', univdir='/data/uboone/processed/rebuild_c50',
          cfgs=('fhc5', 'rhcfull'), gen_obs=None, suffix='_c50'):
    """xlow_open/ylow_open: the first TRUE bin has no lower edge (released proton-tagged convention).
    subdir '' writes the live configs/ family (bin/slice without suffix); gen_obs names the <obs>_fte
    histogram in generator_predictions/gen2d/<gen>[_1p]_ext_<mode>_fte.root (None: no generator lines)."""
    """xedges: list (shared) or list of lists (one per Y slice)."""
    S = 'CC1mu1pi1p' if pfx == 'ccpi1p' else 'CC1mu1piXp'
    xt, xr, xl, xu = [v.format(S=S) if i < 2 else v for i, v in enumerate(VARS[xvar])]
    twoD = yvar is not None
    nsl = len(yedges) - 1 if twoD else 1
    xe = [xedges] * nsl if not isinstance(xedges[0], (list, tuple)) else xedges
    true_lines, reco_lines, slices = [], [], []
    idx = 0
    for j in range(nsl):
        e = xe[j]; nb = len(e) - 1; first = idx
        for b in range(nb):
            last = b == nb - 1
            tx = rng(xt, e[b], e[b + 1], xlow_open and b == 0, xopen and last)
            rx = rng(xr, e[b], e[b + 1], b == 0, xopen and last)
            if twoD:
                yt, yr = VARS[yvar][0].format(S=S), VARS[yvar][1].format(S=S)
                ylast = j == nsl - 1
                tx = rng(yt, yedges[j], yedges[j + 1], ylow_open and j == 0, yopen and ylast) + ' && ' + tx
                rx = rng(yr, yedges[j], yedges[j + 1], j == 0, yopen and ylast) + ' && ' + rx
            true_lines.append(f'0 0 "{S}_MC_Signal && {tx}"')
            reco_lines.append(f'0 0 "{S}_Selected && {rx}"')
            idx += 1
        slices.append((j, e, first))
    nbins = idx
    tag = f'ccpi_{S}_{name}_{"2D" if twoD else "1D"}'
    bin_txt = '\n'.join([tag, 'stv_tree', S, str(nbins + 1)] + true_lines + [f'1 -1 "!{S}_MC_Signal"', str(nbins)] + reco_lines) + '\n'

    # slice config
    out = []
    var_lines = [f'"{xl}" "{xu}" "{xl}" "{xu}"']
    if twoD:
        yl, yu = VARS[yvar][2], VARS[yvar][3]
        var_lines.append(f'"{yl}" "{yu}" "{yl}" "{yu}"')
    var_lines.append('"bin number" "" "bin number" ""')
    out += [str(len(var_lines))] + var_lines
    out.append(str(nsl + 1))
    for j, e, first in slices:
        out += ['"events"', '1', f'0 {len(e)} ' + ' '.join(f'{v:.3f}' for v in e)]
        if twoD:
            # an open top Y slice still divides by its nominal width, as the 1D open-top
            # bins do (SliceBinning uses the last listed edge)
            out.append(f'1 1 {yedges[j]:.3f} {yedges[j + 1]:.3f}')
        else:
            out.append('0')
        nb = len(e) - 1
        out.append(str(nb))
        out += [f'{first + b} 1 {b + 1}' for b in range(nb)]
    bn_idx = 2 if twoD else 1
    out += ['"events"', '1', f'{bn_idx} {nbins + 1} ' + ' '.join(str(i) for i in range(nbins + 1)), '0', str(nbins)]
    out += [f'{i} 1 {i + 1}' for i in range(nbins)]
    slice_txt = '\n'.join(out) + '\n'

    d = os.path.join(os.path.dirname(CFG), subdir) if subdir else os.path.dirname(CFG)
    os.makedirs(d, exist_ok=True)
    rel = 'configs/' + (subdir + '/' if subdir else '')
    open(os.path.join(d, f'{pfx}_{name}_bin_config{suffix}.txt'), 'w').write(bin_txt)
    open(os.path.join(d, f'{pfx}_{name}_slice_config{suffix}.txt'), 'w').write(slice_txt)
    for cfg in cfgs:
        T = TAG[cfg]
        x = (f'# {note}\n'
             f'UnivFile {univdir}/{pfx}_{T}_{name}_univmake.root\n'
             f'SystFile {SYST[pfx]}\nFPFile {FPM[(pfx, cfg)]}\nUnfold WienerSVD 1 second-deriv\nFlux {FLUX[cfg]}\n'
             'Prediction uBTune "MicroBooNE Tune" univ CV\nPrediction FakeData "Fakedata" univ FakeData\n')
        if gen_obs:
            fam = '1p_ext' if pfx == 'ccpi1p' else 'ext'
            for key, label, g in GENS:
                x += f'Prediction {key} "{label}" file {GENROOT}/{g}_{fam}_{GMODE[cfg]}_fte.root {gen_obs}_fte\n'
        open(os.path.join(d, f'{pfx}_xsec_config_numi_{name}{suffix}_{cfg}.txt'), 'w').write(x)
    print(f'wrote {pfx} {name}: {nbins} bins' + (f' ({nsl} Y slices)' if twoD else ''))


if __name__ == '__main__':
    import sys
    os.makedirs(CFG, exist_ok=True)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    mod = sys.argv[1] if len(sys.argv) > 1 else 'c50_candidates'
    CANDS = __import__(mod).CANDS   # list of kwargs dicts, kept next to this script
    for c in CANDS:
        write(**c)
