"""binreview_eval.py -- metrics of every binning of the 0.68-versus-0.50 review, released and candidate,
in FHC, RHC and combined: bins, worst migration diagonal (the unfolder's [DIAGDUMP]), A_C row-sum range,
A_C conditioning s_min/s_1, the Wiener filter factors (eigenvalues of A_C), their sum (degrees of freedom
for signal) and the number >= 0.5 and >= 0.2, bin-averaged data-statistical and total uncertainty (SYSTDUMP), per-bin total
uncertainty range, closure chi2/ndf and p, and the unfolded/truth ratio range. Writes
logs/binreview/eval.tsv.
    python3 report/tools/binreview_eval.py
"""
import os, sys, re, csv, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from binning_eval_c50 import dump, side
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'xsec_analyzer', 'scripts')))
import binreview_candidates as BC
P = '/data/uboone/processed/'; RB = P + 'rebuild_binreview/'
LG = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
OUT = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/binreview/eval.tsv'
CF = {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}
REL = {'incl': {'pmu': 'pmu', 'ppi': 'ppi2bin', 'costhmu': 'costhmu', 'thetamu': 'thetamu', 'costhpi': 'costhpi', 'thmupi': 'thmupi'},
       '1p': {'Whad': 'Whad', 'Wpipr': 'Wpipr', 'dpt': 'dpt2bin', 'dalphat': 'dalphat2bin', 'dphit': 'dphit2bin', 'pn': 'pn2bin',
              'pmu': 'pmu', 'ppi': 'ppi2bin', 'costhmu': 'costhmu', 'costhpi': 'costhpi', 'thmupi': 'thmupi'}}

def diag(log):
    for l in open(log):
        if '[DIAGDUMP]' in l: return [float(x) for x in l.split(':')[-1].split()]
    return None

import numpy as np, uproot
def acmat(side_file):
    return uproot.open(side_file)['h_A_C'].values().T     # values() is [x=j, y=i]; transpose gives A_C[i][j]

def svals(side_file):
    s = np.linalg.svd(acmat(side_file), compute_uv=False); return s / s[0]

def row(fam, obs, label, cfg, log, side_file):
    if not (os.path.exists(log) and os.path.exists(side_file)): return None
    d, chi = dump(log); nb, rows, cond, rel, clo = side(side_file); dg = diag(log); sv = svals(side_file)
    A = acmat(side_file)
    # A_C is similar to the diagonal matrix of Wiener filter factors W_k = S_k^2/(S_k^2+N_k^2): its eigenvalues
    # are those factors (real, in [0,1)) and its trace their sum, the degrees of freedom for signal. A factor
    # >= 0.5 is a component the data determine with signal-to-noise >= 1.
    ev = np.sort(np.linalg.eigvals(A).real)[::-1]
    return dict(family=fam, observable=obs, binning=label, config=cfg, bins=nb, dfs=float(np.trace(A)),
                n50=int((ev >= 0.5).sum()), n20=int((ev >= 0.2).sum()), eig=' '.join(f'{x:.3f}' for x in ev),
                min_diag=min(dg) if dg else float('nan'), diags=' '.join(f'{x:.3f}' for x in dg) if dg else '',
                ac_row_min=min(rows), ac_row_max=max(rows), cond=cond,
                rank10=int((sv >= 0.1).sum()), rank30=int((sv >= 0.3).sum()), svals=' '.join(f'{x:.3f}' for x in sv),
                datastat_pct=d.get('DataStats', float('nan')), total_pct=d.get('total', float('nan')),
                pred_pct=d.get('PredTotal', float('nan')), bin_unc_min=min(rel), bin_unc_max=max(rel),
                chi2=chi[0] if chi else float('nan'), ndf=chi[1] if chi else 0, p=chi[2] if chi else float('nan'),
                clo_min=min(clo) if clo else float('nan'), clo_max=max(clo) if clo else float('nan'))

rows = []
for fam, m in REL.items():
    for obs, rob in m.items():
        for cfg in CF:
            if fam == 'incl':
                lg = LG + f"{cfg}_{'ppi' if rob == 'ppi2bin' else rob}.raw"; sd = P + f'closure_hists_xsec_{CF[cfg]}_{rob}.root'
            else:
                lg = LG + f'1p_{cfg}_{rob}.raw'; sd = P + f'closure_hists_xsec_ccpi1p_{CF[cfg]}_{rob}.root'
            r = row(fam, obs, 'released', cfg, lg, sd)
            if r: rows.append(r)
# proton-angle products (proposed secondary; built 2026-09-17/20 outside the indexed release)
for obs in ('thetap', 'thpipr'):
    for cfg in CF:
        r = row('1p', obs, 'released', cfg, P + f'unfold_ccpi1p_{CF[cfg]}_{obs}.log', P + f'closure_hists_xsec_ccpi1p_{CF[cfg]}_{obs}.root')
        if r: r['note'] = 'study extraction (not in the release index)'; rows.append(r)
for c in BC.CANDS:
    fam = '1p' if c['pfx'] == 'ccpi1p' else 'incl'
    for cfg in c['cfgs']:
        base = f"{c['pfx']}_{CF[cfg]}_{c['name']}"
        r = row(fam, c['xvar'], c['name'], cfg, RB + f'unfold_{base}.log', RB + f'closure_hists_xsec_{base}.root')
        if r: r['note'] = c.get('note', ''); r['edges'] = ' '.join(str(e) for e in c['xedges']); rows.append(r)
        else: print('pending', base)
keys = ['family', 'observable', 'binning', 'config', 'bins', 'dfs', 'n50', 'n20', 'min_diag', 'ac_row_min', 'ac_row_max', 'cond', 'rank10', 'rank30', 'datastat_pct', 'pred_pct', 'total_pct',
        'bin_unc_min', 'bin_unc_max', 'chi2', 'ndf', 'p', 'clo_min', 'clo_max', 'diags', 'eig', 'svals', 'edges', 'note']
with open(OUT, 'w') as o:
    w = csv.DictWriter(o, keys, delimiter='\t', extrasaction='ignore'); w.writeheader()
    for r in rows: w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
print(f'wrote {OUT}: {len(rows)} rows')
