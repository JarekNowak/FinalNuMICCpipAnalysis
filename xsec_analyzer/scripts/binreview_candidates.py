# Candidates of the binning review of 2026-09-25 (0.68 versus 0.50 migration-diagonal criteria, on the
# corrected exposure): the maximin binnings of report/tools/maximin_binning.py (floor 70 selected FHC
# signal events per bin, the minimum of the released schemes; worst diagonal over FHC and RHC), for
# every bin count from one above the released scheme to the finest passing 0.50, and the coarser
# p_mu schemes that pass the criteria. Built for FHC, RHC and combined into rebuild_binreview/.
import json, os
_J = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/binreview/scan_floor{}_{}.json'
_W = dict(subdir='binreview', univdir='/data/uboone/processed/rebuild_binreview', cfgs=('fhc5', 'rhcfull', 'comb'), suffix='_br')
OPEN = {'pmu': True, 'ppi': True, 'costhmu': False, 'thetamu': False, 'costhpi': False, 'thmupi': False,
        'Whad': True, 'Wpipr': True, 'dpt': True, 'dalphat': False, 'dphit': False, 'pn': True, 'thetap': False, 'thpipr': False,
        'pp': True}
# K=5 for the angles: maximin edges at the released bin count, to separate the effect of the edges from the
# effect of the number of bins
KS = {'incl': {'pmu': [4, 5, 6, 7], 'ppi': [3], 'costhmu': [5, 6, 7, 8, 9], 'thetamu': [5, 6, 7, 8, 9, 10],
               'costhpi': [5, 6, 7, 8, 9, 10], 'thmupi': [5, 6, 7, 8, 9, 10, 11]},
      # proton-tagged: floor 50 (the family is 2.6 times smaller); W_had K=3 only exists at floor 30
      '1p': {'Whad': [2, (3, 30)], 'Wpipr': [2], 'dpt': [3], 'dalphat': [2, 3], 'dphit': [3], 'pn': [3],
             'thetap': [3, 4, 5, 6], 'thpipr': [3, 4, 5, 6],
             # proton momentum, added 2026-09-26 (user decision): every bin count up to the finest passing 0.50
             'pp': [2, 3, 4, 5]}}
FLOOR = {'incl': 70, '1p': 50}
CANDS = []
for fam, pfx in (('incl', 'ccpi'), ('1p', 'ccpi1p')):
    for obs, ks in KS[fam].items():
        for K in ks:
            K, fl = (K if isinstance(K, tuple) else (K, FLOOR[fam]))
            f = _J.format(fl, fam)
            if not os.path.exists(f): continue
            v = json.load(open(f))[f'{fam}:{obs}'].get(str(K))
            if not v: continue
            worst, edges = v
            name = f'{obs}K{K}' + ('' if fl == FLOOR[fam] else f'f{fl}')
            CANDS.append(dict(name=name, pfx=pfx, xvar=obs, xedges=edges, xopen=OPEN[obs],
                              note=f'maximin K={K}, worst diag FHC/RHC {worst:.3f} (floor {fl})', **_W))

# The muon angle in cos theta_mu: a partition has the same migration diagonals in either parameterisation,
# and the theta_mu scan grid is finer in the forward region than the cos grid (a cos step of 0.025 is
# 0.22 rad at the first edge), so the best theta_mu partitions, written in cos, are cos theta_mu
# candidates too: K=8 passes 0.68 and K=10 passes 0.50, one bin more than the cos-grid scan allows.
import math
for K in (8, 10):
    v = json.load(open(_J.format(70, 'incl')))['incl:thetamu'].get(str(K))
    if not v: continue
    worst, th = v
    ce = sorted(round(math.cos(t), 4) for t in th); ce[0], ce[-1] = -1.0, 1.0
    CANDS.append(dict(name=f'costhmuT{K}', pfx='ccpi', xvar='costhmu', xedges=ce, xopen=False,
                      note=f'theta_mu maximin K={K} partition in cos, worst diag FHC/RHC {worst:.3f} (floor 70)', **_W))
