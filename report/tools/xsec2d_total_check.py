"""xsec2d_total_check.py -- normalisation check of the two-dimensional extractions.

Summing the flattened sigma per analysis bin of a 2D extraction gives its integrated cross section.
That is not expected to equal the released one-bin total extraction: both are A_C-smeared and the two
A_C are different, so each carries its own loss of normalisation (the row sums are not one). What IS
expected is that the fake data and the fake-data truth lose the SAME amount, so the double ratio

    [ sigma_2D(data) / sigma_2D(truth) ] / [ sigma_total(data) / sigma_total(truth) ]

is one. A gross normalisation error in the 2D chain -- a wrong flux, a double-counted file, a bin
area applied twice -- would move the data and leave the truth, and show up here.

The one-bin total is a differential over cos(theta_mu) in [-1,1], so its histogram value must be
multiplied by that width of 2 to become a cross section.
    python3 report/tools/xsec2d_total_check.py
"""
import os, sys
import numpy as np, uproot
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'xsec_analyzer', 'scripts'))
from newobs_candidates import CANDS
RB, LIVE = '/data/uboone/processed/rebuild_2d', '/data/uboone/processed'
CFG = ['FHC5', 'RHCFULL', 'COMB']

tot = {}
for cfg in CFG:
    for pfx in ('ccpi', 'ccpi1p'):
        f = uproot.open(f'{LIVE}/closure_hists_xsec_{pfx + "_" if pfx == "ccpi1p" else ""}{cfg}_total.root')
        h = f['h_unfolded_nuwro']; w = np.diff(h.axis().edges())
        tot[(pfx, cfg)] = ((h.values() * w).sum(), (h.errors() * w).sum(), (f['h_fakedata_truth'].values() * w).sum())

print('%-24s %-8s %8s %8s %8s %8s %9s' % ('pair', 'cfg', 'sig_2D', 'sig_tot', '2D d/t', 'tot d/t', 'double'))
dev = []
for c in CANDS:
    if not c.get('yvar'): continue
    for cfg in CFG:
        p = f"{RB}/closure_hists_all_xsec_{c['pfx']}_{cfg}_{c['name']}.root"
        if not os.path.exists(p): continue
        f = uproot.open(p)
        s, t = f['h_unfolded_nuwro'].values().sum(), f['h_fakedata_truth'].values().sum()
        st, et, tt = tot[(c['pfx'], cfg)]
        d = (s / t) / (st / tt); dev.append((abs(d - 1), c['name'], cfg, d))
        print('%-24s %-8s %8.3f %8.3f %8.3f %8.3f %9.3f' % (c['name'], cfg, s, st, s / t, st / tt, d))
dev.sort()
print('\nmedian |double ratio - 1| = %.3f over %d extractions; the one-bin total itself carries '
      '%.0f--%.0f%% uncertainty.' % (np.median([d[0] for d in dev]), len(dev),
                                     100 * min(e / s for s, e, _ in tot.values()),
                                     100 * max(e / s for s, e, _ in tot.values())))
print('largest deviations:')
for d in dev[-4:]: print('   %-22s %-8s %.3f' % (d[1], d[2], d[3]))
