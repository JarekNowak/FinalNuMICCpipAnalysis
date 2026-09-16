#!/usr/bin/env python3
"""check_ext_bins.py -- the edges of the 2026-09-16 observables live in two places: the config
candidates (scripts/newobs_candidates.py) and the generator definitions
(generator_predictions/gen2d/obs_ext.h). A mismatch would silently compare the data with a
prediction in a different binning, so check both edges and open-edge flags, and (when the FTE files
exist) that every generator histogram has as many bins as the config has true signal bins."""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from newobs_candidates import CANDS, PI
H = open(os.path.join(REPO, '..', 'generator_predictions', 'gen2d', 'obs_ext.h')).read()
num = lambda s: [PI if t.strip() == 'PI' else float(t) for t in s.split(',')]
bad = 0
def axis(txt):
    m = re.match(r'\{\s*\{([^}]*)\}\s*,\s*(true|false)\s*,\s*(true|false)\s*\}', txt.strip())
    return num(m.group(1)), m.group(2) == 'true', m.group(3) == 'true'
for c in CANDS:
    if c.get('yvar'):
        AX = r'(\{\s*\{[^}]*\},\s*\w+,\s*\w+\s*\})'
        m = re.search(r'Pair \w+ = \{\s*"%s",\s*%s,\s*%s\s*\}' % (c['name'], AX, AX), H)
        if not m: print('MISSING in obs_ext.h:', c['name']); bad += 1; continue
        gx, gy = axis(m.group(1)), axis(m.group(2))
        want = [(c['xedges'], c.get('xlow_open', False), c['xopen']), (c['yedges'], c.get('ylow_open', False), c.get('yopen', False))]
        got = [gx, gy]
    else:
        m = re.search(r'Axis %s_axis\(\) \{ Axis a = (\{\s*\{[^}]*\},\s*\w+,\s*\w+\s*\})' % c['name'], H)
        if not m: print('MISSING in obs_ext.h:', c['name']); bad += 1; continue
        want = [(c['xedges'], c.get('xlow_open', False), c['xopen'])]; got = [axis(m.group(1))]
    for w, g in zip(want, got):
        if [round(v, 6) for v in w[0]] != [round(v, 6) for v in g[0]] or w[1:] != g[1:]:
            print(f"EDGE MISMATCH {c['name']}: config {w} vs generator {g}"); bad += 1
    # FTE bin counts
    try:
        import ROOT
    except ImportError:
        continue
    sub = c['subdir']; d = os.path.join(REPO, 'configs', sub)
    ntrue = sum(1 for l in open(os.path.join(d, f"{c['pfx']}_{c['name']}_bin_config{c['suffix']}.txt")) if l.startswith('0 0 "') and '_MC_Signal &&' in l)
    fam = '1p_ext' if c['pfx'] == 'ccpi1p' else 'ext'
    for g in ('genie', 'gibuu', 'neut', 'nuwro'):
        for mode in ('fhc', 'rhc', 'comb'):
            p = os.path.join(REPO, '..', 'generator_predictions', 'gen2d', f'{g}_{fam}_{mode}_fte.root')
            if not os.path.exists(p): continue
            f = ROOT.TFile.Open(p); h = f.Get(c['gen_obs'] + '_fte')
            if not h or h.GetNbinsX() != ntrue:
                print(f"NBINS {c['name']} {g} {mode}: fte {h.GetNbinsX() if h else None} vs config {ntrue}"); bad += 1
            f.Close()
print('check_ext_bins:', 'OK' if bad == 0 else f'{bad} problem(s)')
sys.exit(1 if bad else 0)
