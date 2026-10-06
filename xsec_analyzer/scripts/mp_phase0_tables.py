#!/usr/bin/env python3
"""Markdown tables of the Phase 0 multi-pion baseline from report/planning/multipion/phase0_baseline.json
(written by scripts/mp_phase0_baseline.py). Prints to stdout.
Usage: python3 scripts/mp_phase0_tables.py > ../report/planning/multipion/phase0_tables.md
"""
import json, os

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
J = json.load(open(REPO + '/report/planning/multipion/phase0_baseline.json'))
SELS = {'CC1mu2pi': 2, 'CC1mu3pi': 3}
FINAL = 6  # index of the final selection stage
STAGES = ['None', 'Vertex in FV', 'Topological', 'MuonCandidate', 'N pions', 'Opening angle', 'Final (selected)',
          'diag: + shower veto', 'diag: + 3-plane mu, pi']
TOPO = ['signal', 'CC0pi', 'CC1pi', 'CC2pi', 'CC3pi+', 'CCpi0', 'CCother', 'NC', 'nueCC', 'OOFV']
CFGS = ['FHC', 'RHC', 'COMB']


def v(cfg, *k):
    return J[cfg].get('|'.join(map(str, k)), 0.0)


def f(x, d=1):
    return ('%.' + str(d) + 'f') % x


print('## 1π cross-check against the released cut-flow (final stage = Selected)\n')
print('| Config | Signal | ν background | Beam-off | Dirt | Purity |')
print('|---|---|---|---|---|---|')
for c in CFGS:
    s, b, e, d = v(c, 'v1pi', 'sig'), v(c, 'v1pi', 'bkg'), v(c, 'v1pi', 'ext'), v(c, 'v1pi', 'dirt')
    print('| %s | %s | %s | %s | %s | %s%% |' % (c, f(s), f(b), f(e), f(d, 2), f(100 * s / (s + b + e + d))))

for s, n in SELS.items():
    print('\n## %s: cut-flow at the full exposure\n' % s)
    for c in CFGS:
        gen = v(c, 'cf', s, 0, 'signal')
        print('\n### %s\n' % c)
        print('| Stage | Signal | ν background | Beam-off | Dirt | Efficiency | Purity |')
        print('|---|---|---|---|---|---|---|')
        for i, name in enumerate(STAGES):
            sig = v(c, 'cf', s, i, 'signal')
            bkg = sum(v(c, 'cf', s, i, t) for t in TOPO[1:])
            ex, di = v(c, 'cf', s, i, 'ext'), v(c, 'cf', s, i, 'dirt')
            pred = sig + bkg + ex + di
            print('| %s | %s | %s | %s | %s | %s%% | %s%% |' % (name, f(sig), f(bkg), f(ex), f(di, 2),
                  f(100 * sig / gen, 2), f(100 * sig / pred, 1) if pred else '0'))
    print('\n### %s: composition of the selected sample\n' % s)
    print('| Class | ' + ' | '.join(CFGS) + ' |')
    print('|---|' + '---|' * len(CFGS))
    rows = TOPO + ['ext', 'dirt']
    for t in rows:
        vals = []
        for c in CFGS:
            tot = sum(v(c, 'cf', s, FINAL, x) for x in rows)
            x = v(c, 'cf', s, FINAL, t)
            vals.append('%s (%s%%)' % (f(x), f(100 * x / tot)))
        print('| %s | %s |' % (t, ' | '.join(vals)))
    print('\n### %s: subsamples, phase space and conventions (selected signal unless stated)\n' % s)
    print('| Quantity | ' + ' | '.join(CFGS) + ' |')
    print('|---|' + '---|' * len(CFGS))
    def row(label, fn):
        print('| %s | %s |' % (label, ' | '.join(fn(c) for c in CFGS)))
    row('≥1p subsample: generated / selected / efficiency',
        lambda c: '%s / %s / %s%%' % (f(v(c, 'sig1p', s, 'gen')), f(v(c, 'sig1p', s, 'sel')),
                                     f(100 * v(c, 'sig1p', s, 'sel') / max(v(c, 'sig1p', s, 'gen'), 1e-9), 2)))
    row('≥1p share of the generated signal',
        lambda c: '%s%%' % f(100 * v(c, 'sig1p', s, 'gen') / max(v(c, 'cf', s, 0, 'signal'), 1e-9)))
    row('all pions above 0.175 GeV/c: share of generated / of selected signal',
        lambda c: '%s%% / %s%%' % (f(100 * v(c, 'thr0175', s, 'gen') / max(v(c, 'cf', s, 0, 'signal'), 1e-9)),
                                   f(100 * v(c, 'thr0175', s, 'sel') / max(v(c, 'cf', s, FINAL, 'signal'), 1e-9))))
    row('removed only by θ(μ, leading π) < 2.6: generated / selected',
        lambda c: '%s / %s' % (f(v(c, 'thetaonly', s, 'gen')), f(v(c, 'thetaonly', s, 'sel'))))
    row('longest candidate is the true leading pion',
        lambda c: '%s%%' % f(100 * v(c, 'lead', s, 'longest') / max(v(c, 'lead', s, 'n'), 1e-9)))
    row('stored candidate (highest LLR) is the true leading pion',
        lambda c: '%s%%' % f(100 * v(c, 'lead', s, 'llrcand') / max(v(c, 'lead', s, 'n'), 1e-9)))
    row('at least N candidates are true π±',
        lambda c: '%s%%' % f(100 * v(c, 'lead', s, 'allNtruepi') / max(v(c, 'lead', s, 'n'), 1e-9)))
    row('selected events with an uncontained pion counted (signal / background)',
        lambda c: '%s%% / %s%%' % (f(100 * v(c, 'unc', s, 'sig') / max(v(c, 'cf', s, FINAL, 'signal'), 1e-9)),
                                   f(100 * (v(c, 'unc', s, 'bkg') + v(c, 'unc', s, 'ext') + v(c, 'unc', s, 'dirt'))
                                     / max(sum(v(c, 'cf', s, FINAL, t) for t in TOPO[1:] + ['ext', 'dirt']), 1e-9))))
    print('\n### %s: overlaps with the inclusive selection (selected sample, all components)\n' % s)
    print('| Overlap | ' + ' | '.join(CFGS) + ' |')
    print('|---|' + '---|' * len(CFGS))
    def ovl(c, what):
        tot = sum(v(c, 'cf', s, FINAL, t) for t in TOPO + ['ext', 'dirt'])
        x = sum(v(c, 'ovl', s, what, k) for k in ('sig', 'bkg', 'ext', 'dirt'))
        return '%s (%s%%)' % (f(x), f(100 * x / tot))
    print('| also selected by the 1π signal region | %s |' % ' | '.join(ovl(c, 'sel1pi') for c in CFGS))
    print('| also in the opened 1π multi-π control region | %s |' % ' | '.join(ovl(c, 'cr_multipi') for c in CFGS))
    print('| truth: also 1π signal (must be 0) | %s |' % ' | '.join(f(v(c, 'excl', s, 'sig1pi'), 3) for c in CFGS))

print('\n## Two-pion and three-pion selections together\n')
print('| Quantity | ' + ' | '.join(CFGS) + ' |')
print('|---|' + '---|' * len(CFGS))
print('| events selected by both (MC + beam-off + dirt) | %s |' % ' | '.join(
    f(sum(v(c, 'ovl', '2pi3pi', k) for k in ('mc', 'ext', 'dirt'))) for c in CFGS))
print('| truth: 2π and 3π signal at once (must be 0) | %s |' % ' | '.join(f(v(c, 'excl', '2pi3pi'), 3) for c in CFGS))
print('| 1π multi-π control region, total prediction | %s |' % ' | '.join(
    f(sum(v(c, 'v1pi', 'cr_multipi', k) for k in ('mc', 'ext', 'dirt'))) for c in CFGS))
