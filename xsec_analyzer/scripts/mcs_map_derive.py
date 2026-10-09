"""mcs_map_derive.py -- the momentum mapping that turns the simulated MCS momentum of exiting muons into the one observed
in beam-on data, from the opened control regions of the single-pion analysis (scripts/cr_mcs_check.py found the
data MCS spectrum of exiting muons shifted down above about 1.2 GeV/c while contained muons agree).

Quantile mapping on the exiting muon candidates of the control-region skim, all periods: f(p) = F_data^-1(F_pred(p)),
with F the weighted cumulative distributions of the MCS momentum (prediction: overlay, beam-off and dirt with the
scales of scripts/cr_mp2pi_check.py). Applied to the MCS-measured muons of the pseudo-data, f makes their spectrum that
of the data in the control regions; the whole difference is attributed to the MCS estimator, which the agreement of
the contained muons (range momentum) supports. Above the last populated quantile the ratio f(p)/p is held constant.

    python3 scripts/mcs_map_derive.py  ->  /data/uboone/processed/mcsmap_study/mcs_map.txt (p, f(p), GeV/c)
"""
import os, sys
import numpy as np
import uproot
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun
import cr_mp2pi_check as C

S = 'CC1mu1piXp'
OUT = '/data/uboone/processed/mcsmap_study/'
G = np.arange(0., 5.0001, 0.01)


def fill(path, scale, weighted, acc):
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_CandidateMuonIndex', f'{S}_CandidateMuonTrackEndContainment', f'{S}_candidate_muon_mom_mcs'] + \
         (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in t.iterate(br, step_size='200 MB', library='np'):
        w = C.safe(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']) if weighted else np.ones(len(a[br[0]]))
        w = w * scale
        m = (a[f'{S}_CandidateMuonIndex'] >= 0) & ~a[f'{S}_CandidateMuonTrackEndContainment'].astype(bool)
        p = a[f'{S}_candidate_muon_mom_mcs']
        m &= np.isfinite(p) & (p > 0)
        acc.append((p[m], w[m]))


def ecdf(acc):
    p = np.concatenate([x for x, _ in acc]); w = np.concatenate([y for _, y in acc])
    o = np.argsort(p); p, w = p[o], w[o]
    c = np.cumsum(w) / w.sum()
    return p, c


def main():
    os.makedirs(OUT, exist_ok=True)
    mc, ext, data = F.samples()
    re = lambda q: C.NEW + os.path.basename(q)
    dirts = dirt_perrun('comb', where=lambda run: C.NEW)
    D, P = [], []
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        for f in data[ip]: fill(re(f), 1.0, False, D)
        for f, s in mc[ip]: fill(re(f), s, True, P)
        for f, s in ext[ip]: fill(re(f), s, False, P)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, True, P)
    pd, cd = ecdf(D); pp, cp = ecdf(P)
    Fp = np.interp(G, pp, cp, left=0., right=1.)
    f = np.interp(Fp, cd, pd)
    # beyond the populated range keep f(p)/p constant (from the 99.5% quantile of the prediction)
    pmax = float(np.interp(0.995, cp, pp))
    k = float(np.interp(pmax, G, f) / pmax)
    f = np.where(G > pmax, k * G, f)
    f = np.maximum.accumulate(f)                         # monotonic
    f[0] = 0.
    with open(OUT + 'mcs_map.txt', 'w') as fp:
        fp.write(f'# MCS momentum mapping of exiting muons, simulation -> beam-on data (scripts/mcs_map_derive.py); '
                 f'{sum(len(x) for x, _ in D)} data and {sum(len(x) for x, _ in P)} simulated candidates; '
                 f'f(p)/p constant ({k:.4f}) above {pmax:.3f} GeV/c\n# p_GeV  f_GeV\n')
        for a, b in zip(G, f): fp.write(f'{a:.3f} {b:.5f}\n')
    for q in (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0):
        print(f'p {q:.1f}  f(p) {np.interp(q, G, f):.3f}  ratio {np.interp(q, G, f) / q:.3f}')
    print('written', OUT + 'mcs_map.txt')


if __name__ == '__main__':
    main()
