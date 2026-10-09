"""mcs_cal_derive.py -- calibration of the simulated MCS momentum of exiting muon candidates to beam-on data, from the
opened control regions of the single-pion analysis. Cause (scripts/cr_mcs_cause.py): an angular noise of reconstructed
tracks on the scale of the MCS segments that the simulation lacks, rising with the drift-direction component and
falling with the visible length.

Central calibration ("nominal"): per cell of |dir x| and visible length of the candidate track, a quantile map g with
which the MCS momentum of the simulated candidates reproduces the data. Only simulated candidates are mapped (overlay
and dirt candidates matched to a simulated particle, backtracked_pdg != 0); beam-off events and overlaid cosmic tracks
are data already, so the target distribution of the simulated part is (F_data - u F_unmapped) / (1 - u), u the
unmapped share of the prediction in the cell. Normalisation free; above the 99.5% quantile g(p)/p is held constant.

Variants for the calibration uncertainty: "shrink" (the parametric model p' = p / sqrt(a^2 + (s p / k)^2), k = 13.6 MeV,
fitted per cell), "coarse" and "fine" cell grids, and the quantile maps derived from period groups alone (old / new
production, FHC / RHC). Cells with fewer than NMIN data candidates use the union of the length bins of their |dir x| row.

Outputs:
  /data/uboone/processed/mcscal/maps/<variant>.txt     per cell: adx_lo adx_hi len_lo len_hi, then g on 0-5 GeV/c
                                                       in 0.01 steps (read by macros/mcs_cal_apply.C)
  report/planning/multipion/mcs_cal_derive.{md,json}   parameters and closure in the control regions

    python3 scripts/mcs_cal_derive.py
"""
import json, os, pickle, sys
import numpy as np
import uproot
import awkward as ak
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun
import cr_mp2pi_check as C

S = 'CC1mu1piXp'
K = 0.0136                                                     # GeV rad
OUT = '/data/uboone/processed/mcscal/'
CACHE = OUT + 'cr_candidates.pkl'
G = np.round(np.arange(0., 5.0001, 0.01), 3)                    # tabulation grid of g
X = np.arange(0., 20.0001, 0.002)                               # fine grid of the CDFs
MB = np.concatenate([np.arange(0., 2.0001, 0.1), [2.25, 2.5, 3.0, 4.0, 50.]])
NMIN = 200
NOMINAL = ([0., 0.3, 0.6, 0.85, 1.01], [0., 100., 200., 1e4])
VARIANTS = {'nominal': ('quantile', NOMINAL, None),
            'shrink': ('shrink', NOMINAL, None),
            'coarse': ('quantile', ([0., 0.45, 1.01], [0., 150., 300., 1e4]), None),
            'fine': ('quantile', ([0., 0.3, 0.6, 0.85, 1.01], [0., 50., 100., 150., 200., 300., 1e4]), None),
            'old production': ('quantile', NOMINAL, ['FHC_R1', 'RHC_R1', 'RHC_R3']),
            'new production': ('quantile', NOMINAL, ['FHC_R4', 'FHC_R5', 'RHC_R4']),
            'FHC': ('quantile', NOMINAL, ['FHC_R1', 'FHC_R4', 'FHC_R5']),
            'RHC': ('quantile', NOMINAL, ['RHC_R1', 'RHC_R3', 'RHC_R4'])}
EB = np.array([0.15, 0.3, 0.425, 0.575, 0.825, 1.25, 50.])    # released p_mu bins
SPLIT = {'abs(dir x)': ('adx', [0, 0.3, 0.6, 0.85, 1.01]), 'dir z': ('dz', [-1.01, 0.5, 0.8, 0.95, 1.01]),
         'length [cm]': ('len', [0, 50, 100, 150, 200, 300, 1e4]), 'start x [cm]': ('sx', [-10, 64, 128, 192, 270])}


def first(arr, idx):
    return ak.to_numpy(ak.fill_none(ak.firsts(arr[ak.local_index(arr) == idx]), np.nan)).astype(float)


def load(path, scale, weighted, kind, period, R):
    """kind 0 data, 1 simulation (overlay, dirt), 2 beam-off"""
    br = [f'{S}_CandidateMuonIndex', f'{S}_CandidateMuonTrackEndContainment', f'{S}_candidate_muon_mom_mcs',
          f'{S}_candidate_muon_mom_reco', f'{S}_candidate_muon_mom_range', 'trk_len_v', 'trk_dir_x_v', 'trk_dir_z_v',
          'trk_sce_start_x_v'] + (['backtracked_pdg', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in uproot.open(path)['stv_tree'].iterate(br, step_size='200 MB', library='ak'):
        if len(a) == 0: continue
        mi = ak.to_numpy(a[f'{S}_CandidateMuonIndex']); ok = mi >= 0
        idx = ak.Array(np.where(ok, mi, 0))
        w = C.safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) if weighted \
            else np.ones(len(a))
        reco = ak.to_numpy(a[f'{S}_candidate_muon_mom_reco']).astype(float)
        rng = ak.to_numpy(a[f'{S}_candidate_muon_mom_range']).astype(float)
        mcsm = ok & (reco > 0) & (reco != rng)                 # MCS-measured, as macros/mcs_cal_apply.C
        sim = (first(a['backtracked_pdg'], idx) != 0) if kind == 1 else np.zeros(len(a), bool)
        R.append(dict(mcs=reco, reco=reco, adx=np.abs(first(a['trk_dir_x_v'], idx)), dz=first(a['trk_dir_z_v'], idx),
                      len=first(a['trk_len_v'], idx), sx=first(a['trk_sce_start_x_v'], idx), w=w * scale, ok=ok,
                      exit=mcsm, sim=sim & mcsm, period=np.full(len(a), period)))


def candidates():
    if os.path.exists(CACHE): return pickle.load(open(CACHE, 'rb'))
    mc, ext, data = F.samples()
    re = lambda q: C.NEW + os.path.basename(q)
    dirts = dirt_perrun('comb', where=lambda run: C.NEW)
    RD, RP = [], []
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        for f in data[ip]: load(re(f), 1.0, False, 0, ip, RD)
        for f, s in mc[ip]: load(re(f), s, True, 1, ip, RP)
        for f, s in ext[ip]: load(re(f), s, False, 2, ip, RP)
        for r, f, s in dirts:
            if r == F.RUN[ip]: load(f, s, True, 1, ip, RP)
        print('==', name, flush=True)
    cat = lambda R: {k: np.concatenate([r[k] for r in R]) for k in R[0]}
    DP = (cat(RD), cat(RP))
    pickle.dump(DP, open(CACHE, 'wb'))
    return DP


def wcdf(p, w, x):
    o = np.argsort(p); p, w = p[o], w[o]
    c = np.cumsum(w); c /= c[-1]
    return np.interp(x, p, c, left=0., right=1.)


def qmap(d, ps, ws, pu, wu):
    """quantile map of the simulated part (ps, ws) given data d and the unmapped prediction (pu, wu)"""
    u = wu.sum() / (wu.sum() + ws.sum()) if len(pu) else 0.
    Fd = wcdf(d, np.ones(len(d)), X)
    T = (Fd - u * wcdf(pu, wu, X)) / (1. - u) if len(pu) else Fd
    T = np.maximum.accumulate(np.clip(T, 0., 1.))
    T = T + np.arange(len(X)) * 1e-12                                  # strictly increasing for the inverse
    Fs = wcdf(ps, ws, G)
    g = np.interp(Fs, T, X)
    pmax = float(np.interp(0.995, wcdf(ps, ws, X), X))
    k = float(np.interp(pmax, G, g) / pmax) if pmax > 0 else 1.
    g = np.where(G > pmax, k * G, g)
    g = np.maximum.accumulate(g); g[0] = 0.
    return g


def shrink(p, a, s):
    return p / np.sqrt(a * a + (s * p / K) ** 2)


def chi2_shrink(par, d, ps, ws, pu, wu):
    a, s = par
    if not (0.7 < a < 1.3 and 0. <= s < 0.05): return 1e12
    x = np.clip(shrink(ps, a, s), 0, 49.9)
    h = np.histogram(x, MB, weights=ws)[0] + np.histogram(np.clip(pu, 0, 49.9), MB, weights=wu)[0]
    h2 = np.histogram(x, MB, weights=ws ** 2)[0] + np.histogram(np.clip(pu, 0, 49.9), MB, weights=wu ** 2)[0]
    hd = np.histogram(np.clip(d, 0, 49.9), MB)[0]
    n = hd.sum() / h.sum(); var = hd + n * n * h2; m = var > 0
    return float(((hd - n * h)[m] ** 2 / var[m]).sum())


def smap(d, ps, ws, pu, wu):
    args = (d, ps, ws, pu, wu)
    r = min((minimize(chi2_shrink, [1.0, s0], args=args, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-4))
             for s0 in (0.003, 0.008)), key=lambda r: r.fun)
    return shrink(G, *r.x), dict(a=float(r.x[0]), s=float(r.x[1]), chi2=float(r.fun))


def cell_mask(X_, edges, i, j):
    ex, el = edges
    return X_['exit'] & (X_['adx'] >= ex[i]) & (X_['adx'] < ex[i + 1]) & (X_['len'] >= el[j]) & (X_['len'] < el[j + 1])


def derive(D, P, method, edges, sel_d, sel_p):
    ex, el = edges
    maps, info = {}, {}
    for i in range(len(ex) - 1):
        for j in range(len(el) - 1):
            md, mp = cell_mask(D, edges, i, j) & sel_d, cell_mask(P, edges, i, j) & sel_p
            n = int(md.sum()); merged = n < NMIN
            if merged:                                                    # union of the length bins of the row
                md = sum(cell_mask(D, edges, i, jj) for jj in range(len(el) - 1)).astype(bool) & sel_d
                mp = sum(cell_mask(P, edges, i, jj) for jj in range(len(el) - 1)).astype(bool) & sel_p
            ms, mu = mp & P['sim'], mp & ~P['sim']
            args = (D['mcs'][md], P['mcs'][ms], P['w'][ms], P['mcs'][mu], P['w'][mu])
            if method == 'quantile':
                g, extra = qmap(*args), {}
            else:
                g, extra = smap(*args)
            maps[(i, j)] = g
            info[(i, j)] = dict(n=n, merged=merged, u=float(P['w'][mu].sum() / max(P['w'][mp].sum(), 1e-9)),
                                g_over_p={str(q): float(np.interp(q, G, g) / q) for q in (0.5, 1.0, 1.5, 2.0, 3.0)}, **extra)
    return maps, info


def apply(P, edges, maps):
    ex, el = edges
    out = P['reco'].copy()
    for (i, j), g in maps.items():
        m = cell_mask(P, edges, i, j) & P['sim']
        p = P['reco'][m]
        out[m] = np.where(p <= G[-1], np.interp(p, G, g), p * g[-1] / G[-1])
    return out


def write_map(name, edges, maps):
    os.makedirs(OUT + 'maps', exist_ok=True)
    ex, el = edges
    with open(OUT + 'maps/' + name.replace(' ', '_') + '.txt', 'w') as f:
        f.write(f'# MCS calibration "{name}" (scripts/mcs_cal_derive.py): per cell adx_lo adx_hi len_lo len_hi, then g(p) for '
                f'p = 0, 0.01, ..., {G[-1]:g} GeV/c; above {G[-1]:g} GeV/c g(p)/p is constant\n')
        for (i, j), g in maps.items():
            f.write(f'{ex[i]:g} {ex[i + 1]:g} {el[j]:g} {el[j + 1]:g} ' + ' '.join(f'{x:.5f}' for x in g) + '\n')


def closure(D, P, pm):
    """released p_mu bins (all candidates) and the high-MCS fraction of exiting candidates, data / prediction"""
    r = {}
    hd = np.histogram(np.clip(D['reco'][D['ok']], 0, 49.9), EB)[0]
    hp = np.histogram(np.clip(pm[P['ok']], 0, 49.9), EB, weights=P['w'][P['ok']])[0]
    r['pmu'] = ((hd / hp) / (hd.sum() / hp.sum())).tolist()
    ed, ep = D['exit'], P['exit']
    for nm, (key, bins) in SPLIT.items():
        row = []
        for b in range(len(bins) - 1):
            a_ = ed & (D[key] >= bins[b]) & (D[key] < bins[b + 1]); b_ = ep & (P[key] >= bins[b]) & (P[key] < bins[b + 1])
            fd = (a_ & (D['reco'] > 1.2)).sum() / max(a_.sum(), 1)
            fp = P['w'][b_ & (pm > 1.2)].sum() / max(P['w'][b_].sum(), 1e-9)
            row.append(fd / fp if fp > 0 else float('nan'))
        r[nm] = row
    rows = []
    for ip, name in enumerate(F.PNAME):
        a_, b_ = ed & (D['period'] == ip), ep & (P['period'] == ip)
        if a_.sum() == 0: continue
        fd = (a_ & (D['reco'] > 1.2)).sum() / a_.sum(); fp = P['w'][b_ & (pm > 1.2)].sum() / P['w'][b_].sum()
        rows.append(f'{name} {fd / fp:.2f}')
    r['periods'] = rows
    hd = np.histogram(np.clip(D['reco'][ed], 0, 49.9), MB)[0]
    hp = np.histogram(np.clip(pm[ep], 0, 49.9), MB, weights=P['w'][ep])[0]
    hp2 = np.histogram(np.clip(pm[ep], 0, 49.9), MB, weights=P['w'][ep] ** 2)[0]
    n = hd.sum() / hp.sum(); var = hd + n * n * hp2; m = var > 0
    r['chi2_mcs'] = [float(((hd - n * hp)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)]
    return r


def main():
    D, P = candidates()
    out = {}
    L = ['# MCS calibration of exiting muon candidates: simulation to beam-on data', '',
         'Generated by `scripts/mcs_cal_derive.py` from the control-region skim. Central calibration: per cell of |dir x| '
         'and visible length, a quantile map of the MCS momentum of simulated candidates (overlay and dirt candidates matched '
         'to a simulated particle) to the data; beam-off and overlaid cosmic tracks are not mapped. Variants for the '
         'uncertainty: the parametric model p\' = p / sqrt(a^2 + (s p / k)^2) ("shrink"), coarser and finer grids, and maps '
         'from period groups alone.', '',
         f'Simulated share of the exiting prediction: {P["w"][P["sim"]].sum() / P["w"][P["exit"]].sum():.3f}.', '']
    pre = closure(D, P, P['reco'])
    L += ['## Closure in the control regions (all periods)', '',
          'Released p_mu bins: data / prediction, shape. High-MCS fraction: exiting candidates with MCS > 1.2 GeV/c, data / '
          'prediction. chi2: MCS spectrum of exiting candidates, shape.', '',
          '| calibration | p_mu bins | abs(dir x) | dir z | length | start x | chi2 MCS |', '|---|---|---|---|---|---|---|']
    fmt = lambda v: ' '.join(f'{x:.2f}' for x in v)
    L.append(f'| none | {fmt(pre["pmu"])} | {fmt(pre["abs(dir x)"])} | {fmt(pre["dir z"])} | {fmt(pre["length [cm]"])} | '
             f'{fmt(pre["start x [cm]"])} | {pre["chi2_mcs"][0]:.0f} / {pre["chi2_mcs"][1]} |')
    out['none'] = dict(closure=pre)
    period_rows = [f'| none | ' + ', '.join(pre['periods']) + ' |']
    for name, (method, edges, groups) in VARIANTS.items():
        ids = [F.PNAME.index(n) for n in groups] if groups else list(range(len(F.PNAME)))
        maps, info = derive(D, P, method, edges, np.isin(D['period'], ids), np.isin(P['period'], ids))
        write_map(name, edges, maps)
        pm = apply(P, edges, maps)
        c = closure(D, P, pm)
        out[name] = dict(method=method, edges=[list(map(float, e)) for e in edges], periods=groups,
                         cells={f'{i},{j}': v for (i, j), v in info.items()}, closure=c)
        L.append(f'| {name} | {fmt(c["pmu"])} | {fmt(c["abs(dir x)"])} | {fmt(c["dir z"])} | {fmt(c["length [cm]"])} | '
                 f'{fmt(c["start x [cm]"])} | {c["chi2_mcs"][0]:.0f} / {c["chi2_mcs"][1]} |')
        period_rows.append(f'| {name} | ' + ', '.join(c['periods']) + ' |')
        print('==', name, flush=True)
    L += ['', 'Bins: p_mu ' + ', '.join(f'{EB[b]:g}-{EB[b + 1]:g}' for b in range(len(EB) - 1)) + '; ' +
          '; '.join(f'{nm} ' + ', '.join(f'{bins[b]:g}-{bins[b + 1]:g}' for b in range(len(bins) - 1)) for nm, (_, bins) in SPLIT.items()), '',
          '## High-MCS fraction per period, data / prediction', '', '| calibration | periods |', '|---|---|'] + period_rows
    L += ['', '## Central calibration per cell: g(p) / p (data candidates, unmapped share u)', '',
          '| abs(dir x) | length [cm] | data | u | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 GeV/c | shrink a | shrink s [mrad] |',
          '|---|---|---|---|---|---|---|---|---|---|---|']
    ex, el = NOMINAL
    for key, v in out['nominal']['cells'].items():
        i, j = map(int, key.split(',')); sv = out['shrink']['cells'][key]
        L.append(f'| {ex[i]:g}-{ex[i + 1]:g} | {el[j]:g}-{el[j + 1]:g} | {v["n"]}{" (row)" if v["merged"] else ""} | {v["u"]:.3f} | '
                 + ' | '.join(f'{x:.3f}' for x in v['g_over_p'].values()) + f' | {sv["a"]:.3f} | {1e3 * sv["s"]:.2f} |')
    json.dump(out, open(os.path.join(C.OUTDIR, 'mcs_cal_derive.json'), 'w'), indent=1, default=float)
    open(os.path.join(C.OUTDIR, 'mcs_cal_derive.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
