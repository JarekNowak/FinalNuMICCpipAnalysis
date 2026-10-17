"""cr_mp2pi_syst_check.py -- covariance-aware comparison of the two-pion classifier inputs and output with beam-on data
in the opened control regions of the single-pion analysis: the shape chi2 of scripts/cr_mp2pi_check.py, which carries
the data and simulation statistics only, against the same test with the systematic covariance of the prediction.

Skim: /data/uboone/processed/mp2pi_cr_skim_shr (MCS-calibrated classifier inputs; the two-pion signal region is absent).
Quantities: the event score and the summed P(pi) of the two-pion domain (regions 1-4), and over the tracks of the
classifier's pool (events with one slice and the vertex in the fiducial volume, as scripts/cr_mp2pi_check.py) P(mu),
P(pi), the MCS momentum and the LLR score. Covariance of the prediction (overlay, beam-off,
dirt): flux (weight_ppfx_all, 600 universes, replacing the PPFX CV weight), cross section (weight_All_UBGenie, 600,
replacing the tune), reinteraction (weight_reint_all, 1000), the eight Run-4 detector variations taken fractionally per
horn mode (variation over variation CV in the same control-region events, applied to the nominal overlay of that horn
mode; the statistics of the variation samples are included in this term), simulation statistics; data statistics.
The normalisation is removed (shape of each universe and variation normalised to the CV), and the last bin is dropped
for the inversion. Detector ratios: event-level quantities from the two-pion production (mp2pi_mcscal/detvar_fluxw,
flux-weighted), track-level ones from the LLR-and-geometry dump (mp_evt_llronly_mcscal/detvar).

    python3 scripts/cr_mp2pi_syst_check.py  ->  report/planning/multipion/cr_mp2pi_syst_check.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun

SK = '/data/uboone/processed/mp2pi_cr_skim_shr/'
DVE = '/data/uboone/processed/mp2pi_mcscal/detvar_fluxw/'
DVT = '/data/uboone/processed/mp_evt_llronly_mcscal/detvar/'
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
S, B, BE = 'CC1mu1piXp', 'CC1mu2piBDTEvt_', 'CC1mu2piEvt_'
CR = ['sb_cc0pi', 'sb_multipi', 'sb_pi0', 'sb_cosmic']
KNOBS = ['LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']
UNIV = {'flux': ('weight_ppfx_all', 'ppfx'), 'cross section': ('weight_All_UBGenie', 'tune'), 'reinteraction': ('weight_reint_all', 'cv')}
EVQ = {'event score': np.linspace(0., 0.92, 13), 'summed P(pi)': np.linspace(0., 6., 13)}
TRQ = {'P(mu)': ('etrk_pid_mu', np.linspace(0, 1, 11)), 'P(pi)': ('etrk_pid_pi', np.linspace(0, 1, 11)),
       'MCS momentum': ('etrk_mom_mcs', np.array([0, .2, .4, .6, .8, 1., 1.2, 1.5, 2., 50.])),
       'LLR score': ('etrk_llr', np.linspace(-1, 1, 11))}


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def hist(x, bins, w):
    return np.histogram(np.clip(x, bins[0], bins[-1] - 1e-9), bins, weights=w)[0]


def crmask(a, S_=S):
    m = np.zeros(len(a), bool)
    for c in CR: m |= ak.to_numpy(a[f'{S_}_{c}']).astype(bool)
    return m


def fill(path, scale, kind, H):
    """kind: data, mc (universes), ext, dirt (CV only)"""
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_{c}' for c in CR] + [B + x for x in ('evt_region', 'evt_score', 'evt_features', 'evt_presel')] + \
         [B + v for v, _ in TRQ.values()] + [B + 'etrk_idx']
    if kind in ('mc', 'dirt'): br += ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
    if kind == 'mc': br += [u for u, _ in UNIV.values()]
    br = [b for b in br if b in t.keys()]
    for a in t.iterate(br, step_size='150 MB', library='ak'):
        a = a[crmask(a)]
        if len(a) == 0: continue
        dom = ak.to_numpy(a[B + 'evt_region']) >= 1
        pre = ak.to_numpy(a[B + 'evt_presel']).astype(bool)
        if kind in ('mc', 'dirt'):
            tc, pc, nw = (ak.to_numpy(a[x]).astype(float) for x in ('tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'))
            w = safe(tc * pc * nw) * scale
        else:
            w = np.full(len(a), scale)
        ev = {'event score': ak.to_numpy(a[B + 'evt_score']), 'summed P(pi)': ak.to_numpy(ak.fill_none(ak.pad_none(a[B + 'evt_features'], 15)[:, 14], -1.))}
        ntr = ak.to_numpy(ak.num(a[B + 'etrk_idx'])) * pre
        tr = {q: ak.to_numpy(ak.flatten(a[B + v][pre])) for q, (v, _) in TRQ.items()}
        wt = np.repeat(w, ntr)
        for q, bins in EVQ.items():
            H[(kind, q, 'cv')] = H.get((kind, q, 'cv'), 0.) + hist(ev[q][dom], bins, w[dom])
            H[(kind, q, 'w2')] = H.get((kind, q, 'w2'), 0.) + hist(ev[q][dom], bins, w[dom] ** 2)
        for q, (_, bins) in TRQ.items():
            H[(kind, q, 'cv')] = H.get((kind, q, 'cv'), 0.) + hist(tr[q], bins, wt)
            H[(kind, q, 'w2')] = H.get((kind, q, 'w2'), 0.) + hist(tr[q], bins, wt * wt)
        if kind != 'mc': continue
        for term, (ub, conv) in UNIV.items():
            if ub not in a.fields: continue
            U = a[ub]; n = int(ak.max(ak.num(U)))
            full = ak.to_numpy(ak.num(U)) == n
            Um = np.ones((len(a), n)); Um[full] = ak.to_numpy(U[full]).astype(float)
            base = (pc * nw) if conv == 'tune' else (tc * nw) if conv == 'ppfx' else (tc * pc * nw)
            Wu = safe(Um * base[:, None]) * scale; Wu[~full] = w[~full, None]
            for q, bins in EVQ.items():
                idx = np.digitize(np.clip(ev[q][dom], bins[0], bins[-1] - 1e-9), bins) - 1
                h = np.zeros((len(bins) - 1, n)); np.add.at(h, idx, Wu[dom])
                H[('mc', q, term)] = H.get(('mc', q, term), 0.) + h
            Wt = np.repeat(Wu, ntr, axis=0)
            for q, (_, bins) in TRQ.items():
                idx = np.digitize(np.clip(tr[q], bins[0], bins[-1] - 1e-9), bins) - 1
                h = np.zeros((len(bins) - 1, n)); np.add.at(h, idx, Wt)
                H[('mc', q, term)] = H.get(('mc', q, term), 0.) + h


def detvar_ratios():
    """per horn mode and variation: histogram of each quantity in control-region events, variation over CV"""
    R = {}
    for horn in ('fhc', 'rhc'):
        hv = {}
        for v in ['CV'] + KNOBS:
            h = {}
            t = uproot.open(f'{DVE}xsec-ana-detvar_run4{horn}_{v}.root')['stv_tree']
            br = [f'{S}_{c}' for c in CR] + ['CC1mu2piBDT_evt_region', 'CC1mu2piBDT_evt_score', 'CC1mu2piBDT_evt_features',
                                             'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight', 'numi_flux_rw_weight']
            for a in t.iterate([b for b in br if b in t.keys()], step_size='200 MB', library='ak'):
                m = crmask(a) & (ak.to_numpy(a['CC1mu2piBDT_evt_region']) >= 1); a = a[m]
                if len(a) == 0: continue
                w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']).astype(float))
                if 'numi_flux_rw_weight' in a.fields: w = w * ak.to_numpy(a['numi_flux_rw_weight']).astype(float)
                ev = {'event score': ak.to_numpy(a['CC1mu2piBDT_evt_score']), 'summed P(pi)': ak.to_numpy(a['CC1mu2piBDT_evt_features'][:, 14])}
                for q, bins in EVQ.items(): h[q] = h.get(q, 0.) + hist(ev[q], bins, w)
            t = uproot.open(f'{DVT}xsec-ana-detvar_run4{horn}_{v}.root')['stv_tree']
            br = [f'{S}_{c}' for c in CR] + [BE + 'evt_presel', BE + 'etrk_idx'] + [BE + x for x, _ in TRQ.values()] + \
                 ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
            for a in t.iterate([b for b in br if b in t.keys()], step_size='200 MB', library='ak'):
                m = crmask(a) & ak.to_numpy(a[BE + 'evt_presel']).astype(bool); a = a[m]
                if len(a) == 0: continue
                w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']).astype(float))
                wt = np.repeat(w, ak.to_numpy(ak.num(a[BE + 'etrk_idx'])))
                for q, (x, bins) in TRQ.items(): h[q] = h.get(q, 0.) + hist(ak.to_numpy(ak.flatten(a[BE + x])), bins, wt)
            hv[v] = h
        for v in KNOBS:
            R[(horn, v)] = {q: np.nan_to_num(hv[v][q] / np.where(hv['CV'][q] > 0, hv['CV'][q], np.nan), nan=1.) for q in hv['CV']}
    return R


def shape_dev(p_var, p):
    """shape part of a variation: normalised to the CV total"""
    return p_var * (p.sum() / p_var.sum()) - p


def test(d, p, C):
    """shape chi2: bins with neither data nor predicted variance are left out, then the last bin"""
    r = d.sum() / p.sum()
    V = r * r * C + np.diag(d)
    k = np.flatnonzero((np.diag(V) > 0) | (d > 0))[:-1]
    x = (d - r * p)[k]
    c = float(x @ np.linalg.lstsq(V[np.ix_(k, k)], x, rcond=None)[0]); ndf = len(k)
    return c, ndf, float(chi2dist.sf(c, ndf))


def main():
    import pickle
    cache = '/data/uboone/processed/mp2pi_cr_skim_shr/syst_check_hist.pkl'
    if os.path.exists(cache):
        HP, R = pickle.load(open(cache, 'rb'))
    else:
        HP, R = build()
        pickle.dump((HP, R), open(cache, 'wb'))
    report(HP, R)


def build():
    mc, ext, data = F.samples()
    re = lambda q: SK + os.path.basename(q)
    dirts = dirt_perrun('comb', where=lambda run: SK)
    HP = {'FHC': {}, 'RHC': {}}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        H = HP['FHC' if name.startswith('FHC') else 'RHC']
        for f in data[ip]: fill(re(f), 1., 'data', H)
        for f, s in mc[ip]: fill(re(f), s, 'mc', H)
        for f, s in ext[ip]: fill(re(f), s, 'ext', H)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, 'dirt', H)
        print('==', name, flush=True)
    return HP, detvar_ratios()


def report(HP, R):
    res = {}
    L = ['# Two-pion classifier inputs and output in the opened control regions: covariance-aware shape tests', '',
         'Generated by `scripts/cr_mp2pi_syst_check.py` (see its docstring for the covariance). Shape chi2/ndf (p-value) with the '
         'normalisation removed: statistics only (data and simulation), and with the full systematic covariance. FHC and RHC '
         'periods summed; the detector variations act per horn mode.', '',
         '| quantity | stat. only | full covariance | flux | cross section | reinteraction | detector | data / prediction (total) |',
         '|---|---|---|---|---|---|---|---|']
    for q in list(EVQ) + list(TRQ):
        d = sum(HP[h][('data', q, 'cv')] for h in HP)
        p_h = {h: HP[h][('mc', q, 'cv')] + HP[h].get(('ext', q, 'cv'), 0.) + HP[h].get(('dirt', q, 'cv'), 0.) for h in HP}
        p = sum(p_h.values())
        w2 = sum(HP[h][('mc', q, 'w2')] + HP[h].get(('ext', q, 'w2'), 0.) + HP[h].get(('dirt', q, 'w2'), 0.) for h in HP)
        Cs = np.diag(w2); terms = {}
        for term in UNIV:
            Pu = sum(HP[h][('mc', q, term)] + (p_h[h] - HP[h][('mc', q, 'cv')])[:, None] for h in HP)
            dev = np.stack([shape_dev(Pu[:, u], p) for u in range(Pu.shape[1])], axis=1)
            terms[term] = dev @ dev.T / dev.shape[1]
        Cd = 0.
        for v in KNOBS:
            pv = sum(p_h[h] + HP[h][('mc', q, 'cv')] * (R[(h.lower(), v)][q] - 1.) for h in HP)
            dv = shape_dev(pv, p); Cd = Cd + np.outer(dv, dv)
        terms['detector'] = Cd
        stat = test(d, p, Cs); full = test(d, p, Cs + sum(terms.values()))
        one = {t: test(d, p, Cs + C)[0] for t, C in terms.items()}
        res[q] = dict(stat=stat, full=full, single_term_chi2=one, data=list(d), pred=list(p), norm=float(d.sum() / p.sum()))
        f = lambda r: f'{r[0]:.1f}/{r[1]} ({r[2]:.3f})'
        L.append(f'| {q} | {f(stat)} | {f(full)} | ' + ' | '.join(f'{one[t]:.1f}' for t in ('flux', 'cross section', 'reinteraction', 'detector'))
                 + f' | {d.sum() / p.sum():.3f} |')
    L += ['', 'Columns flux to detector: chi2 with the statistics and that term alone.']
    json.dump(res, open(os.path.join(OUTDIR, 'cr_mp2pi_syst_check.json'), 'w'), indent=1, default=float)
    open(os.path.join(OUTDIR, 'cr_mp2pi_syst_check.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
