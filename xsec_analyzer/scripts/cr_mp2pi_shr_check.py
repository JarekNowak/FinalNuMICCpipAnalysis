"""cr_mp2pi_shr_check.py -- the shower features of the two-pion event classifier (set shr of scripts/mp_evt.py),
beam-on data against the prediction (overlay, beam-off, dirt) in the opened control regions of the single-pion
analysis, before they are used. Input: the control-region skim with the shower dump and the MCS-calibrated classifier
inputs (slurm/slurm_mp2pi_cr_skim.sbatch into /data/uboone/processed/mp2pi_cr_skim_shr; CC1mu1piXp and
CC1mu2piBDTEvt). Events of the two-pion domain (regions 1-4; the signal region is absent from the skim); the features
use the tracks the selection assigned (evt_mu_idx, evt_pi0_idx, evt_pi1_idx). Scales as scripts/cr_mp2pi_check.py.
Per feature: data / prediction of the normalisation and the shape chi2 with the normalisation removed, per period and
summed; and the share of events with a shower.

    python3 scripts/cr_mp2pi_shr_check.py [skim_dir]  ->  report/planning/multipion/cr_mp2pi_shr_check.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun
import mp_evt as M

SK = (sys.argv[1].rstrip('/') + '/') if len(sys.argv) > 1 else '/data/uboone/processed/mp2pi_cr_skim_shr/'
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
S, B = 'CC1mu1piXp', 'CC1mu2piBDTEvt_'
CR = ['sb_cc0pi', 'sb_multipi', 'sb_pi0', 'sb_cosmic']
BINS = {'shr_E1': np.array([-1, 1, 30, 60, 100, 150, 200, 300, 500, 1e9]), 'shr_E2': np.array([-1, 1, 20, 50, 100, 1e9]),
        'shr_Etot': np.array([-1, 1, 30, 60, 100, 150, 200, 300, 500, 1e9]), 'shr_n30': np.arange(-0.5, 4.5),
        'shr_d1': np.array([-2, 0, 2, 5, 10, 20, 40, 80, 1e9]), 'shr_d2': np.array([-2, 0, 5, 20, 1e9]),
        'shr_ts1': np.array([-2, 0, 0.1, 0.2, 0.3, 0.4, 0.5]), 'shr_dedx1': np.array([-2, 0, 1.5, 2.5, 3.5, 5, 8, 1e9]),
        'shr_cos1': np.array([-3, -1, 0, 0.5, 0.8, 0.95, 1.01]), 'shr_mgg': np.array([-2, 0, 50, 100, 150, 250, 1e9]),
        'shr_pi0mass': np.array([-2, 0, 50, 100, 150, 250, 1e9])}


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def fill(path, scale, weighted, H):
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_{c}' for c in CR] + [B + x for x in ('evt_region', 'evt_mu_idx', 'evt_pi0_idx', 'evt_pi1_idx',
                                                     'etrk_idx', 'etrk_dirx', 'etrk_diry', 'etrk_dirz', 'evt_pi0_mass')] \
        + [B + 'eshr_' + x for x in M.SHR] + (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    a = t.arrays(br, library='ak')
    reg = ak.to_numpy(a[B + 'evt_region'])
    cr = np.zeros(len(a), bool)
    for c in CR: cr |= ak.to_numpy(a[f'{S}_{c}']).astype(bool)
    a = a[(reg >= 1) & cr]
    if len(a) == 0: return
    incr = {c: ak.to_numpy(a[f'{S}_{c}']).astype(bool) for c in CR}
    w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) * scale if weighted \
        else np.full(len(a), scale)
    # directions of the assigned tracks, from the dump entries whose index is the assigned PFParticle
    T = {}
    for k in ('dirx', 'diry', 'dirz'):
        cols = []
        for role in ('evt_mu_idx', 'evt_pi0_idx', 'evt_pi1_idx'):
            m = a[B + 'etrk_idx'] == a[B + role]
            cols.append(ak.to_numpy(ak.fill_none(ak.firsts(a[B + 'etrk_' + k][m]), 0.)).astype(float))
        T[k] = np.stack(cols, axis=1)
    ev = M.padded_showers(a, B)
    rows = np.arange(len(a))
    f = M.shower_features(ev, T, rows, [np.zeros(len(a), int), np.ones(len(a), int), np.full(len(a), 2)])
    for k, bins in BINS.items():
        x = np.clip(f[k], bins[0], bins[-1] * (1 - 1e-12) if bins[-1] > 0 else bins[-1])
        for key, m in [(k, np.ones(len(a), bool))] + [((k, c), incr[c]) for c in CR]:
            h = H.setdefault(key, [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)])
            h[0] += np.histogram(x[m], bins, weights=w[m])[0]; h[1] += np.histogram(x[m], bins, weights=w[m] ** 2)[0]


def shape_chi2(d, p, p2):
    Dn, Pn = d.sum(), p.sum()
    if Dn <= 0 or Pn <= 0: return float('nan'), 0
    r = Dn / Pn
    var = d + p2 * r * r
    m = var > 0
    return float(((d - p * r)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)


def main():
    mc, ext, data = F.samples()
    re = lambda p: SK + os.path.basename(p)
    dirts = dirt_perrun('comb', where=lambda run: SK)
    out, HD, HP = {}, {}, {}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        Hd, Hp = {}, {}
        for f_ in data[ip]: fill(re(f_), 1., False, Hd)
        for f_, s in mc[ip]: fill(re(f_), s, True, Hp)
        for f_, s in ext[ip]: fill(re(f_), s, False, Hp)
        for r, f_, s in dirts:
            if r == F.RUN[ip]: fill(f_, s, True, Hp)
        out[name] = {k: dict(zip(('chi2', 'ndf'), shape_chi2(Hd[k][0], Hp[k][0], Hp[k][1])),
                             data=list(Hd[k][0]), pred=list(Hp[k][0])) for k in BINS if k in Hd and k in Hp}
        for H, T_ in ((Hd, HD), (Hp, HP)):
            for k, v in H.items():
                T_[k] = [x + y for x, y in zip(T_[k], v)] if k in T_ else [x.copy() for x in v]
        print(f'== {name} done', flush=True)
    out['all'] = {k: dict(zip(('chi2', 'ndf'), shape_chi2(HD[k][0], HP[k][0], HP[k][1])), data=list(HD[k][0]),
                          pred=list(HP[k][0])) for k in BINS}
    json.dump(out, open(os.path.join(OUTDIR, 'cr_mp2pi_shr_check.json'), 'w'), indent=1, default=float)
    periods = list(out)
    d, p = np.array(out['all']['shr_E1']['data']), np.array(out['all']['shr_E1']['pred'])
    L = ['# Shower features of the two-pion event classifier in the opened control regions', '',
         f'Generated by `scripts/cr_mp2pi_shr_check.py` from `{SK}`. Two-pion domain (regions 1-4) of control-region '
         'events; prediction = overlay (MCS-calibrated classifier inputs), beam-off and dirt. Shape chi2/ndf with the '
         'normalisation removed; the last columns give data over prediction per bin, all periods, each normalised to one.',
         '', f'Events with a shower: data {1 - d[0] / d.sum():.3f}, prediction {1 - p[0] / p.sum():.3f}.', '',
         '| feature | ' + ' | '.join(periods) + ' | bin edges | data / prediction (shape) |', '|---|' + '---|' * (len(periods) + 2)]
    for k, bins in BINS.items():
        dd, pp = np.array(out['all'][k]['data']), np.array(out['all'][k]['pred'])
        r = (dd / dd.sum()) / np.where(pp > 0, pp / pp.sum(), np.nan)
        L.append(f'| {k} | ' + ' | '.join(f"{out[q][k]['chi2']:.1f}/{out[q][k]['ndf']}" if k in out[q] else 'n/a' for q in periods)
                 + ' | ' + ' '.join(f'{b:g}' for b in bins[1:-1]) + ' | ' + ' '.join(f'{x:.2f}' for x in r) + ' |')
    # each control region against its own normalisation (the common normalisation above mixes the regions, whose
    # data excesses differ: the pi0 region carries the largest)
    L += ['', '## Per control region, all periods: shape chi2/ndf against the region\'s own normalisation', '',
          '| feature | sb_cc0pi: chi2/ndf (data / prediction, shape) | sb_multipi | sb_pi0: chi2/ndf (data / prediction, shape) | sb_cosmic |',
          '|---|---|---|---|---|']
    per = {}
    for k in BINS:
        cells = []
        for c in CR:
            dd, pp, p2 = HD[(k, c)][0], HP[(k, c)][0], HP[(k, c)][1]
            chi, ndf = shape_chi2(dd, pp, p2); per[f'{k} | {c}'] = dict(chi2=chi, ndf=ndf, data=list(dd), pred=list(pp))
            r = (dd / dd.sum()) / np.where(pp > 0, pp / pp.sum(), np.nan) if dd.sum() > 0 else np.full(len(dd), np.nan)
            cells.append(f'{chi:.1f}/{ndf} (' + ' '.join(f'{x:.2f}' for x in r) + ')' if c in ('sb_cc0pi', 'sb_pi0') else f'{chi:.1f}/{ndf}')
        L.append(f'| {k} | ' + ' | '.join(cells) + ' |')
    for c in CR:
        d0, p0 = HD[('shr_E1', c)][0], HP[('shr_E1', c)][0]
        L.append(f'') if c == CR[0] else None
        L.append(f'- {c}: data / prediction {d0.sum() / p0.sum():.3f}; events with a shower: data {1 - d0[0] / d0.sum():.3f}, prediction {1 - p0[0] / p0.sum():.3f}')
    out['per_region'] = per
    json.dump(out, open(os.path.join(OUTDIR, 'cr_mp2pi_shr_check.json'), 'w'), indent=1, default=float)
    open(os.path.join(OUTDIR, 'cr_mp2pi_shr_check.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
