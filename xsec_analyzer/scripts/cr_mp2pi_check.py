"""cr_mp2pi_check.py -- inputs and outputs of the two-pion selection CC1mu2piBDT, beam-on data against the prediction
in the opened control regions of the single-pion analysis (report/other_notes/multipion_note.tex, status), per period
and summed.

Input: /data/uboone/processed/mp2pi_cr_skim (slurm/slurm_mp2pi_cr_skim.sbatch): every sample processed as the
control-region skim with CC1mu1piXp and CC1mu2piBDTEvt. The skim holds only events in the four control regions of
the single-pion selection and none in the signal region of either selection, so the two-pion signal region (event
score > 0.92) is absent from data and prediction alike. Overlay and beam-off scales of bkgfit/farsb.py, dirt scales
of bkgfit/dirt.py, as scripts/cr_pid_diag_check.py.

Track level: every track of the classifier's pool in events with one slice and the vertex in the fiducial volume
(etrk_*): the four class probabilities of the LLR-and-geometry classifier and its inputs that the dump keeps (LLR,
length, vertex distance, track score, range and MCS momenta). Event level, in the two-pion domain (regions 1-4):
the event score (0 to 0.92), the populations of the four sidebands and the leading event features. Per quantity:
the normalisation ratio data / prediction and the shape chi2 with the normalisation removed.

    python3 scripts/cr_mp2pi_check.py [skim_dir] [tag]  ->  report/planning/multipion/cr_mp2pi_check[_tag].{json,md}
(e.g. the skim with the MCS-calibrated simulation, /data/uboone/processed/mp2pi_cr_skim_mcscal, tag mcscal)
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun

NEW = (sys.argv[1].rstrip('/') + '/') if len(sys.argv) > 1 else '/data/uboone/processed/mp2pi_cr_skim/'
TAG = ('_' + sys.argv[2]) if len(sys.argv) > 2 else ''
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
S, B = 'CC1mu1piXp', 'CC1mu2piBDTEvt'
CR = ['sb_cc0pi', 'sb_multipi', 'sb_pi0', 'sb_cosmic']
P01 = np.linspace(0, 1, 21)
TRK = {'pid_mu': P01, 'pid_pi': P01, 'pid_p': P01, 'pid_other': P01, 'llr': np.linspace(-1, 1, 21),
       'len': np.linspace(0, 150, 31), 'dist': np.linspace(0, 10, 21), 'ts': np.linspace(0.5, 1, 11),
       'mom_range_mu': np.linspace(0, 1.5, 31), 'mom_mcs': np.linspace(0, 2, 21)}
# event features (CC1mu1piXp::classify_event order): L2-L1, L2, summed P(pi), tracks with P(pi) > 0.5, pool size,
# topological score
FEAT = {'L2_minus_L1': (0, np.linspace(-10, 10, 21)), 'L2': (3, np.linspace(-20, 0, 21)),
        'sum_P_pi': (14, np.linspace(0, 6, 25)), 'n_P_pi_gt_0p5': (18, np.arange(-0.5, 8.5, 1)),
        'pool_size': (21, np.arange(2.5, 9.5, 1)), 'topological_score': (25, np.linspace(0, 1, 21))}
SCORE = np.linspace(0, 0.92, 24)
REGIONS = {1: 'three-pion-like', 2: 'proton-like', 3: 'score 0.5-0.92', 4: 'score < 0.5'}


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def add(H, key, x, w, bins):
    ok = np.isfinite(x)
    xc = np.clip(x[ok], bins[0], bins[-1] - 1e-9)
    cur = H.get(key) or [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1), 0., 0.]
    H[key] = [cur[0] + np.histogram(xc, bins, weights=w[ok])[0], cur[1] + np.histogram(xc, bins, weights=w[ok] ** 2)[0],
              cur[2] + float((x[ok] * w[ok]).sum()), cur[3] + float(w[ok].sum())]


def fill(path, scale, weighted, H):
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_{c}' for c in CR] + [f'{B}_etrk_{v}' for v in TRK] + \
         [f'{B}_evt_region', f'{B}_evt_score', f'{B}_evt_features'] + \
         (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in t.iterate(br, step_size='200 MB', library='ak'):
        if len(a) == 0: continue
        w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) if weighted \
            else np.ones(len(a))
        w = w * scale
        reg = ak.to_numpy(a[f'{B}_evt_region'])
        for ic, c in enumerate(CR + ['all']):
            m = np.ones(len(a), bool) if c == 'all' else ak.to_numpy(a[f'{S}_{c}'])
            # event populations: in the two-pion domain or not, and per sideband
            for r in (-1, 1, 2, 3, 4):
                add(H, ('events', c, f'region {r}'), np.zeros((m & (reg == r)).sum()), w[m & (reg == r)], np.array([-1., 1.]))
            if c != 'all': continue
            dom = m & (reg >= 1)
            add(H, ('event', 'score'), ak.to_numpy(a[f'{B}_evt_score'])[dom], w[dom], SCORE)
            fe = a[f'{B}_evt_features'][dom]
            for k, (i, bins) in FEAT.items():
                add(H, ('event', k), ak.to_numpy(fe[:, i]), w[dom], bins)
            tr = {v: a[f'{B}_etrk_{v}'] for v in TRK}
            ww = ak.to_numpy(ak.flatten(ak.broadcast_arrays(w, tr['ts'])[0]))
            for v, bins in TRK.items():
                add(H, ('track', v), ak.to_numpy(ak.flatten(tr[v])), ww, bins)


def shape_chi2(d, d2, p, p2):
    Dn, Pn = d.sum(), p.sum()
    if Dn <= 0 or Pn <= 0: return float('nan'), 0
    r = Dn / Pn
    var = d2 + p2 * r * r
    m = var > 0
    return float(((d - p * r)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)


def summarise(Hd, Hp):
    res = {}
    for key in sorted(set(Hd) | set(Hp)):
        d = Hd.get(key, [np.zeros(1), np.zeros(1), 0., 0.]); p = Hp.get(key, [np.zeros(1), np.zeros(1), 0., 0.])
        dn, pn = d[0].sum(), p[0].sum()
        ratio = dn / pn if pn > 0 else float('nan')
        err = ratio * np.sqrt(d[1].sum() / dn ** 2 + p[1].sum() / pn ** 2) if dn > 0 and pn > 0 else float('nan')
        chi2, ndf = shape_chi2(d[0], d[1], p[0], p[1]) if len(d[0]) > 1 else (float('nan'), 0)
        res[' | '.join(key)] = dict(data=dn, pred=pn, ratio=ratio, ratio_err=err, chi2=chi2, ndf=ndf,
                                    data_mean=d[2] / d[3] if d[3] else float('nan'),
                                    pred_mean=p[2] / p[3] if p[3] else float('nan'),
                                    data_hist=list(d[0]), pred_hist=list(p[0]))
    return res


def main():
    mc, ext, data = F.samples()
    re = lambda path: NEW + os.path.basename(path)
    dirts = dirt_perrun('comb', where=lambda run: NEW)
    out, HD, HP = {}, {}, {}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        Hd, Hp = {}, {}
        for f in data[ip]: fill(re(f), 1.0, False, Hd)
        for f, s in mc[ip]: fill(re(f), s, True, Hp)
        for f, s in ext[ip]: fill(re(f), s, False, Hp)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, True, Hp)
        out[name] = summarise(Hd, Hp)
        for H, T in ((Hd, HD), (Hp, HP)):
            for k, v in H.items():
                T[k] = [a + b for a, b in zip(T[k], v)] if k in T else [np.array(x) if isinstance(x, np.ndarray) else x for x in v]
        print(f'== {name} done', flush=True)
    out['all'] = summarise(HD, HP)
    json.dump(out, open(os.path.join(OUTDIR, f'cr_mp2pi_check{TAG}.json'), 'w'), indent=1, default=float)
    periods = list(out)
    head = ' | '.join(periods)
    L = ['# Two-pion selection in the opened control regions: data against the prediction', '',
         f'Generated by `scripts/cr_mp2pi_check.py` from the control-region skim of `slurm/slurm_mp2pi_cr_skim.sbatch` ({NEW})'
         ' (CC1mu1piXp and CC1mu2piBDTEvt). Only events in the four opened control regions of the single-pion analysis,'
         ' none in the signal region of either selection: the two-pion signal region (score > 0.92) is absent. '
         'Prediction: overlay, beam-off and dirt. Cells: data / prediction (ratio).', '',
         '## Event populations by control region and two-pion region (region -1: outside the two-pion domain)', '',
         f'| control region | two-pion region | {head} |', '|---|---|' + '---|' * len(periods)]
    for c in CR + ['all']:
        for r in (-1, 1, 2, 3, 4):
            k = f'events | {c} | region {r}'
            cells = []
            for p in periods:
                x = out[p].get(k)
                cells.append(f'{x["data"]:.0f} / {x["pred"]:.1f} ({x["ratio"]:.2f} ± {x["ratio_err"]:.2f})'
                             if x and x['pred'] > 0 else 'n/a')
            L.append(f'| {c} | {REGIONS.get(r, "outside")} | ' + ' | '.join(cells) + ' |')
    L += ['', '## Shapes: chi2/ndf with the normalisation removed (mean data - prediction)', '',
          f'| level | quantity | {head} |', '|---|---|' + '---|' * len(periods)]
    for lev, keys in (('event', ['score'] + list(FEAT)), ('track', list(TRK))):
        for v in keys:
            k = f'{lev} | {v}'
            cells = []
            for p in periods:
                x = out[p].get(k)
                cells.append(f'{x["chi2"]:.1f}/{x["ndf"]} ({x["data_mean"] - x["pred_mean"]:+.3g})' if x and x['ndf'] > 0 else 'n/a')
            L.append(f'| {lev} | {v} | ' + ' | '.join(cells) + ' |')
    open(os.path.join(OUTDIR, f'cr_mp2pi_check{TAG}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
