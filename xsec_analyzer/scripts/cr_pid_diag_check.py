"""cr_pid_diag_check.py -- the Bragg likelihoods, BDT scores and particle-classifier probabilities of
the pion identification, beam-on data against the prediction in the opened control regions, per
period. Step 2 of the check whether the periods of the older ntuple production (FHC Run 1, RHC Runs 1
and 3b, whose beam-on, beam-off and dirt ntuples lack trk_bragg_pion_v) differ from the newer ones
(report/planning/multipion/PHASE2_SUMMARY.md); step 1 is scripts/cr_pid_production_check.py.

Input: /data/uboone/processed/pid_diag_skim (slurm/slurm_pid_diag_skim.sbatch): every sample of the
comparison processed as the control-region skim with CC1mu1piXp and its diagnostic twin
CC1mu1piXpPIDDiag (pdiag_* per pion-pool candidate, mudiag_* for the muon candidate). Overlay and
beam-off scales of bkgfit/farsb.py, dirt scales of bkgfit/dirt.py.

Populations: "pool", the pion-pool candidates with track score >= 0.5; "eligible", the pool
candidates the released pion identification can count (contained, vertex distance <= 4 cm, length
> 20 cm); "muon", the muon candidate. Per variable and period: the shape chi2 (normalisation
removed) and the mean difference. Per cut and period: the pass fraction of the eligible candidates
in data and prediction.

    python3 scripts/cr_pid_diag_check.py  ->  report/planning/multipion/cr_pid_diag_check.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun

NEW = '/data/uboone/processed/pid_diag_skim/'
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
S, D = 'CC1mu1piXp', 'CC1mu1piXpPIDDiag'
OLD = ('FHC_R1', 'RHC_R1', 'RHC_R3')
BR, BD, BP = np.linspace(0, 1.25, 26), np.linspace(-0.4, 0.45, 18), np.linspace(0, 1, 21)
POOLV = {'bragg_p': BR, 'bragg_mu': BR, 'bragg_mip': BR, 'bragg_pion': BR, 'mip_bdt': BD, 'pi_bdt': BD,
         'llr': np.linspace(-1, 1, 21), 'pid_pi': BP, 'pid_p': BP, 'pid_mu': BP}
MUV = {'bragg_p': BR, 'bragg_mu': BR, 'bragg_mip': BR, 'bragg_pion': BR, 'mip_bdt': BD, 'muon_bdt': BD,
       'pid_mu': BP}
# physical ranges; anything else is a sentinel (Bragg-pion -1: branch absent; classifier -1: not evaluated)
VALID = {'bragg_p': (0., 500.), 'bragg_mu': (0., 500.), 'bragg_mip': (0., 500.), 'bragg_pion': (0., 500.),
         'mip_bdt': (-1., 1.), 'pi_bdt': (-1., 1.), 'muon_bdt': (-1., 1.), 'llr': (-1., 2.),
         'pid_pi': (0., 1.), 'pid_p': (0., 1.), 'pid_mu': (0., 1.)}
PV = sorted(set(POOLV) | {'ts', 'len', 'dist', 'contained', 'counted'})
MV = sorted(MUV)


def cuts(c):
    """{cut: (applicable, pass)} on the eligible candidates c (dict of numpy arrays)"""
    one = np.ones(len(c['llr']), bool)
    llr, mip, pib = c['llr'] > 0.1, c['mip_bdt'] > -0.1, c['pi_bdt'] > -0.1
    has = c['bragg_pion'] >= 0
    brg = c['bragg_pion'] >= 0.08
    return {'LLR > 0.1': (one, llr), 'MIP BDT > -0.1': (one, mip), 'pion BDT > -0.1': (one, pib),
            'Bragg-pion >= 0.08': (has, brg),
            'Bragg-pion >= 0.08, given the other three': (has & llr & mip & pib, brg),
            'released identification': (one, c['counted'] == 1),
            'released identification without Bragg-pion': (one, llr & mip & pib),
            'P(pi) > 0.20': (one, c['pid_pi'] > 0.20), 'P(pi) > 0.3087': (one, c['pid_pi'] > 0.3087)}


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def add_hist(H, key, x, w, bins):
    lo, hi = VALID[key[1]]
    ok = np.isfinite(x) & (x >= lo) & (x <= hi)
    xc = np.clip(x[ok], bins[0], bins[-1] - 1e-9)
    cur = H.get(key) or [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1), 0., 0., 0.]
    H[key] = [cur[0] + np.histogram(xc, bins, weights=w[ok])[0], cur[1] + np.histogram(xc, bins, weights=w[ok] ** 2)[0],
              cur[2] + float((x[ok] * w[ok]).sum()), cur[3] + float(w[ok].sum()), cur[4] + float(w[~ok].sum())]


def fill(path, scale, weighted, H, C):
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_sb_cc0pi', f'{S}_sb_pi0', f'{S}_sb_multipi', f'{S}_sb_cosmic'] + \
         [f'{D}_pdiag_{v}' for v in PV] + [f'{D}_mudiag_{v}' for v in MV] + \
         (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in t.iterate(br, step_size='200 MB', library='ak'):
        a = a[a[f'{S}_sb_cc0pi'] | a[f'{S}_sb_pi0'] | a[f'{S}_sb_multipi'] | a[f'{S}_sb_cosmic']]
        if len(a) == 0: continue
        w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) if weighted \
            else np.ones(len(a))
        w = w * scale
        # muon candidate
        m = ak.to_numpy(a[f'{D}_mudiag_muon_bdt']) > -998
        for v, bins in MUV.items():
            add_hist(H, ('muon', v), ak.to_numpy(a[f'{D}_mudiag_{v}'])[m], w[m], bins)
        # pion pool
        p = {v: a[f'{D}_pdiag_{v}'] for v in PV}
        ww = ak.broadcast_arrays(w, p['ts'])[0]
        pool = p['ts'] >= 0.5
        elig = pool & (p['contained'] == 1) & (p['dist'] <= 4) & (p['len'] > 20)
        for role, sel in (('pool', pool), ('eligible', elig)):
            wf = ak.to_numpy(ak.flatten(ww[sel]))
            for v, bins in POOLV.items():
                add_hist(H, (role, v), ak.to_numpy(ak.flatten(p[v][sel])), wf, bins)
        c = {v: ak.to_numpy(ak.flatten(p[v][elig])) for v in PV}
        wf = ak.to_numpy(ak.flatten(ww[elig]))
        for k, (app, ok) in cuts(c).items():
            cur = C[k] if k in C else np.zeros(4)
            C[k] = cur + np.array([(wf * (app & ok)).sum(), (wf ** 2 * (app & ok)).sum(), (wf * app).sum(), (wf ** 2 * app).sum()])


def frac(c):
    """pass fraction and its variance; weighted binomial, sum w^2 (pass - p)^2 / (sum w)^2"""
    pw, pw2, tw, tw2 = c
    if tw <= 0: return float('nan'), float('nan')
    f = pw / tw
    return f, ((1 - f) ** 2 * pw2 + f ** 2 * (tw2 - pw2)) / tw ** 2


def shape_chi2(d, d2, p, p2):
    Dn, Pn = d.sum(), p.sum()
    if Dn <= 0 or Pn <= 0: return float('nan'), 0
    r = Dn / Pn
    var = d2 + p2 * r * r
    m = var > 0
    return float(((d - p * r)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)


def main():
    mc, ext, data = F.samples()
    re = lambda path: NEW + os.path.basename(path)
    dirts = dirt_perrun('comb', where=lambda run: NEW)
    out = {}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        Hd, Hp, Cd, Cp = {}, {}, {}, {}
        for f in data[ip]: fill(re(f), 1.0, False, Hd, Cd)
        for f, s in mc[ip]: fill(re(f), s, True, Hp, Cp)
        for f, s in ext[ip]: fill(re(f), s, False, Hp, Cp)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, True, Hp, Cp)
        res = {'variables': {}, 'cuts': {}}
        for key in Hd:
            d1, d2, dsx, dn, dbad = Hd[key]; p1, p2, psx, pn, pbad = Hp[key]
            chi2, ndf = shape_chi2(d1, d2, p1, p2)
            res['variables'][f'{key[0]} | {key[1]}'] = dict(
                chi2=chi2, ndf=ndf, data_n=dn, pred_n=pn, data_mean=dsx / dn if dn else float('nan'),
                pred_mean=psx / pn if pn else float('nan'),
                data_invalid_frac=dbad / (dn + dbad) if dn + dbad else float('nan'),
                pred_invalid_frac=pbad / (pn + pbad) if pn + pbad else float('nan'),
                data_hist=list(d1), pred_hist=list(p1))
        for k in Cd:
            fd, vd = frac(Cd[k]); fp, vp = frac(Cp[k])
            res['cuts'][k] = dict(data=fd, data_err=np.sqrt(vd), pred=fp, pred_err=np.sqrt(vp), data_n=Cd[k][2],
                                  pull=(fd - fp) / np.sqrt(vd + vp) if vd + vp > 0 else float('nan'))
        out[name] = res
        print(f'== {name} done', flush=True)
    json.dump(out, open(os.path.join(OUTDIR, 'cr_pid_diag_check.json'), 'w'), indent=1, default=float)
    periods = list(out)
    head = ' | '.join(f'{p}{" (old)" if p in OLD else ""}' for p in periods)
    L = ['# Pion identification in the opened control regions: Bragg likelihoods, BDT scores and classifier', '',
         'Beam-on data against the prediction (overlay, beam-off, dirt), per period. Older-production periods: '
         + ', '.join(OLD) + '. Populations: pool = pion-pool candidates with track score >= 0.5; eligible = pool '
         'candidates the released pion identification can count (contained, vertex distance <= 4 cm, length > 20 cm); '
         'muon = the muon candidate.', '',
         '## Shape chi2/ndf (mean data - prediction; sentinel fraction data/prediction)', '',
         f'| population | variable | {head} |', '|---|---|' + '---|' * len(periods)]
    for k in sorted({k for r in out.values() for k in r['variables']}):
        role, var = k.split(' | ')
        cells = []
        for p in periods:
            r = out[p]['variables'].get(k)
            cells.append(f'{r["chi2"]:.1f}/{r["ndf"]} ({r["data_mean"] - r["pred_mean"]:+.3g}; inv {100 * r["data_invalid_frac"]:.0f}/{100 * r["pred_invalid_frac"]:.0f}%)'
                         if r and r['ndf'] > 0 else 'n/a')
        L.append(f'| {role} | {var} | ' + ' | '.join(cells) + ' |')
    L += ['', '## Pass fractions of the eligible candidates: data / prediction (pull)', '',
          f'| cut | {head} |', '|---|' + '---|' * len(periods)]
    for k in out[periods[0]]['cuts']:
        cells = []
        for p in periods:
            r = out[p]['cuts'][k]
            cells.append(f'{100 * r["data"]:.1f} / {100 * r["pred"]:.1f}% ({r["pull"]:+.1f})' if np.isfinite(r['data']) else
                         f'n/a / {100 * r["pred"]:.1f}%' if np.isfinite(r['pred']) else 'n/a')
        L.append(f'| {k} | ' + ' | '.join(cells) + ' |')
    open(os.path.join(OUTDIR, 'cr_pid_diag_check.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
