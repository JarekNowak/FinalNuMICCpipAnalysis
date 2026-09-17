"""bkgfit.studies -- the work plan of BACKGROUND_FIT_PLAN.md, runnable step by step.

  python3 -m bkgfit.studies composition     step 0a: class composition and nu/nubar fractions per beam
  python3 -m bkgfit.studies sensitivity     step 0b: Asimov Fisher information, identifiability, expected
                                                     gain on the one-bin total (conditional and fitted)
  python3 -m bkgfit.studies robustness [N]  step 0c: decorrelation variants, protocol-style propagation, leave-one-out closure
  python3 -m bkgfit.studies ensembles [N] [P] large toy ensembles (pull mean convergence, bias diagnosis), P processes
  python3 -m bkgfit.studies validate [N]    step 2 : toys, injections, signal injection, FHC-vs-RHC
                                                     parameter goodness of fit, A/B/MCMC cross-checks
  python3 -m bkgfit.studies data            exploratory fits to the beam-on control regions (POST-INSPECTION)
  python3 -m bkgfit.studies figures
Results: ../logs/bkgfit/<step>.json ; figures: ../report/figures/bkgfit/
"""
import json, os, sys, time
import numpy as np
from . import model as M
from . import fit as F
from . import gain as Gn

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.abspath(os.path.join(HERE, '..', '..', 'logs', 'bkgfit')); os.makedirs(LOG, exist_ok=True)
FIG = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'bkgfit')); os.makedirs(FIG, exist_ok=True)
FIT_BINS = M.CC0PI + M.PI0 + [M.COSMIC]
ALL = M.PERIODS

PARAMS = {
    'none': {},
    'CC0pi': {'CC0pi': [2, 3]},
    'CC0pi+pi0': {'CC0pi': [2, 3], 'pi0': [4, 5]},
    'CC0pi nu/nubar': {'CC0pi_nu': [2], 'CC0pi_nubar': [3]},
    'CC0pi+pi0 nu/nubar': {'CC0pi_nu': [2], 'CC0pi_nubar': [3], 'pi0_nu': [4], 'pi0_nubar': [5]},
    'all classes': {'CC0pi': [2, 3], 'pi0': [4, 5], 'multipi': [6, 7], 'other': [8, 9]},
}


def dump(name, obj):
    def conv(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.floating, np.integer)): return o.item()
        if isinstance(o, np.bool_): return bool(o)
        raise TypeError(type(o))
    with open(os.path.join(LOG, name + '.json'), 'w') as f: json.dump(obj, f, indent=1, default=conv)
    print('wrote', os.path.join(LOG, name + '.json'))


def T():
    if not hasattr(T, 'c'): T.c = M.Templates()
    return T.c


# --------------------------------------------------------------------------------------------- step 0a
def composition():
    t = T(); out = {}
    regions = {'SR': [M.SR], 'CC0pi': M.CC0PI, 'pi0': M.PI0, 'multipi': [M.MULTI], 'cosmic': [M.COSMIC]}
    for expo, per in (('full exposure', ALL), ('data periods', M.DATA_PERIODS)):
        for beam in ('FHC', 'RHC'):
            idx = [M.PERIODS.index(p) for p in per if p.startswith(beam)]
            row = {}
            for rn, bins in regions.items():
                cv = t.cv[idx][:, :, bins].sum((0, 2)); ext = t.ext[idx][:, bins].sum(); dirt = t.dirt[idx][:, bins].sum()
                tot = cv.sum() + ext + dirt
                row[rn] = dict(total=float(tot), ext=float(ext / tot), dirt=float(dirt / tot),
                               classes={M.CLASS_NAMES[c]: dict(frac=float((cv[2 * c] + cv[2 * c + 1]) / tot),
                                                                nubar=float(cv[2 * c + 1] / max(cv[2 * c] + cv[2 * c + 1], 1e-9)))
                                        for c in range(5)},
                               nubar_all=float(cv[1::2].sum() / cv.sum()))
            out[f'{expo} {beam}'] = row
    dump('composition', out)
    for k, row in out.items():
        print(f'\n{k}')
        for rn, r in row.items():
            cl = '  '.join(f"{n} {v['frac'] * 100:4.1f}% (nubar {v['nubar'] * 100:3.0f}%)" for n, v in r['classes'].items())
            print(f"  {rn:8s} N={r['total']:8.1f} EXT {r['ext'] * 100:4.1f}% nubar(all nu) {r['nubar_all'] * 100:4.1f}% | {cl}")
    return out


# --------------------------------------------------------------------------------------------- step 0b
def sensitivity():
    t = T(); out = {'fisher': {}, 'gain': {}}
    for expo, per in (('full exposure', ALL), ('data periods', M.DATA_PERIODS)):
        for beams in (('FHC',), ('RHC',), ('FHC', 'RHC')):
            periods = [p for p in per if p[:3] in beams]
            for pn, pars in PARAMS.items():
                if not pars: continue
                m = M.Model(t, periods, FIT_BINS, params=pars, merge=True)
                I = m.fisher(np.ones(len(m.names)))
                C = np.linalg.pinv(I); sig = np.sqrt(np.clip(np.diag(C), 0, None))
                corr = C / np.outer(sig, sig) if np.all(sig > 0) else C
                key = f"{expo} | {'+'.join(beams)} | {pn}"
                out['fisher'][key] = dict(names=m.names, sigma=sig, corr=corr, cond=float(np.linalg.cond(I)))
                print(f"{key:55s} " + '  '.join(f'{n}={s:.3f}' for n, s in zip(m.names, sig)) +
                      (f"  max|corr|={np.max(np.abs(corr - np.eye(len(sig)))):.2f}" if len(sig) > 1 else ''))
        # expected gain on the one-bin total
        G = Gn.build(t, per, FIT_BINS)
        free = {k: [v for v in pars.values()] for k, pars in PARAMS.items() if pars}
        s = Gn.summary(G, free)
        out['gain'][expo] = s
        print(f'\n== expected relative uncertainty on the one-bin total, {expo}')
        for b, row in s.items():
            print(f'  {b}: ' + '  '.join(f'{k}={v:.3f}' for k, v in row.items()))
    dump('sensitivity', out)
    return out


# --------------------------------------------------------------------------------------------- step 0c
def robustness(nclos=30):
    """Where the expected gain comes from, and whether it survives: decorrelation variants, the frozen
    procedure's propagation rule, and leave-one-out closure over flux universes, GENIE universes and knobs."""
    t = T(); out = {'variants': {}, 'closure': {}}
    per = ALL
    free = {k: [v for v in pars.values()] for k, pars in PARAMS.items() if pars}
    variants = [('A: all correlations (matched detector, noise subtracted)', dict()),
                ('A, detector POT-normalised, noise NOT subtracted', dict(det='pot')),
                ('A, detector as the sideband macros (unscaled) + independent noise', dict(det='raw')),
                ('A, flux decorrelated SR-CR', dict(decorrelate=('flux',))),
                ('A, detector decorrelated SR-CR', dict(decorrelate=('detector',))),
                ('A, flux and detector decorrelated', dict(decorrelate=('flux', 'detector'))),
                ('A, flux, detector, GENIE multisim decorrelated', dict(decorrelate=('flux', 'detector', 'xsec_multi')))]
    for name, kw in variants:
        G = Gn.build(t, per, FIT_BINS, **kw); sm = Gn.summary(G, {k: free[k] for k in ('CC0pi', 'CC0pi+pi0')})
        out['variants'][name] = {b: {k: v for k, v in row.items() if not k.startswith('prior ') or k == 'prior total'} for b, row in sm.items()}
        print(f"{name:50s} " + '  '.join(f"{b}: prior {r['prior total']:.3f} A {r['conditional (A)']:.3f} B(CC0pi+pi0) {r['fit CC0pi+pi0 (B)']:.3f}" for b, r in sm.items()))
    dn = {}
    for mode in ('run4fhc', 'run4rhc'):
        f = t.dvm[mode]['f']; C = t.dvm[mode]['noise']; bins = FIT_BINS + [M.TOT_CC0PI, M.TOT_PI0, M.SR]
        raw = np.sqrt((f[:, bins] ** 2).sum(0)); noise = np.sqrt(np.array([C[:, b, b].sum() for b in bins]))
        corr = np.sqrt(np.clip((f[:, bins] ** 2).sum(0) - np.array([C[:, b, b].sum() for b in bins]), 0, None))
        dn[mode] = dict(bins=[M.bin_label(b) for b in bins], shift_quad=raw, noise_quad=noise, corrected=corr)
        print(f"{mode}: detector shift (quadrature over knobs) per 2D bin: median {np.median(raw[:42]):.3f}, noise {np.median(noise[:42]):.3f}, corrected {np.median(corr[:42]):.3f}; "
              f"totals CC0pi {raw[-3]:.3f}/{noise[-3]:.3f}/{corr[-3]:.3f} pi0 {raw[-2]:.3f}/{noise[-2]:.3f}/{corr[-2]:.3f} SR {raw[-1]:.3f}/{noise[-1]:.3f}/{corr[-1]:.3f}")
    out['detector_noise'] = dn
    G = Gn.build(t, per, FIT_BINS)
    for fs in ('CC0pi', 'CC0pi+pi0'):
        ps = Gn.protocol_style(G, free[fs]); out['protocol-style ' + fs] = ps
        print(f"protocol-style propagation, {fs}: " + '  '.join(f"{b}: prior {r['prior']:.3f} -> {r['with_fit']:.3f} (sigma theta {np.round(r['sigma_theta'], 3)})" for b, r in ps.items()))
    rng = np.random.default_rng(7)
    for fam, idx in (('flux', rng.choice(600, nclos, replace=False)), ('xsec_multi', rng.choice(600, nclos, replace=False)),
                     ('detector', range(len(M.KNOBS)))):
        rows = Gn.closure(t, per, FIT_BINS, fam, [int(i) for i in idx])
        summ = {}
        for b in ('FHC', 'RHC'):
            a = {k: np.array([r[b][k] for r in rows]) for k in ('prior_rel', 'post_rel', 'prior_pull', 'post_pull')}
            summ[b] = {k: dict(mean=float(v.mean()), rms=float(np.sqrt((v ** 2).mean())), max=float(np.abs(v).max())) for k, v in a.items()}
        out['closure'][fam] = dict(rows=rows, summary=summ)
        print(f"closure {fam:10s} " + '  '.join(f"{b}: |bias| rms prior {s['prior_rel']['rms']:.3f} post {s['post_rel']['rms']:.3f}; pull rms prior {s['prior_pull']['rms']:.2f} post {s['post_pull']['rms']:.2f} (max {s['post_pull']['max']:.2f})" for b, s in summ.items()))
    dump('robustness', out)
    return out


# --------------------------------------------------------------------------------------------- step 2
BINSETS = {'totals': [M.TOT_CC0PI, M.TOT_PI0, M.COSMIC], '2D': FIT_BINS}


def _asimov(m, cs=None, beam_factor=None, knob=None, univ=None):
    cs = np.ones(M.NC) if cs is None else cs
    mc = np.einsum('gcb,c->gb', m.cv, cs)
    if univ is not None:
        s_, u = univ; mc = mc + m.univ[s_][u].sum(1)
    if knob is not None: mc = mc * (1 + m.dvfrac[:, knob, :])
    if beam_factor is not None:
        g, f = beam_factor; mc = mc.copy(); mc[g] = mc[g] * f
    return (mc + m.fixed).ravel()


def validate(ntoys=200):
    t = T(); rng = np.random.default_rng(20260916); out = {}
    per = M.DATA_PERIODS
    classes = lambda cl, f: np.where(np.isin(np.arange(M.NC), cl), f, 1.0)
    for bs, bins in BINSETS.items():
        o = out[bs] = {}
        mk = lambda pars: M.Model(t, per, bins, params=pars, merge=True)
        o['asimov'] = {}
        for pn, pars in PARAMS.items():
            if not pars: continue
            m = mk(pars); r = F.fit(m, m.predict(np.ones(len(m.names))), minos=True)
            o['asimov'][pn] = r; print(bs, 'asimov', pn, {k: round(v, 4) for k, v in r['values'].items()}, {k: round(v, 3) for k, v in r['errors'].items()})
        o['toys'] = {}
        for pn in ('CC0pi', 'CC0pi+pi0'):
            m = mk(PARAMS[pn]); th = np.ones(len(m.names)); pulls, cover, chi = [], [], []; t0 = time.time()
            for i in range(ntoys):
                d = F.toy(m, th, rng); r = F.fit(m, d)
                if not r['valid']: continue
                v = np.array([r['values'][n] for n in m.names]); e = np.array([r['errors'][n] for n in m.names])
                pulls.append((v - 1) / e); cover.append(np.abs(v - 1) < e); chi.append(r['chi2'])
            pulls = np.array(pulls); cover = np.array(cover); ndf = len(m.data) - len(m.names)
            o['toys'][pn] = dict(n=len(pulls), pull_mean=pulls.mean(0), pull_width=pulls.std(0), coverage=cover.mean(0),
                                chi2_mean=float(np.mean(chi)), ndf=ndf, pulls=pulls)
            print(f'{bs} toys {pn}: n={len(pulls)} pull mean {np.round(pulls.mean(0), 3)} width {np.round(pulls.std(0), 3)} '
                  f'coverage {np.round(cover.mean(0), 3)} <chi2> {np.mean(chi):.1f}/{ndf} ({time.time() - t0:.0f}s)')
        m = mk(PARAMS['CC0pi+pi0']); inj = {}
        cases = [('CC0pi x1.25', dict(cs=classes([2, 3], 1.25)), {'CC0pi': 1.25, 'pi0': 1.0}),
                 ('pi0 x0.80', dict(cs=classes([4, 5], 0.8)), {'CC0pi': 1.0, 'pi0': 0.8}),
                 ('signal x1.30', dict(cs=classes([0, 1], 1.3)), {'CC0pi': 1.0, 'pi0': 1.0}),
                 ('other x1.30 (not fitted)', dict(cs=classes([8, 9], 1.3)), {'CC0pi': 1.0, 'pi0': 1.0}),
                 ('RHC neutrino MC x1.15 only', dict(beam_factor=(1, 1.15)), {'CC0pi': 1.0, 'pi0': 1.0}),
                 ('detector WMX +1 sigma', dict(knob=M.KNOBS.index('WMX')), {'CC0pi': 1.0, 'pi0': 1.0}),
                 ('detector Recomb2 +1 sigma', dict(knob=M.KNOBS.index('Recomb2')), {'CC0pi': 1.0, 'pi0': 1.0})] \
            + [(f'flux universe {u}', dict(univ=('flux', u)), {'CC0pi': 1.0, 'pi0': 1.0}) for u in (0, 1, 2)] \
            + [(f'GENIE universe {u}', dict(univ=('xsec_multi', u)), {'CC0pi': 1.0, 'pi0': 1.0}) for u in (0, 1, 2)]
        for label, kw, truth in cases:
            d = _asimov(m, **kw); r = F.fit(m, d); g = F.gof(m, r, d)
            inj[label] = dict(fit=r, gof=g, truth=truth, bias_over_sigma={k: (r['values'][k] - truth[k]) / r['errors'][k] for k in m.names})
            print(f"{bs} inject {label:28s} " + '  '.join(f"{k}={r['values'][k]:.3f}+-{r['errors'][k]:.3f} (truth {truth[k]})" for k in m.names) + f"  chi2={r['chi2']:.2f}/{g['ndf']}")
        o['injections'] = inj
        o['beam_split'] = {}
        split = {'CC0pi_FHC': {'classes': [2, 3], 'beam': 'FHC'}, 'CC0pi_RHC': {'classes': [2, 3], 'beam': 'RHC'},
                 'pi0_FHC': {'classes': [4, 5], 'beam': 'FHC'}, 'pi0_RHC': {'classes': [4, 5], 'beam': 'RHC'}}
        for label, fR in (('asimov', 1.0), ('RHC neutrino MC x1.15 only', 1.15)):
            d = _asimov(m, beam_factor=(1, fR)) if fR != 1.0 else _asimov(m)
            lr = F.lr_test(m, M.Model(t, per, bins, params=split, merge=True), d)
            o['beam_split'][label] = lr; print(f"{bs} FHC-vs-RHC LR (shared vs beam-specific norms) [{label}]: delta chi2 {lr['delta_chi2']:.2f}/{lr['ndf']} p={lr['p']:.3g}")
        o['mcmc'] = {}
        for label, cs in (('asimov', np.ones(M.NC)), ('CC0pi x1.25', classes([2, 3], 1.25))):
            d = _asimov(m, cs=cs); r = F.fit(m, d, minos=True); mc = F.mcmc(m, d, n=4000 if bs == 'totals' else 2500, burn=1000, chains=4)
            o['mcmc'][label] = dict(minuit=r, mean=mc['mean'], std=mc['std'], q16=mc['q16'], q84=mc['q84'], rhat=mc['rhat'], ess=mc['ess'])
            np.save(os.path.join(LOG, f"mcmc_{bs}_{label.replace(' ', '_')}.npy"), mc['samples'])
            print(f"{bs} MCMC [{label}]: " + '  '.join(f"{k}: minuit {r['values'][k]:.3f}+-{r['errors'][k]:.3f} | mcmc {mc['mean'][k]:.3f}+-{mc['std'][k]:.3f} rhat {mc['rhat'][k]:.3f}" for k in m.names))
    dump('validate', out)
    return out


# --------------------------------------------------------------------------------------------- large ensembles
_ENS = {}


def _ens_worker(args):
    """One chunk of toys in a worker process. Toy variants isolate the source of any bias:
    'standard'  systematic throw from the covariance, then Poisson, fitted with the CNP variance (the 2026-09-16 default);
    'poisson'   Poisson around the nominal prediction, no systematic throw;
    'gauss'     systematic throw, then Gaussian statistical fluctuation (sqrt(mu)) instead of Poisson;
    'pearson'   standard toys, fitted with the Pearson data variance (mu at nominal) instead of CNP."""
    bs, pn, variant, seed, n = args
    if 'T' not in _ENS: _ENS['T'] = M.Templates()
    t = _ENS['T']; bins = BINSETS[bs]
    m = M.Model(t, M.DATA_PERIODS, bins, params=PARAMS[pn], merge=True, stat='pearson' if variant == 'pearson' else 'cnp')
    th = np.ones(len(m.names)); mu = m.predict(th); V = m.cov_syst(th)
    w, Q = np.linalg.eigh(V); L = Q * np.sqrt(np.clip(w, 0, None))
    rng = np.random.default_rng(seed); out = []
    for i in range(n):
        shift = np.zeros_like(mu) if variant == 'poisson' else L @ rng.standard_normal(len(w))
        lam = np.clip(mu + shift, 0, None)
        d = (lam + np.sqrt(lam) * rng.standard_normal(len(lam))) if variant == 'gauss' else rng.poisson(lam).astype(float)
        r = F.fit(m, d)
        if not r['valid']: out.append([np.nan] * (2 * len(th) + 1)); continue
        out.append([r['values'][k] for k in m.names] + [r['errors'][k] for k in m.names] + [r['chi2']])
    return np.array(out)


def ensembles(ntoys=10000, nproc=8, variants=('standard', 'poisson', 'gauss', 'pearson')):
    """Large toy ensembles for the pull mean, width and coverage, with the running mean against the number of toys."""
    from multiprocessing import Pool
    out = {}; chunk = 250
    jobs = []
    for bs in BINSETS:
        for pn in ('CC0pi', 'CC0pi+pi0'):
            for v in variants:
                n = ntoys if v == 'standard' else max(ntoys // 2, 1000)
                for c in range(0, n, chunk):
                    jobs.append((bs, pn, v, 1000003 * (len(jobs) + 1), min(chunk, n - c)))
    t0 = time.time()
    with Pool(nproc) as pool:
        res = pool.map(_ens_worker, jobs, chunksize=1)
    groups = {}
    for (bs, pn, v, _, _), r in zip(jobs, res): groups.setdefault((bs, pn, v), []).append(r)
    for (bs, pn, v), rs in groups.items():
        a = np.concatenate(rs); a = a[np.all(np.isfinite(a), axis=1)]; k = len(PARAMS[pn])
        vals, errs, chi = a[:, :k], a[:, k:2 * k], a[:, 2 * k]
        pulls = (vals - 1) / errs; n = len(pulls)
        steps = np.unique(np.round(np.logspace(2, np.log10(n), 25)).astype(int))
        run = [dict(n=int(s), mean=pulls[:s].mean(0).tolist(), err=(pulls[:s].std(0) / np.sqrt(s)).tolist()) for s in steps]
        # bias on theta itself and the average reported error, to separate a shift of theta from an error-scale effect
        out[f'{bs} | {pn} | {v}'] = dict(names=list(PARAMS[pn]), n=n, pull_mean=pulls.mean(0), pull_mean_err=pulls.std(0) / np.sqrt(n),
                                        pull_width=pulls.std(0), pull_width_err=pulls.std(0) / np.sqrt(2 * (n - 1)),
                                        coverage=(np.abs(vals - 1) < errs).mean(0), theta_mean=vals.mean(0), theta_rms=vals.std(0),
                                        err_mean=errs.mean(0), chi2_mean=float(chi.mean()), running=run)
        o = out[f'{bs} | {pn} | {v}']
        print(f"{bs:6s} {pn:10s} {v:8s} n={n:6d} pull mean " + ' '.join(f'{x:+.4f}+-{e:.4f}' for x, e in zip(o['pull_mean'], o['pull_mean_err']))
              + '  width ' + ' '.join(f'{x:.3f}' for x in o['pull_width']) + '  coverage ' + ' '.join(f'{x:.3f}' for x in o['coverage'])
              + '  <theta> ' + ' '.join(f'{x:.4f}' for x in o['theta_mean']) + f'  <chi2> {o["chi2_mean"]:.2f}')
    print(f'ensembles: {len(jobs)} chunks in {time.time() - t0:.0f}s')
    dump('ensembles', out)
    return out


# --------------------------------------------------------------------------------------------- data
def data():
    """POST-INSPECTION, EXPLORATORY: fits to the beam-on control regions already opened (cr_data_note).
    Not a release input and not a frozen-procedure response."""
    t = T(); out = {}
    allp = M.DATA_PERIODS; noR1 = [p for p in allp if p != 'FHC_R1']
    def report(key, m, r, g):
        out[key] = dict(fit=r, gof=g, labels=m.labels)
        print(f"{key:48s} chi2 {g['chi2']:8.2f}/{g['ndf']:3d} p={g['p']:.3g}  " + '  '.join(f"{k}={v:.3f}+-{r['errors'][k]:.3f}" for k, v in r['values'].items()))
    for bs, bins in BINSETS.items():
        for pername, per in (('all periods', allp), ('no FHC R1', noR1)):
            for pn in ('none', 'CC0pi', 'CC0pi+pi0', 'CC0pi+pi0 nu/nubar', 'all classes'):
                if bs == 'totals' and pn in ('CC0pi+pi0 nu/nubar', 'all classes'): continue
                m = M.Model(t, per, bins, params=PARAMS[pn], merge=True); r = F.fit(m, minos=bool(PARAMS[pn])); g = F.gof(m, r)
                report(f'{bs} | {pername} | {pn}', m, r, g)
                if PARAMS[pn]:
                    ri = F.fit_iterated(m); out[f'{bs} | {pername} | {pn}']['iterated'] = ri
            # likelihood-ratio tests against the unfitted prediction
            for pn in ('CC0pi', 'CC0pi+pi0'):
                lr = F.lr_test(M.Model(t, per, bins, merge=True), M.Model(t, per, bins, params=PARAMS[pn], merge=True))
                out[f'{bs} | {pername} | LR {pn} vs none'] = lr; print(f"   LR {pn} vs none: delta chi2 {lr['delta_chi2']:.2f}/{lr['ndf']} p={lr['p']:.3g}")
        # exposure hypothesis: per period, free FHC Run-1 factor
        for pn in ('none', 'CC0pi+pi0'):
            m = M.Model(t, allp, bins, params=PARAMS[pn], extra={'f_FHC_R1': ['FHC_R1']}, merge=False)
            r = F.fit(m, minos=True); g = F.gof(m, r); report(f'{bs} | per period + f_FHC_R1 | {pn}', m, r, g)
            prof = F.profile(m, 'f_FHC_R1', np.linspace(0.3, 1.2, 37)); out[f'{bs} | per period + f_FHC_R1 | {pn}']['profile'] = dict(grid=np.linspace(0.3, 1.2, 37), chi2=prof)
        # FHC vs RHC consistency (no FHC R1): shared against beam-specific normalisations, inside the joint covariance
        split = {'CC0pi_FHC': {'classes': [2, 3], 'beam': 'FHC'}, 'CC0pi_RHC': {'classes': [2, 3], 'beam': 'RHC'},
                 'pi0_FHC': {'classes': [4, 5], 'beam': 'FHC'}, 'pi0_RHC': {'classes': [4, 5], 'beam': 'RHC'}}
        lr = F.lr_test(M.Model(t, noR1, bins, params=PARAMS['CC0pi+pi0'], merge=True), M.Model(t, noR1, bins, params=split, merge=True))
        out[f'{bs} | beam-split LR (no FHC R1)'] = lr
        print(f"{bs} FHC-vs-RHC LR (no R1): delta chi2 {lr['delta_chi2']:.2f}/{lr['ndf']} p={lr['p']:.3g}; split " + str({k: round(v, 3) for k, v in lr['alt']['values'].items()}))
    # the frozen global normalisation test recomputed: as frozen (multisims only, unscaled detector ratios) and with
    # the corrected detector covariance (POT-normalised, event-matched, noise subtracted), same eight totals
    TOT8 = [M.TOT_CC0PI, M.MULTI, M.TOT_PI0, M.COSMIC]
    for (label, kw), (pername, per) in [(lk, pp) for pp in (('all periods', allp), ('no FHC R1', noR1)) for lk in (
                      ('as frozen: multisims, unscaled detector', dict(det='raw', dvstat=False, systs=['flux', 'xsec_multi', 'reint'], pot=False)),
                      ('multisims, corrected detector', dict(det='matched', systs=['flux', 'xsec_multi', 'reint'], pot=False)),
                      ('full model, corrected detector', dict(det='matched')))]:
        label = f'{label} | {pername}'
        m = M.Model(t, per, TOT8, merge=True, stat='neyman', **kw)
        mu = m.predict([]); V = m.cov_syst([]) + np.diag(np.maximum(m.data, 1.0)); r = m.data - mu
        chi = float(r @ np.linalg.solve(V, r)); Vsd = m.cov_detector(m.mc([])) + np.diag(np.maximum(m.data, 1.0)) + np.diag((np.einsum('gcb->gb', m.w2) + m.fixed2).ravel())
        pulls = r / np.sqrt(np.diag(V)); pulls_sd = r / np.sqrt(np.diag(Vsd))
        from scipy import stats as _st
        out[f'global test | {label}'] = dict(labels=m.labels, data=m.data, pred=mu, chi2=chi, p=float(_st.chi2.sf(chi, len(r))), pull_full=pulls, pull_statdet=pulls_sd,
                                            det_frac=np.sqrt(np.diag(m.cov_detector(m.mc([])))) / mu)
        print(f"global test [{label}]: chi2 {chi:.2f}/{len(r)} p={_st.chi2.sf(chi, len(r)):.3f}; pulls full {np.round(pulls, 2)}; st+det {np.round(pulls_sd, 2)}; detector frac {np.round(np.sqrt(np.diag(m.cov_detector(m.mc([])))) / mu, 3)}")
    # pre-fit and post-fit (totals fit) residuals in the 2D bins, no FHC R1
    m2 = M.Model(t, noR1, FIT_BINS, params=PARAMS['CC0pi+pi0'], merge=True)
    rt = out['totals | no FHC R1 | CC0pi+pi0']['fit']; th = np.array([rt['values'][k] for k in m2.names])
    mu = m2.predict(th); V = m2.cov_syst(th) + m2.cov_stat(mu, m2.data)
    out['residuals'] = dict(labels=m2.labels, data=m2.data, prefit=m2.predict(np.ones(2)), postfit=mu, sigma=np.sqrt(np.diag(V)))
    # conditional constraint of the signal-region total prediction, from totals and from 2D bins (model-only; no SR data)
    for bs, bins in BINSETS.items():
        ball = bins + [M.SR]; mall = M.Model(t, noR1, ball, merge=True); nb = len(ball)
        cr = [g * nb + i for g in range(len(mall.groups)) for i in range(nb - 1)]; sr = [g * nb + nb - 1 for g in range(len(mall.groups))]
        meanA, covA, _ = F.conditional(mall, cr, sr, mall.data[cr])
        pri = mall.predict([])[sr]; psig = np.sqrt(np.diag(mall.cov_syst([]))[sr])
        out[f'conditional SR total ({bs}, no FHC R1)'] = dict(prior=pri, prior_sigma=psig, posterior=meanA, posterior_sigma=np.sqrt(np.diag(covA)))
        print(f'conditional SR total prediction from {bs} (FHC, RHC): prior {np.round(pri, 1)} +- {np.round(psig, 1)} -> posterior {np.round(meanA, 1)} +- {np.round(np.sqrt(np.diag(covA)), 1)}')
    dump('data', out)
    return out


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'composition'
    if step == 'validate': validate(int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    elif step == 'robustness': robustness(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
    elif step == 'ensembles': ensembles(int(sys.argv[2]) if len(sys.argv) > 2 else 10000, int(sys.argv[3]) if len(sys.argv) > 3 else 8)
    else: globals()[step]()
