"""mp2pi_figures_cr.py -- figures of the two-pion selection against beam-on data in the opened control regions of the
single-pion analysis, and of the pseudo-experiment ensembles (report/other_notes/multipion_note.tex, sec:crdata).

  cr_mcs.pdf       MCS momentum of the pool tracks, data against the prediction before and after the MCS calibration
                   (cr_mp2pi_check.json and cr_mp2pi_check_mcscal.json of scripts/cr_mp2pi_check.py)
  cr_llr.pdf       LLR score of pool tracks at least 40 cm long: data, the calibrated prediction, and the prediction with
                   the overlay changed by the recombination variation (Recomb2 over its CV, per horn mode, from the
                   track dump of the detector-variation samples, mp_evt_llronly_mcscal/detvar); histograms cached in
                   /data/uboone/processed/mp2pi_cr_skim_mcscal/fig_llr_hist.pkl
  cr_score.pdf     event score of the two-pion domain (regions 1-4) and
  cr_ppi.pdf       summed P(pi) of the event: data against the calibrated prediction with the shape uncertainty of the
                   covariance of scripts/cr_mp2pi_syst_check.py (its cache syst_check_hist.pkl)
  ens_pulls.pdf    pulls of the ensembles, combined configuration, without and with the sideband constraint
  ens_altmodels.pdf  extracted over member truth for the alternative and stress models, per configuration and constraint
  ens_coverage.pdf   coverage of the quoted data-statistics interval against the confidence level (hashed-seed ensemble)

Every prediction is normalised to the data (shape comparison). The skims hold only events in the four control regions
and none in the signal region of either selection.

    python3 scripts/mp2pi_figures_cr.py [--only cr_mcs,cr_llr,...] [--ens nominal,nominal2]
        ->  report/figures/multipion/{cr_mcs,cr_llr,cr_score,cr_ppi,ens_pulls,ens_altmodels}.pdf
"""
import glob, json, os, pickle, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
PLAN = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
FIG = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'multipion'))
SKIM = '/data/uboone/processed/mp2pi_cr_skim_mcscal/'
DVT = '/data/uboone/processed/mp_evt_llronly_mcscal/detvar/'
LLR_CACHE = SKIM + 'fig_llr_hist.pkl'
SYST_CACHE = ('/data/uboone/processed/mp2pi_cr_skim_2f/' if '--prod' in sys.argv and sys.argv[sys.argv.index('--prod') + 1] == '2f'
              else '/data/uboone/processed/mp2pi_cr_skim_shr/') + 'syst_check_hist.pkl'
PROD = sys.argv[sys.argv.index('--prod') + 1] if '--prod' in sys.argv else ''   # 2f: ensembles of the two-fold production
ENS = f'/data/uboone/processed/mp2pi_{PROD or "mcscal"}/ens/'
FSUF = f'_{PROD}' if PROD else ''
LABEL = 'MicroBooNE NuMI data, control regions, internal'
ENS_NAME = {'nominal': 'consecutive seeds', 'nominal2': 'hashed seeds', 'alt545': 'GENIE universe 545',
            'altdelta': r'$\Delta\to N\pi$ angular', 'stress_band': 'stress: score-band signal +20%',
            'stress_mult': r'stress: $\geq3\pi$ +30%, $1\pi$ -20%'}
C_DATA, C_PRED, C_OLD, C_VAR = 'black', '#1f5fa8', '#c8553d', '#2a9d58'
C_SRC = {'overlay': '#8fb3de', 'beam-off': '#e9c46a', 'dirt': '#b5838d'}

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13, 'legend.fontsize': 10.5, 'xtick.labelsize': 11,
                     'ytick.labelsize': 11, 'axes.linewidth': 0.9, 'pdf.fonttype': 42})


def frame(ratio_label='data / pred.'):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(6.4, 4.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.25], hspace=0.06))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.93, bottom=0.14)
    fig.text(0.13, 0.945, LABEL, color='0.45', fontsize=9.5, ha='left', va='bottom')
    rx.axhline(1., color='0.4', lw=0.8, ls='--')
    rx.set_ylabel(ratio_label, fontsize=11)
    return fig, ax, rx


def steps(ax, edges, y, **kw):
    return ax.stairs(y, edges, **kw)


def data_points(ax, edges, d, **kw):
    c = 0.5 * (edges[1:] + edges[:-1])
    return ax.errorbar(c, d, yerr=np.sqrt(d), fmt='o', ms=3.5, color=C_DATA, lw=1, capsize=0, zorder=5, **kw)


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name), dpi=150)
    plt.close(fig)
    print('wrote', os.path.join(FIG, name))


# ----------------------------------------------------------------------------------------------------------- MCS
def fig_mcs():
    a = json.load(open(os.path.join(PLAN, 'cr_mp2pi_check.json')))['all']['track | mom_mcs']
    b = json.load(open(os.path.join(PLAN, 'cr_mp2pi_check_mcscal.json')))['all']['track | mom_mcs']
    edges = np.linspace(0, 2, 21)              # TRK['mom_mcs'] of scripts/cr_mp2pi_check.py, last bin with overflow
    d = np.array(b['data_hist']); assert np.allclose(d, a['data_hist'])
    pu, pc = (np.array(x['pred_hist']) * d.sum() / np.sum(x['pred_hist']) for x in (a, b))
    fig, ax, rx = frame()
    steps(ax, edges, pu, color=C_OLD, lw=1.6, label=f"prediction, uncalibrated ($\\chi^2$ = {a['chi2']:.0f}/{a['ndf']})")
    steps(ax, edges, pc, color=C_PRED, lw=1.8, label=f"prediction, MCS-calibrated ($\\chi^2$ = {b['chi2']:.1f}/{b['ndf']})")
    data_points(ax, edges, d, label=f'beam-on data ({d.sum():,.0f} tracks)')
    ax.set_ylabel('tracks / 0.1 GeV/$c$'); ax.set_ylim(0, 1.25 * max(d.max(), pu.max(), pc.max()))
    ax.legend(loc='upper right', frameon=False)
    ax.text(0.98, 0.52, 'pool tracks, all periods\nshape: prediction normalised to data\n$\\chi^2$: data and simulation statistics',
            transform=ax.transAxes, ha='right', va='top', fontsize=9.5, color='0.3')
    c = 0.5 * (edges[1:] + edges[:-1]); e = np.sqrt(d)
    for p, col, off in ((pu, C_OLD, -0.012), (pc, C_PRED, 0.012)):
        rx.errorbar(c + off, d / p, yerr=e / p, fmt='o', ms=3, color=col, lw=1, capsize=0)
    rx.set_ylim(0.2, 1.4); rx.set_xlim(0, 2)
    rx.set_xlabel('MCS momentum [GeV/$c$] (last bin: overflow)')
    save(fig, 'cr_mcs.pdf')


# ----------------------------------------------------------------------------------------------------------- LLR
LLR_EDGES = np.linspace(-1, 1, 41)
LONG = 40.


def llr_build():
    import uproot, awkward as ak
    from bkgfit import farsb as F
    from bkgfit.dirt import dirt_perrun
    B = 'CC1mu2piBDTEvt_etrk_'

    def h_tracks(path, scale, weighted):
        t = uproot.open(path)['stv_tree']
        br = [B + 'len', B + 'llr'] + (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
        h = np.zeros(len(LLR_EDGES) - 1); h2 = np.zeros_like(h)
        for a in t.iterate(br, step_size='200 MB', library='ak'):
            if len(a) == 0: continue
            w = ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']).astype(float) if weighted \
                else np.ones(len(a))
            w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.) * scale
            m = a[B + 'len'] >= LONG
            x = ak.to_numpy(ak.flatten(a[B + 'llr'][m])).astype(float)
            ww = ak.to_numpy(ak.flatten(ak.broadcast_arrays(w, a[B + 'llr'])[0][m]))
            x = np.clip(x, LLR_EDGES[0], LLR_EDGES[-1] - 1e-9)
            h += np.histogram(x, LLR_EDGES, weights=ww)[0]; h2 += np.histogram(x, LLR_EDGES, weights=ww ** 2)[0]
        return h, h2

    def h_detvar(horn, v):
        """long pool tracks of control-region events in the track dump of a detector-variation sample, CV-weighted"""
        BE = 'CC1mu2piEvt_'
        t = uproot.open(f'{DVT}xsec-ana-detvar_run4{horn}_{v}.root')['stv_tree']
        cr = [f'CC1mu1piXp_{c}' for c in ('sb_cc0pi', 'sb_multipi', 'sb_pi0', 'sb_cosmic')]
        h = np.zeros(len(LLR_EDGES) - 1)
        for a in t.iterate(cr + [BE + 'evt_presel', BE + 'etrk_len', BE + 'etrk_llr', 'tuned_cv_weight', 'ppfx_cv_weight',
                                 'normalisation_weight'], step_size='200 MB', library='ak'):
            m = np.zeros(len(a), bool)
            for c in cr: m |= ak.to_numpy(a[c]).astype(bool)
            a = a[m & ak.to_numpy(a[BE + 'evt_presel']).astype(bool)]
            if len(a) == 0: continue
            w = ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']).astype(float)
            w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)
            sel = a[BE + 'etrk_len'] >= LONG
            x = np.clip(ak.to_numpy(ak.flatten(a[BE + 'etrk_llr'][sel])).astype(float), LLR_EDGES[0], LLR_EDGES[-1] - 1e-9)
            h += np.histogram(x, LLR_EDGES, weights=ak.to_numpy(ak.flatten(ak.broadcast_arrays(w, a[BE + 'etrk_llr'])[0][sel])))[0]
        return h

    mc, ext, data = F.samples()
    dirts = dirt_perrun('comb', where=lambda run: SKIM)
    H = {h: {k: np.zeros(len(LLR_EDGES) - 1) for k in ('data', 'overlay', 'beam-off', 'dirt', 'w2')} for h in ('FHC', 'RHC')}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        Hh = H[name[:3]]
        for f in data[ip]: Hh['data'] += h_tracks(SKIM + os.path.basename(f), 1., False)[0]
        for kind, lst, wtd in (('overlay', mc[ip], True), ('beam-off', ext[ip], False)):
            for f, s in lst:
                h, h2 = h_tracks(SKIM + os.path.basename(f), s, wtd); Hh[kind] += h; Hh['w2'] += h2
        for r, f, s in dirts:
            if r == F.RUN[ip]:
                h, h2 = h_tracks(f, s, True); Hh['dirt'] += h; Hh['w2'] += h2
        print('==', name, flush=True)
    for horn in ('fhc', 'rhc'):
        cv, rc = h_detvar(horn, 'CV'), h_detvar(horn, 'Recomb2')
        H[horn.upper()]['ratio_Recomb2'] = np.where(cv > 0, (rc / rc.sum()) / np.where(cv > 0, cv / cv.sum(), 1.), 1.)
        H[horn.upper()]['detvar_cv'] = cv; H[horn.upper()]['detvar_Recomb2'] = rc
        print('== detvar', horn, flush=True)
    return H


def fig_llr():
    if os.path.exists(LLR_CACHE):
        H = pickle.load(open(LLR_CACHE, 'rb'))
    else:
        H = llr_build(); pickle.dump(H, open(LLR_CACHE, 'wb'))
    d = sum(H[h]['data'] for h in H)
    p = sum(H[h]['overlay'] + H[h]['beam-off'] + H[h]['dirt'] for h in H)
    # recombination variation: shape of Recomb2 over its CV applied to the overlay of each horn mode
    pr = sum(H[h]['overlay'] * H[h]['ratio_Recomb2'] + H[h]['beam-off'] + H[h]['dirt'] for h in H)
    pn, prn = p * d.sum() / p.sum(), pr * d.sum() / pr.sum()
    syst = json.load(open(os.path.join(PLAN, 'cr_mp2pi_syst_check.json')))['LLR score']
    fig, ax, rx = frame()
    steps(ax, LLR_EDGES, pn, color=C_PRED, lw=1.8, label='prediction (MCS-calibrated)')
    steps(ax, LLR_EDGES, prn, color=C_VAR, lw=1.6, ls='--', label='with recombination variation')
    data_points(ax, LLR_EDGES, d, label=f'beam-on data ({d.sum():,.0f})')
    sp = LLR_EDGES[1:] < -0.9 + 1e-9
    nd, npd = d[sp].sum(), p[sp].sum()
    ax.set_yscale('log'); ax.set_ylim(0.5, 3e3 * max(d.max(), pn.max()))
    ax.set_ylabel('tracks / 0.05')
    ax.legend(loc='upper left', frameon=False, fontsize=9.5)
    ax.text(0.98, 0.97, f'pool tracks $\\geq${LONG:.0f} cm, all periods\nprediction normalised to data\n'
            f'LLR $<-0.9$: data {nd:.0f}, pred. {npd:.0f}\n'
            f"covariance test: $\\chi^2$ = {syst['full'][0]:.1f}/{syst['full'][1]}",
            transform=ax.transAxes, ha='right', va='top', fontsize=9, color='0.3')
    c = 0.5 * (LLR_EDGES[1:] + LLR_EDGES[:-1]); e = np.sqrt(d)
    with np.errstate(divide='ignore', invalid='ignore'):
        rx.errorbar(c, np.where(pn > 0, d / pn, np.nan), yerr=np.where(pn > 0, e / pn, np.nan), fmt='o', ms=3, color=C_DATA, lw=1,
                    capsize=0, zorder=5)
    rx.stairs(np.where(pn > 0, prn / pn, np.nan), LLR_EDGES, color=C_VAR, lw=1.6, ls='--')
    rx.set_yscale('log'); rx.set_ylim(0.3, 30); rx.set_xlim(-1, 1)
    rx.set_ylabel('over pred.', fontsize=11)
    rx.set_xlabel('LLR score (proton $-1$, muon $+1$)')
    save(fig, 'cr_llr.pdf')
    return dict(data_below_m09=float(nd), pred_below_m09=float(npd), pred_norm_below_m09=float(pn[sp].sum()),
                recomb_norm_below_m09=float(prn[sp].sum()), norm=float(d.sum() / p.sum()))


# ------------------------------------------------------------------------------------------ score, summed P(pi)
def syst_pieces(q):
    import cr_mp2pi_syst_check as SC
    HP, R = pickle.load(open(SYST_CACHE, 'rb'))
    d = sum(HP[h][('data', q, 'cv')] for h in HP)
    src = {'overlay': sum(HP[h][('mc', q, 'cv')] for h in HP), 'beam-off': sum(HP[h].get(('ext', q, 'cv'), 0.) for h in HP),
           'dirt': sum(HP[h].get(('dirt', q, 'cv'), 0.) for h in HP)}
    p_h = {h: HP[h][('mc', q, 'cv')] + HP[h].get(('ext', q, 'cv'), 0.) + HP[h].get(('dirt', q, 'cv'), 0.) for h in HP}
    p = sum(p_h.values())
    w2 = sum(HP[h][('mc', q, 'w2')] + HP[h].get(('ext', q, 'w2'), 0.) + HP[h].get(('dirt', q, 'w2'), 0.) for h in HP)
    Cs = np.diag(w2); C = 0.
    for term in SC.UNIV:
        Pu = sum(HP[h][('mc', q, term)] + (p_h[h] - HP[h][('mc', q, 'cv')])[:, None] for h in HP)
        dev = np.stack([SC.shape_dev(Pu[:, u], p) for u in range(Pu.shape[1])], axis=1)
        C = C + dev @ dev.T / dev.shape[1]
    for v in SC.KNOBS:
        pv = sum(p_h[h] + HP[h][('mc', q, 'cv')] * (R[(h.lower(), v)][q] - 1.) for h in HP)
        dv = SC.shape_dev(pv, p); C = C + np.outer(dv, dv)
    return d, p, src, Cs, C, SC.test(d, p, Cs), SC.test(d, p, Cs + C), (SC.EVQ[q])


def fig_event(q, name, xlabel, ylab, xmax=None):
    d, p, src, Cs, C, stat, full, edges = syst_pieces(q)
    r = d.sum() / p.sum()
    fig, ax, rx = frame()
    bottom = np.zeros_like(p)
    for k in ('beam-off', 'dirt', 'overlay'):
        ax.stairs(bottom + r * src[k], edges, baseline=bottom, fill=True, color=C_SRC[k], alpha=0.9, label=k)
        bottom = bottom + r * src[k]
    sig = r * np.sqrt(np.diag(C + Cs))
    ax.stairs(r * p + sig, edges, baseline=r * p - sig, fill=True, facecolor='none', edgecolor='0.35', hatch='////', lw=0,
              label='shape uncertainty (syst. + MC stat.)')
    steps(ax, edges, r * p, color=C_PRED, lw=1.4)
    data_points(ax, edges, d, label=f'beam-on data ({d.sum():,.0f} events)')
    ax.set_ylabel(ylab); ax.set_ylim(0, 1.75 * max(d.max(), (r * p + sig).max()))
    h, l = ax.get_legend_handles_labels()
    order = [l.index(x) for x in l if x.startswith('beam-on')] + [l.index(x) for x in l if not x.startswith('beam-on')]
    ax.legend([h[i] for i in order], [l[i] for i in order], loc='upper right', frameon=False, ncol=1, fontsize=9.5)
    ax.text(0.97, 0.50, f'two-pion domain, score $\\leq$ 0.92\nprediction normalised to data ($\\times${r:.2f})\n'
            f'$\\chi^2$, full covariance: {full[0]:.1f}/{full[1]}\n$\\chi^2$, statistics only: {stat[0]:.1f}/{stat[1]}',
            transform=ax.transAxes, ha='right', va='top', fontsize=9.5, color='0.3')
    c = 0.5 * (edges[1:] + edges[:-1])
    rel = np.where(p > 0, sig / (r * p), 0.)
    rx.stairs(1 + rel, edges, baseline=1 - rel, fill=True, color='0.82', lw=0)
    rx.errorbar(c, np.where(p > 0, d / (r * p), np.nan), yerr=np.where(p > 0, np.sqrt(d) / (r * p), np.nan), fmt='o', ms=3.5,
                color=C_DATA, lw=1, capsize=0)
    rx.axhline(1., color='0.4', lw=0.8, ls='--')
    rx.set_ylim(0.5, 1.5); rx.set_xlim(edges[0], xmax if xmax else edges[-1]); rx.set_xlabel(xlabel)
    save(fig, name)
    return dict(stat=stat, full=full, norm=float(r))


# ------------------------------------------------------------------------------------------------------ ensembles
def ens_members(ens, cfg, mode):
    import mp2pi_ens_summary as ES
    dirs = sorted(glob.glob(f'{ENS}{ens}/m[0-9][0-9][0-9]'))
    return [m for m in (ES.member(x, cfg, mode) for x in dirs) if m]


def fig_pulls(names):
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 4.4), sharey=True)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.87, bottom=0.14, wspace=0.05)
    fig.text(0.11, 0.955, 'MicroBooNE NuMI simulation, pseudo-experiments, internal', color='0.45', fontsize=9.5)
    bins = np.linspace(-4, 4, 25); xs = np.linspace(-4, 4, 400); bw = bins[1] - bins[0]
    cols = ['#1f5fa8', '#c8553d', '#2a9d58']; out = {}
    for ax, (mode, title) in zip(axs, (('', 'no constraint'), ('_constr', 'sideband constraint'))):
        ymax = 0.
        for i, ens in enumerate(names):
            M = ens_members(ens, 'comb', mode)
            if len(M) < 5: continue
            x = np.array([m['x'] for m in M]); ex = np.array([m['expected'] for m in M]); st = np.array([m['stat'] for m in M])
            pull = (x - ex) / (st * x); mu, sd, n = pull.mean(), pull.std(ddof=1), len(pull)
            hh = ax.hist(np.clip(pull, bins[0], bins[-1] - 1e-9), bins, histtype='step', lw=1.6, color=cols[i],
                         label=f'{ENS_NAME.get(ens, ens)} ({n})\nmean {mu:+.2f}, width {sd:.2f}')[0]
            ax.plot(xs, n * bw * np.exp(-0.5 * ((xs - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)), color=cols[i], lw=1.1, ls='-.',
                    alpha=0.8)
            ymax = max(ymax, hh.max(), n * bw / np.sqrt(2 * np.pi) / min(sd, 1.))
            out[f'{ens} {mode or "none"}'] = dict(n=n, mean=float(mu), width=float(sd))
            if i == 0:
                ax.plot(xs, n * bw * np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi), color='0.25', lw=1.2, ls='--',
                        label='unit Gaussian')
        ax.set_title(f'combined, {title}', fontsize=11)
        ax.set_ylim(0, 1.75 * ymax); ax.legend(loc='upper left', frameon=False, fontsize=8.8)
        ax.set_xlim(-4, 4); ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    axs[0].set_ylabel('pseudo-experiments')
    fig.supxlabel('(extracted $-$ expected) / quoted data statistics', fontsize=12, y=0.02)
    save(fig, f'ens_pulls{FSUF}.pdf')
    return out


def fig_coverage(ens='nominal2'):
    """coverage of the quoted data-statistics interval against the confidence level, FHC, RHC and combined, without and
    with the constraint, with the binomial band of the ensemble size"""
    from scipy.stats import norm
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.6), sharey=True)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.86, bottom=0.16, wspace=0.05)
    fig.text(0.11, 0.945, f'MicroBooNE NuMI simulation, pseudo-experiments ({ENS_NAME.get(ens, ens)}), internal',
             color='0.45', fontsize=9)
    cl = np.linspace(0.02, 0.995, 120); out = {}
    for ax, (mode, title) in zip(axs, (('', 'no constraint'), ('_constr', 'sideband constraint'))):
        n_ = 0
        for cfg, nm, col in (('fhc5', 'FHC', '#c8553d'), ('rhcfull', 'RHC', '#2a9d58'), ('comb', 'combined', '#1f5fa8')):
            M = ens_members(ens, cfg, mode)
            if len(M) < 5: continue
            x = np.array([m['x'] for m in M]); ex = np.array([m['expected'] for m in M]); st = np.array([m['stat'] for m in M])
            pull = (x - ex) / (st * x); n_ = len(pull)
            cov = np.array([np.mean(np.abs(pull) < norm.ppf(0.5 + c / 2)) for c in cl])
            ax.plot(100 * cl, 100 * cov, color=col, lw=1.6, label=f'{nm}: width {pull.std(ddof=1):.2f}')
            out[f'{cfg} {mode or "none"}'] = {f'{c:g}': float(np.mean(np.abs(pull) < norm.ppf(0.5 + c / 2))) for c in (0.68, 0.95)}
        band = np.sqrt(cl * (1 - cl) / max(n_, 1))
        ax.fill_between(100 * cl, 100 * (cl - band), 100 * (cl + band), color='0.85', lw=0, label=f'nominal $\\pm$ binomial ({n_})')
        ax.plot(100 * cl, 100 * cl, color='0.3', lw=1, ls='--')
        ax.set_title(title, fontsize=11); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.legend(fontsize=8.5, loc='upper left', frameon=False)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel('members inside the interval [%]')
    fig.supxlabel('confidence level of the quoted data-statistics interval [%]', fontsize=11, y=0.02)
    save(fig, f'ens_coverage{FSUF}.pdf')
    return out


def fig_altmodels():
    models = [m for m in ('alt545', 'altdelta', 'stress_band', 'stress_mult', 'stress_plead', 'stress_psoft')
              if len(glob.glob(f'{ENS}{m}/m[0-9][0-9][0-9]')) >= 5]
    lab = {'alt545': 'GENIE multisim universe 545', 'altdelta': '$\\Delta\\to N\\pi$ angular variation',
           'stress_band': 'stress: score-band signal +20%', 'stress_mult': 'stress: $\\geq3\\pi$ +30%, 1$\\pi$ $-$20%',
           'stress_plead': 'stress: signal $\\times(1+0.5\\sin p_\\pi^{\\rm lead})$',
           'stress_psoft': 'stress: signal $\\times(1+0.5\\cos\\pi p_\\pi^{\\rm soft}/0.5)$'}
    cats = [(c, m) for c in ('fhc5', 'rhcfull', 'comb') for m in ('', '_constr')]
    names = {'fhc5': 'FHC', 'rhcfull': 'RHC', 'comb': 'comb.'}
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True, gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.80, bottom=0.12)
    fig.text(0.14, 0.965, 'MicroBooNE NuMI simulation, pseudo-experiments, internal', color='0.45', fontsize=9.5)
    cols = ['#1f5fa8', '#c8553d', '#2a9d58', '#7b5ea7', '#e9a23b', '#6b6b6b']; out = {}
    bx.axhspan(-1, 1, color='0.88', lw=0, label='$\\pm1$ total uncertainty')
    for i, mdl in enumerate(models):
        off = (i - 0.5 * (len(models) - 1)) * 0.11
        for j, (c, m) in enumerate(cats):
            M = ens_members(mdl, c, m)
            if len(M) < 5: continue
            rt = np.array([q['x'] / q['truth'] for q in M]); tot = np.mean([q['total'] for q in M])
            mu, er = rt.mean(), rt.std(ddof=1) / np.sqrt(len(rt))
            kw = dict(fmt='o' if m else 's', ms=4.5, color=cols[i], lw=1.1, capsize=0, mfc=cols[i] if m else 'white')
            ax.errorbar(j + off, mu, yerr=er, label=lab[mdl] if j == 0 else None, **kw)
            bx.errorbar(j + off, (mu - 1) / tot, yerr=er / tot, **kw)
            out[f'{mdl} {c}{m}'] = dict(n=len(rt), ratio=float(mu), err=float(er), bias_sigma=float((mu - 1) / tot), total=float(tot))
    ax.axhline(1., color='0.4', lw=0.8, ls='--'); bx.axhline(0., color='0.4', lw=0.8, ls='--')
    ax.text(0.01, 0.04, 'mean of each ensemble and its error; open: no constraint, filled: constrained',
            transform=ax.transAxes, fontsize=8.5, color='0.3')
    ax.set_ylabel('extracted / truth'); ax.set_ylim(0.86, 1.12)
    bx.set_ylabel('bias / total unc.'); bx.set_ylim(-1.3, 1.3)
    bx.set_xticks(range(len(cats))); bx.set_xticklabels([f"{names[c]}\n{'constr.' if m else 'none'}" for c, m in cats], fontsize=10.5)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.42), ncol=2, frameon=False, fontsize=8.5)
    bx.legend(loc='lower right', frameon=False, fontsize=9)
    save(fig, f'ens_altmodels{FSUF}.pdf')
    return out


def main():
    only = sys.argv[sys.argv.index('--only') + 1].split(',') if '--only' in sys.argv else None
    if '--ens' in sys.argv:
        ens = sys.argv[sys.argv.index('--ens') + 1].split(',')
    else:
        ens = ['nominal'] + (['nominal2'] if len(ens_members('nominal2', 'comb', '')) >= 20 else [])
    run = lambda k: only is None or k in only
    res = {}
    if run('cr_mcs'): fig_mcs()
    if run('cr_llr'): res['cr_llr'] = fig_llr()
    if run('cr_score'): res['cr_score'] = fig_event('event score', f'cr_score{FSUF}.pdf', 'event score', 'events / 0.077')
    if run('cr_ppi'): res['cr_ppi'] = fig_event('summed P(pi)', 'cr_ppi.pdf', 'summed $P(\\pi)$ of the event', 'events / 0.5', xmax=4.)
    if run('ens_pulls'): res['ens_pulls'] = fig_pulls(ens)
    if run('ens_coverage'): res['ens_coverage'] = fig_coverage()
    if run('ens_altmodels'): res['ens_altmodels'] = fig_altmodels()
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
