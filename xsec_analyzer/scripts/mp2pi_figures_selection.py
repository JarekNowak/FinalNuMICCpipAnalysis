"""mp2pi_figures_selection.py -- figures of the two-pion selection (particle classifier, event classifier, working
point, regions) for the slides and the note, from the dumps and results of the Phase 1-4 studies.

  pid_roc             pion-vs-rest ROC on the pion-candidate pool of the held-out runs (FHC Run 5, RHC Run 3): every
                      track with track score >= 0.5 and length >= 5 cm that is not the muon candidate and lies within
                      9.5 cm of the vertex, unmatched tracks labelled "other" (the population of
                      phase1_classifier_unmatched_common.json); the LLR-and-geometry classifier, the 43-input
                      classifier, the first-draft pion BDT (split at 20 cm), the pion BDT of the single-pion analysis
                      and the LLR score (dumps /data/uboone/processed/mp_pid/pidtrk-*.root, models model_xgb_*.json)
  pid_probs           P(pi) of the LLR-and-geometry classifier on the same tracks, by true class
  evt_score           event score of the adopted event classifier (evt_clf_n2_llronly_pid) on the MCS-calibrated dump
                      (/data/uboone/processed/mp_evt_llronly_mcscal), two-pion domain, test periods projected to the
                      full exposure (combined), stacked by true class
  evt_effpur          efficiency and purity of the signal region against the score threshold, same events
  wp_scan             expected uncertainty of the one-bin combined total against the threshold, by term: without the
                      constraint (phase2_wp_llronly.json, flux 17% x (1 + B/S)) and with the four sidebands
                      (phase4_constraint_llronly_t*.json, the ConstrainedCalculator formalism)
  region_composition  predicted composition of the signal region and the four sidebands (phase4_constraint_llronly_t0.92)
  region_corr         correlation of the predicted counts of the five regions under the systematic covariance (flux,
                      cross section, reinteraction, detector; constraint_cov_llronly_t0.92.npz of
                      `scripts/mp_evt.py constraint --tag llronly --threshold 0.92 --dump-cov`)

    python3 scripts/mp2pi_figures_selection.py [name ...]   ->  report/figures/multipion/<name>.pdf
The track-level predictions are cached in /data/uboone/processed/mp_pid/figcache_pool_test.npz (--recompute to redo).
"""
import glob, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
FIGDIR = os.path.join(REPO, 'report', 'figures', 'multipion')
PLAN = os.path.join(REPO, 'report', 'planning', 'multipion')
PID = '/data/uboone/processed/mp_pid/'
EVT_CAL = '/data/uboone/processed/mp_evt_llronly_mcscal/'
EVT_ADOPT = '/data/uboone/processed/mp_evt_llronly/'
THR = 0.92

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 13, 'axes.labelsize': 14, 'axes.titlesize': 14, 'legend.fontsize': 11.5,
                     'xtick.labelsize': 12, 'ytick.labelsize': 12, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e',
                     'ytick.color': '#52514e', 'grid.color': '#e4e3df', 'grid.linewidth': 0.8, 'legend.frameon': False,
                     'pdf.fonttype': 42, 'lines.linewidth': 2.})
# categorical slots in fixed order (validated palette); the true classes keep their colour in every figure
SLOT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#8a8984'
CLASSES = [('Signal (CC2$\\pi^\\pm$)', SLOT[0]), ('CC1$\\pi^\\pm$', SLOT[1]), ('CC$\\geq$3$\\pi^\\pm$', SLOT[2]),
           ('CC$\\pi^0$', SLOT[3]), ('CC0$\\pi$', SLOT[4]), ('NC', SLOT[5]), ('Other $\\nu$, dirt', SLOT[6]),
           ('Beam-off', SLOT[7])]
# mp_evt.TOPO -> index in CLASSES (CC2pi outside the signal, CC other, nue CC and out-of-FV are "other")
TOPO_CLASS = {'CC2pi signal': 0, 'CC1pi': 1, 'CC3pi+': 2, 'CCpi0': 3, 'CC0pi': 4, 'NC': 5, 'CC2pi': 6, 'CCother': 6,
              'nueCC': 6, 'OOFV': 6, 'dirt': 6, 'beam-off': 7}
REGION_LABELS = ['Signal region', 'Three-pion-like', 'Proton-like', 'Score 0.5$-$0.92', 'Score < 0.5']


def tag(fig):
    fig.text(0.01, 0.99, 'MicroBooNE simulation, internal', ha='left', va='top', fontsize=11, color=MUTED)


def save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, name + '.pdf'), bbox_inches='tight', pad_inches=0.06)
    plt.close(fig)
    print(f'  wrote {name}.pdf', flush=True)


# ------------------------------------------------------------------ particle classifier
def pool_tracks():
    """P(mu, pi, p, other) of the two classifiers and the comparison scores on the pool of the held-out runs."""
    cache = PID + 'figcache_pool_test.npz'
    if os.path.exists(cache) and '--recompute' not in sys.argv:
        return dict(np.load(cache))
    import uproot, xgboost as xgb
    argv = sys.argv; sys.argv = [argv[0], '--common']
    import mp_pid_train as T                      # FEATURES of the 43-input models (--common at import)
    sys.argv = argv
    cols = sorted(set(T.FEATURES + T.OLD + T.META))
    parts = []
    for f in sorted(glob.glob(PID + 'pidtrk-*.root')):
        a = uproot.open(f)['trk'].arrays(['period'], library='np')
        if not np.isin(a['period'], T.TEST_PERIODS).any(): continue
        parts.append(uproot.open(f)['trk'].arrays(cols, library='np'))
        print(f'  {os.path.basename(f)}', flush=True)
    d = {c: np.concatenate([p[c] for p in parts]) for c in cols}
    lab = d['label'].copy(); lab[~(d['purity'] > 0.5)] = 3
    txt = open(os.path.join(REPO, 'xsec_analyzer', 'include', 'XSecAnalyzer', 'Selections', 'CC1mu2pi.hh')).read()
    dcut = float(re.search(r'pion_vtx_distance_cut\(\)\s*const\s*override\s*\{\s*return\s*([-0-9.]+)', txt).group(1))
    m = (np.isin(d['period'], T.TEST_PERIODS) & (d['ts'] >= 0.5) & (d['len'] >= 5.) & (d['is_mu_cand'] == 0)
         & (d['dist'] <= dcut))
    X = np.stack([d[c][m].astype(np.float32) for c in T.FEATURES], axis=1)
    X[~np.isfinite(X) | (X <= -9998) | (np.abs(X) > 1e30)] = T.MISSING
    out = dict(label=lab[m], horn=d['horn'][m], dist_cut=np.array(dcut))
    for key, model in (('P43', 'model_xgb_unmatched_common.json'), ('Pllr', 'model_xgb_unmatched_common_llronly.json')):
        b = xgb.Booster(); b.load_model(PID + model)
        out[key] = b.predict(xgb.DMatrix(X, feature_names=T.FEATURES, missing=np.nan)).astype(np.float32)
        print(f'  predicted {model}', flush=True)
    soft = d['len'][m] < 20.
    out['s_first'] = np.where(soft, d['s_mppi_soft'][m], d['s_mppi_hard'][m]).astype(np.float32)
    out['s_pion'] = d['s_pion'][m].astype(np.float32); out['llr'] = d['llr'][m].astype(np.float32)
    np.savez_compressed(cache, **out)
    return out


def roc(y, s):
    """Pion efficiency against the fraction of other tracks passing, and the area (ties split evenly)."""
    o = np.argsort(-s, kind='stable'); y = y[o]; s = s[o]
    last = np.r_[np.flatnonzero(np.diff(s)), len(s) - 1]      # one point per distinct score
    tp = np.cumsum(y)[last] / y.sum(); fp = np.cumsum(~y)[last] / (~y).sum()
    tp = np.r_[0., tp]; fp = np.r_[0., fp]
    return fp, tp, float(getattr(np, "trapezoid", None)(tp, fp) if hasattr(np, "trapezoid") else np.trapz(tp, fp))


def fig_pid_roc(T):
    y = T['label'] == 1
    curves = [('LLR and geometry (adopted)', T['Pllr'][:, 1], SLOT[0], '-'),
              ('43 inputs, calorimetric', T['P43'][:, 1], SLOT[1], '-'),
              ('First-draft pion BDT', T['s_first'], SLOT[2], '--'),
              ('Single-pion pion BDT', T['s_pion'], SLOT[3], '--'),
              ('LLR score', -T['llr'], SLOT[4], ':')]          # pions are muon-like: low LLR (proton high)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    res = {}
    for name, s, c, ls in curves:
        ok = np.isfinite(s); s = np.where(ok, s, -np.inf)
        fp, tp, a = roc(y, s)
        if a < 0.5:                                              # orientation of a score is a convention
            fp, tp, a = roc(y, -s)
        res[name] = a
        ax.plot(fp, tp, color=c, ls=ls, label=f'{name}  ({a:.3f})')
    ax.plot([0, 1], [0, 1], color='#c9c8c3', lw=1, zorder=0)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
    ax.set_xlabel('Other tracks passing (fraction)'); ax.set_ylabel('Pion efficiency')
    ax.grid(True); ax.legend(loc='lower right', title='ROC area', title_fontsize=10.5, alignment='left', fontsize=10.5,
                             bbox_to_anchor=(1.0, 0.0), borderaxespad=0.2)
    ax.set_title(f'Pion-candidate pool, FHC Run 5 + RHC Run 3: {len(y):,} tracks, {int(y.sum()):,} pions', loc='left',
                 fontsize=11, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); tag(fig); save(fig, 'pid_roc')
    return res


def fig_pid_probs(T):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    b = np.linspace(0, 1, 41)
    for k, (name, c, ls) in enumerate([('True muon', SLOT[1], '-'), ('True pion', SLOT[0], '-'),
                                       ('True proton', SLOT[2], '--'), ('Other (unmatched)', SLOT[7], ':')]):
        m = T['label'] == k
        ax.hist(T['Pllr'][m, 1], bins=b, histtype='step', density=True, color=c, ls=ls, lw=2,
                label=f'{name} ({100 * m.mean():.0f}%)')
    ax.set_yscale('log'); ax.set_xlim(0, 1)
    ax.set_xlabel('$P(\\pi)$, LLR-and-geometry classifier'); ax.set_ylabel('Tracks (area normalised)')
    ax.grid(True, axis='y'); ax.legend(loc='upper right', ncol=2, columnspacing=1.0, handlelength=1.6)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi * 400)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); tag(fig); save(fig, 'pid_probs')


# ------------------------------------------------------------------ event classifier
def event_scores():
    """Adopted event score on the MCS-calibrated dump; weights projected to the full exposure (test periods)."""
    argv = sys.argv; sys.argv = [argv[0]]
    os.environ['MP_EVT_DIR'] = EVT_CAL
    import mp_evt as M
    sys.argv = argv
    import xgboost as xgb
    D = dict(np.load(EVT_CAL + 'events.npz')); sums = json.load(open(EVT_CAL + 'sums.json'))
    pot = {int(r): v for r, v in sums['pot'].items()}
    names = M.feat_names(2, 'pid')
    b = xgb.Booster(); b.load_model(EVT_ADOPT + 'evt_clf_n2_llronly_pid.json')
    te = np.isin(D['run'], M.TEST_RUNS) & D['n2_d6'].astype(bool)
    sc = b.predict(xgb.DMatrix(np.stack([D[f'n2_{k}'][te] for k in names], 1), feature_names=names))
    horn = lambda r: M.FHC_RUNS if r in M.FHC_RUNS else M.RHC_RUNS
    proj = {r: sum(pot[x] for x in horn(r)) / pot[r] for r in M.TEST_RUNS}
    run = D['run'][te]
    w = D['w'][te] * np.vectorize(proj.get)(run)
    kind, sig, topo = D['kind'][te], D['sig2'][te].astype(bool), D['topology'][te]
    tc = np.array([TOPO_CLASS[t] for t in M.TOPO] + [6])          # topology 10: beam-off and dirt (no truth)
    cls = np.where(kind == 1, 7, np.where(kind == 2, 6, tc[np.clip(topo, 0, 10)]))
    cls = np.where(sig, 0, np.where(cls == 0, 6, cls))
    G = sum(sums['generated'][str(r)]['2pi'] * proj[r] for r in M.TEST_RUNS)
    return dict(score=sc, w=w, cls=cls, G=G)


def fig_evt_score(E):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    b = np.linspace(0, 1, 26)
    h = [np.histogram(E['score'][E['cls'] == k], bins=b, weights=E['w'][E['cls'] == k])[0] for k in range(len(CLASSES))]
    order = [7, 6, 5, 4, 3, 2, 1, 0]                          # signal on top
    bottom = np.zeros(len(b) - 1)
    for k in order:
        ax.bar(b[:-1], h[k], width=np.diff(b), bottom=bottom, align='edge', color=CLASSES[k][1], edgecolor='white',
               linewidth=0.6, label=CLASSES[k][0])
        bottom += h[k]
    ax.set_yscale('log'); ax.set_xlim(0, 1); ax.set_ylim(5, bottom.max() * 4)
    for x, t in ((0.5, '0.5'), (THR, '0.92')):
        ax.axvline(x, color=INK, ls='--', lw=1.2)
    ax.text(THR - 0.01, 0.95, 'signal region $\\rightarrow$', transform=ax.get_xaxis_transform(), ha='right', va='top',
            fontsize=11, color=INK2)
    ax.set_xlabel('Event score'); ax.set_ylabel('Events / 0.04')
    ax.set_title('Two-pion domain, combined (test periods projected to the full exposure)', loc='left', fontsize=11,
                 color=INK2)
    hd, lb = ax.get_legend_handles_labels(); hd, lb = hd[::-1], lb[::-1]; o = [0, 4, 1, 5, 2, 6, 3, 7]
    ax.legend([hd[i] for i in o], [lb[i] for i in o], ncol=4, loc='upper center', bbox_to_anchor=(0.45, -0.27),
              fontsize=10, columnspacing=0.8, handlelength=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); tag(fig); save(fig, 'evt_score')
    sr = E['score'] > THR
    S = E['w'][sr & (E['cls'] == 0)].sum(); N = E['w'][sr].sum()
    return dict(sr_events=float(N), sr_signal=float(S), purity=float(S / N), efficiency=float(S / E['G']))


def fig_evt_effpur(E):
    t = np.round(np.arange(0.50, 0.981, 0.005), 3)
    s, n = E['w'] * (E['cls'] == 0), E['w']
    eff = np.array([s[E['score'] > x].sum() for x in t]) / E['G']
    pur = np.array([s[E['score'] > x].sum() / n[E['score'] > x].sum() for x in t])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(t, 100 * pur, color=SLOT[0], label='Purity')
    ax.plot(t, 100 * eff, color=SLOT[1], label='Efficiency')
    ax.axvline(THR, color=INK, ls='--', lw=1.2)
    i = int(np.argmin(np.abs(t - THR)))
    for v, c, dy, va in ((pur[i], SLOT[0], 8, 'bottom'), (eff[i], SLOT[1], -8, 'top')):
        ax.plot(THR, 100 * v, 'o', color=c, ms=7, mec='white', mew=1.5)
        ax.annotate(f'{100 * v:.1f}%', (THR, 100 * v), xytext=(-8, dy), textcoords='offset points', ha='right',
                    va=va, fontsize=12, color=INK, bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none'))
    ax.set_xlim(0.5, 0.98); ax.set_ylim(0, 100)
    ax.set_xlabel('Event-score threshold'); ax.set_ylabel('Signal region (%)')
    ax.grid(True); ax.legend(loc='upper left')
    fig.tight_layout(rect=(0, 0, 1, 0.96)); tag(fig); save(fig, 'evt_effpur')


# ------------------------------------------------------------------ working point
TERMS = [('stat', 'data statistics', 'Data statistics', SLOT[0]), ('flux', 'flux', 'Flux', SLOT[1]),
         ('xsec', 'cross section', 'Cross section', SLOT[2]), ('reint', 'reinteraction', 'Reinteraction', SLOT[3]),
         ('detector', 'detector', 'Detector', SLOT[4])]


def fig_wp_scan():
    W = json.load(open(os.path.join(PLAN, 'phase2_wp_llronly.json')))
    t = np.array(W['thresholds']); c = W['sets']['pid']['COMB']
    keep = (t >= 0.80) & (t <= 0.965)
    files = sorted(glob.glob(os.path.join(PLAN, 'phase4_constraint_llronly_t0.*.json')))
    files = [f for f in files if re.search(r'_t0\.\d+\.json$', f)]
    tc, rows, none = [], [], []
    for f in files:
        r = json.load(open(f)); tc.append(r['threshold'])
        rows.append(r['constraint']['all sidebands']); none.append(r['constraint']['none'])
    o = np.argsort(tc); tc = np.array(tc)[o]; rows = [rows[i] for i in o]; none = [none[i] for i in o]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for key, _, lab, col in TERMS:
        a1.plot(t[keep], 100 * np.array(c[key])[keep], color=col, label=lab)
    a1.plot(t[keep], 100 * np.array(c['total'])[keep], color=INK, lw=3, label='Total')
    for key, k2, lab, col in TERMS:
        a2.plot(tc, [100 * r[k2] for r in rows], color=col, marker='o', ms=5, label=lab)
    a2.plot(tc, [100 * r['total'] for r in rows], color=INK, lw=3, marker='o', ms=6, label='Total')
    a1.set_ylim(0, 75); a2.set_ylim(0, 35)
    a1.set_title('No constraint', loc='left'); a2.set_title('Four sidebands', loc='left')
    for a in (a1, a2):
        a.axvline(THR, color=INK, ls='--', lw=1.2); a.grid(True)
        a.set_xlim(0.80 if a is a1 else 0.84, 0.965 if a is a1 else 0.965); a.set_xlabel('Event-score threshold')
    a1.set_ylabel('Uncertainty of the total (%)')
    i1 = int(np.argmin(np.abs(t - THR))); i2 = int(np.argmin(np.abs(tc - THR)))
    box = dict(boxstyle='round,pad=0.15', fc='white', ec='none')
    a1.annotate(f'{100 * c["total"][i1]:.1f}%', (THR, 100 * c['total'][i1]), xytext=(4, 10), textcoords='offset points',
                ha='left', fontsize=12, bbox=box)
    a2.annotate(f'{100 * rows[i2]["total"]:.1f}%', (THR, 100 * rows[i2]['total']), xytext=(5, -8),
                textcoords='offset points', ha='left', va='top', fontsize=12, bbox=box)
    h, l = a2.get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=6, bbox_to_anchor=(0.5, 0.94), fontsize=11.5, handlelength=1.5,
               columnspacing=1.2)
    fig.tight_layout(rect=(0, 0, 1, 0.86)); tag(fig); save(fig, 'wp_scan')
    return dict(none_at_092=c['total'][i1], constrained_at_092=rows[i2]['total'],
                constrained=dict(zip([float(x) for x in tc], [r['total'] for r in rows])))


# ------------------------------------------------------------------ regions
def fig_region_composition():
    r = json.load(open(os.path.join(PLAN, 'phase4_constraint_llronly_t0.92.json')))
    names = r['regions']; comp = r['composition']
    F = np.zeros((len(names), len(CLASSES)))
    for g, nm in enumerate(names):
        for k, v in comp[nm].items(): F[g, TOPO_CLASS[k]] += v
    N = F.sum(1); Ft = F / N[:, None]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    y = np.arange(len(names))[::-1]
    left = np.zeros(len(names))
    for k, (lab, col) in enumerate(CLASSES):
        ax.barh(y, Ft[:, k], left=left, height=0.62, color=col, edgecolor='white', linewidth=1.2, label=lab)
        for g in range(len(names)):
            if Ft[g, k] >= 0.09:
                ax.text(left[g] + Ft[g, k] / 2, y[g], f'{100 * Ft[g, k]:.0f}', ha='center', va='center', fontsize=10.5,
                        color='white' if k in (0, 5, 6, 7, 1) else INK)
        left += Ft[:, k]
    for g in range(len(names)):
        ax.text(1.02, y[g], f'{N[g]:,.0f}', ha='left', va='center', fontsize=11.5, color=INK)
    ax.text(1.02, len(names) - 0.35, 'events', ha='left', va='center', fontsize=10.5, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels(REGION_LABELS); ax.set_xlim(0, 1)
    ax.set_xlabel('Fraction of predicted events (%)'); ax.set_xticks([0, .25, .5, .75, 1])
    ax.set_xticklabels(['0', '25', '50', '75', '100'])
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    hd, lb = ax.get_legend_handles_labels(); o = [0, 4, 1, 5, 2, 6, 3, 7]      # rows read in class order
    ax.legend([hd[i] for i in o], [lb[i] for i in o], ncol=4, loc='upper center', bbox_to_anchor=(0.42, -0.2),
              fontsize=10, columnspacing=0.8, handlelength=1.0)
    fig.tight_layout(rect=(0, 0, 0.93, 0.96)); tag(fig); save(fig, 'region_composition')
    return dict(events=N.tolist(), signal=F[:, 0].tolist())


def fig_region_corr():
    f = EVT_ADOPT + 'constraint_cov_llronly_t0.92.npz'
    Z = np.load(f)
    keys = [k for k in Z.files if k.startswith('cov_') and (k in ('cov_flux', 'cov_reint', 'cov_detector') or k.startswith('cov_xsec'))]
    C = sum(Z[k] for k in keys)[:5, :5]
    R = C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
    from matplotlib.colors import LinearSegmentedColormap
    # every correlation is positive: one-hue sequential scale from 0 to 1
    cmap = LinearSegmentedColormap.from_list('seq', ['#f4f8fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf',
                                                      '#184f95', '#0d366b'])
    assert R.min() > 0., R.min()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    im = ax.imshow(R, cmap=cmap, vmin=0, vmax=1)
    for i in range(5):
        for j in range(5):
            ax.text(j, i, f'{R[i, j]:.2f}', ha='center', va='center', fontsize=12,
                    color='white' if R[i, j] > 0.6 else INK)
    lab = ['Signal\nregion', '3$\\pi$-like', 'p-like', 'Score\n0.5$-$0.92', 'Score\n$<$0.5']
    ax.set_xticks(range(5)); ax.set_xticklabels(lab, fontsize=11); ax.set_yticks(range(5)); ax.set_yticklabels(lab, fontsize=11)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.set_label('Correlation'); cb.outline.set_visible(False)
    ax.set_title('Flux, cross section, reinteraction, detector', fontsize=12, color=INK2, loc='left')
    fig.tight_layout(rect=(0, 0, 1, 0.96)); tag(fig); save(fig, 'region_corr')
    return dict(keys=keys, corr=R.round(3).tolist())


def main():
    want = [a for a in sys.argv[1:] if not a.startswith('--')] or ['pid_roc', 'pid_probs', 'evt_score', 'evt_effpur',
                                                                  'wp_scan', 'region_composition', 'region_corr']
    out = {}
    if {'pid_roc', 'pid_probs'} & set(want):
        T = pool_tracks()
        if 'pid_roc' in want: out['pid_roc'] = fig_pid_roc(T)
        if 'pid_probs' in want: fig_pid_probs(T)
    if {'evt_score', 'evt_effpur'} & set(want):
        E = event_scores()
        if 'evt_score' in want: out['evt_score'] = fig_evt_score(E)
        if 'evt_effpur' in want: fig_evt_effpur(E)
    if 'wp_scan' in want: out['wp_scan'] = fig_wp_scan()
    if 'region_composition' in want: out['region_composition'] = fig_region_composition()
    if 'region_corr' in want: out['region_corr'] = fig_region_corr()
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
