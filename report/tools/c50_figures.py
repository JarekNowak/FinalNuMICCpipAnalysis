# c50_figures.py -- figures for the 2026-09-16 studies (0.50 criterion, proton angles, 2D pairs).
# Everything here is PRE-EXTRACTION: migration diagonals, response matrices, angular resolution by
# proton-candidate truth, and generator predictions. No unfolded cross sections.
#   python3 report/tools/c50_figures.py            -> report/figures/c50/*.png
# Event arrays are cached in logs/c50/cache/*.npz (delete to re-read the ntuples).
import os, csv, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
LOG = os.path.join(REPO, 'logs', 'c50'); CACHE = os.path.join(LOG, 'cache')
OUT = os.path.join(REPO, 'report', 'figures', 'c50'); os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
GEN = os.path.join(REPO, 'generator_predictions', 'gen2d')
PI = 3.1416

# reference categorical palette (fixed order), text inks, recessive grid
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
MARK = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BLUES = LinearSegmentedColormap.from_list('blues', ['#fcfcfb', '#cde2fb', '#86b6ef', '#2a78d6', '#184f95', '#0d366b'])
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'text.color': INK, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'font.size': 10, 'legend.frameon': False,
                     'lines.linewidth': 2})

def crit_lines(ax):
    ax.axhline(0.68, color=INK, lw=1, ls='--'); ax.axhline(0.50, color=INK2, lw=1, ls=':')
    x1 = ax.get_xlim()[1]
    ax.text(x1, 0.68, ' 0.68', va='center', ha='left', color=INK, fontsize=8, clip_on=False)
    ax.text(x1, 0.50, ' 0.50', va='center', ha='left', color=INK2, fontsize=8, clip_on=False)

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight'); plt.close(fig)
    print('wrote', name)

# ------------------------------------------------------------------------------------------------
# 1. migration diagonals per true bin, from the screen TSVs
def read_tsv(mode, sample):
    rows = []
    for f in [f'screen_{mode}_{sample}.tsv', f'screen_{mode}_{sample}_extra.tsv']:
        p = os.path.join(LOG, f)
        if os.path.exists(p): rows += list(csv.DictReader(open(p), delimiter='\t'))
    return rows

def diag_panels(mode, sample, panels, name, title, ncol):
    rows = read_tsv(mode, sample)
    n = len(panels); nrow = (n + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.1 * nrow), squeeze=False)
    for ax, (obs, labels, ptitle) in zip(axs.flat, panels):
        for k, lab in enumerate(labels):
            r = [x for x in rows if x['obs'] == obs and x['label'] == lab]
            if not r: continue
            y = [float(x['diag']) for x in r]; xs = np.arange(1, len(y) + 1)
            released = lab.startswith('released')
            ax.plot(xs, y, marker=MARK[k], ms=6, color=INK2 if released else SER[k % 8],
                    ls='-' if released else '--', lw=2 if released else 1.5, label=lab)
        ax.set_ylim(0, 1.02); ax.set_title(ptitle, fontsize=10, loc='left')
        ax.set_xlabel('true bin'); ax.set_ylabel('migration diagonal')
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        crit_lines(ax); ax.legend(fontsize=7, loc='lower left')
    for ax in list(axs.flat)[n:]: ax.set_visible(False)
    fig.suptitle(title, x=0.01, ha='left', fontsize=12)
    fig.tight_layout(); save(fig, name)

def fig_diagonals():
    for mode in ('fhc', 'rhc'):
        M = mode.upper()
        diag_panels(mode, 'incl', [
            ('pmu', ['released 7', 'cand 8', 'cand 9'], 'p_mu'),
            ('costhmu', ['released 5', 'cand 6', 'cand 7', 'cand 8'], 'cos theta_mu'),
            ('thetamu', ['released 5', 'cand 7', 'cand 8'], 'theta_mu'),
            ('costhpi', ['released 5 (fwd split)', 'cand 6', 'cand 6 fwd', 'cand 7'], 'cos theta_pi'),
            ('thmupi', ['released 5', 'cand 6', 'cand 7', 'cand 8'], 'theta_mu-pi'),
            ('ppi', ['released 2', 'cand 3', 'cand 3b', 'cand 4'], 'p_pi')],
            f'diag_incl_{mode}', f'Inclusive: column-normalised migration diagonal per true bin, {M}', 3)
        diag_panels(mode, '1p', [
            ('deltaPt', ['released 2', 'c50 3'], 'delta p_T'),
            ('dalphat', ['released 2', 'c50 3'], 'delta alpha_T'),
            ('dphit', ['released 2', 'c50 4'], 'delta phi_T'),
            ('pn', ['released 2', 'c50 4'], 'p_n'),
            ('Wpipr', ['released 2', 'c50 2'], 'W_pi-p'),
            ('costhp', ['c68 4', 'c50 6'], 'cos theta_p'),
            ('thetap', ['c68 4', 'c68 5', 'c50 7'], 'theta_p (4 bins built)'),
            ('thpipr', ['c68 3', 'c50 7'], 'theta_pi-p (3 bins built)')],
            f'diag_1p_{mode}', f'Proton-tagged: column-normalised migration diagonal per true bin, {M}', 4)

# 2D: stacked diagonal / X-leak / Y-leak per flattened bin
PAIRS2D = [('incl', 'costhpi|costhmu', '3x2', 'cos theta_pi | cos theta_mu', 2),
           ('incl', 'costhmu|pmu', '2x2', 'cos theta_mu | p_mu', 2),
           ('incl', 'costhpi|thmupi', '3x2', 'cos theta_pi | theta_mu-pi', 2),
           ('incl', 'costhpi|ppi', '2x2', 'cos theta_pi | p_pi', 2),
           ('1p', 'thetap|deltaPt', '2x2', 'theta_p | delta p_T', 2),
           ('1p', 'thpipr|deltaPt', '2x2', 'theta_pi-p | delta p_T', 2),
           ('1p', 'thetap|Wpipr', '2x2', 'theta_p | W_pi-p', 2),
           ('1p', 'thpipr|Wpipr', '2x2', 'theta_pi-p | W_pi-p', 2)]

def fig_2d_leak():
    for mode in ('fhc', 'rhc'):
        fig, axs = plt.subplots(2, 4, figsize=(16, 6.4))
        for ax, (sample, obs, lab, title, nx) in zip(axs.flat, PAIRS2D):
            r = [x for x in read_tsv(mode, sample) if x['obs'] == obs and x['label'] == lab]
            d = np.array([float(x['diag']) for x in r]); xl = np.array([float(x['xleak']) for x in r])
            yl = np.array([float(x['yleak']) for x in r]); n = np.array([float(x['sel_sig']) for x in r])
            sl = np.array([int(x['slice']) for x in r]); xs = np.arange(1, len(d) + 1)
            gap = dict(edgecolor=SURF, linewidth=2)
            ax.bar(xs, d, color=SER[0], label='diagonal', **gap)
            ax.bar(xs, xl, bottom=d, color=SER[1], label='X-leak (same Y slice)', **gap)
            ax.bar(xs, yl, bottom=d + xl, color=SER[2], label='Y-leak (other Y slice)', **gap)
            for b in range(1, len(d)):
                if sl[b] != sl[b - 1]: ax.axvline(b + 0.5, color=INK2, lw=1)
            for x, v, nn in zip(xs, d, n):
                ax.text(x, 0.02, f'{nn:.0f}', ha='center', va='bottom', color=SURF, fontsize=7)
            ax.set_ylim(0, 1.05); ax.set_xticks(xs); ax.set_title(f'{title}  ({lab})', fontsize=10, loc='left')
            ax.set_xlabel('true bin (Y slice | X bin)'); ax.set_ylabel('fraction of selected signal')
            crit_lines(ax)
        h, l = axs.flat[0].get_legend_handles_labels()
        fig.suptitle(f'2D pairs: where selected signal is reconstructed, per true bin, {mode.upper()} '
                     '(numbers in bars: expected selected signal at data POT)', x=0.01, ha='left', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.legend(h, l, loc='upper left', bbox_to_anchor=(0.01, 0.955), ncol=3, fontsize=9)
        save(fig, f'diag_2d_{mode}')

# ------------------------------------------------------------------------------------------------
# event-level arrays (same weights and per-run scales as macros/binning_screen_2d.C)
FILES = {
    'fhc': [('Run1_fhc_new_numi_flux_fhc_pandora_ntuple', 0.14101), ('Run2_fhc_new_numi_flux_fhc_pandora_ntuple', 0.05085),
            ('Run4_fhc_new_numi_flux_fhc_pandora_ntuple', 0.07323), ('reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc', 0.11560)],
    'rhc': [(f'{r}_rhc_new_numi_flux_rhc_pandora_ntuple', s) for r, s in
            [('Run1', 0.06728), ('Run2', 0.04478), ('Run4a', 0.08847), ('Run4b', 0.08847), ('Run4c', 0.08847)]]
           + [(f'Run3_rhc_new_numi_flux_rhc_pandora_ntuple_{x}', 0.09066) for x in ('aa', 'ab', 'ac', 'ad', 'ae')]}

def load(mode, sample):
    import uproot
    cache = os.path.join(CACHE, f'{mode}_{sample}.npz')
    if os.path.exists(cache): return dict(np.load(cache))
    S = 'CC1mu1pi1p' if sample == '1p' else 'CC1mu1piXp'
    base = '/data/uboone/processed/' + ('w/' if sample == '1p' else '')
    names = ['Selected', 'MC_Signal', 'candidate_muon_mom_true', 'candidate_muon_mom_reco', 'candidate_muon_costh_true',
             'candidate_muon_costh_reco', 'candidate_pion_costh_true', 'candidate_pion_costh_reco',
             'true_mu_pi_opening_angle', 'mu_pi_opening_angle', 'candidate_pion_mom_true', 'candidate_pion_mom_reco']
    if sample == '1p':
        names += ['proton_costh_true', 'proton_costh_reco', 'deltaPt_true', 'deltaPt_reco', 'W_pipr_true', 'W_pipr_reco']
    out = {}
    for fn, scale in FILES[mode]:
        path = f'{base}xsec-ana-{fn}.root'
        t = uproot.open(path)['stv_tree']
        a = t.arrays([f'{S}_{n}' for n in names] + ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'], library='np')
        keep = a[f'{S}_Selected'].astype(bool) | a[f'{S}_MC_Signal'].astype(bool)
        cv = a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']
        cv = np.where(np.isfinite(cv) & (cv >= 0) & (cv <= 30), cv, 1.0)
        chunk = {n: a[f'{S}_{n}'][keep].astype(float) for n in names}
        chunk['w'] = (cv * scale)[keep]
        if sample == '1p':
            f = uproot.open(f'{base}friends_pa/xsec-ana-{fn}.pa.root')['pa'].arrays(['thpipr_reco', 'thpipr_true', 'prmatch'], library='np')
            for k in f: chunk[k] = f[k][keep].astype(float)
        for k, v in chunk.items(): out.setdefault(k, []).append(v)
        print('  read', fn, keep.sum())
    out = {k: np.concatenate(v) for k, v in out.items()}
    np.savez(cache, **out); return out

def binidx(v, e, top_open):
    e = np.asarray(e); b = np.searchsorted(e, v, side='right') - 1
    n = len(e) - 1
    b = np.where(v < e[0], -1, b)
    b = np.where(b >= n, n - 1 if top_open else -1, b)
    return b

def ppi_reco(p):
    return np.sqrt((np.sqrt(p ** 2 + 0.011164) - 0.10566 + 0.13957) ** 2 - 0.019480)

def response(tb, rb, w, nb):
    """column-normalised (true columns), normalised to ALL selected signal in the true bin"""
    M = np.zeros((nb, nb)); den = np.zeros(nb)
    for i in range(nb):
        m = tb == i; den[i] = w[m].sum()
        for j in range(nb): M[j, i] = w[m & (rb == j)].sum()
    return np.divide(M, den, out=np.zeros_like(M), where=den > 0)

def heat(ax, R, title, ticks=None, seps=()):
    nb = R.shape[0]
    ax.imshow(R, origin='lower', cmap=BLUES, vmin=0, vmax=1, extent=(0.5, nb + 0.5, 0.5, nb + 0.5))
    for i in range(nb):
        for j in range(nb):
            ax.text(i + 1, j + 1, f'{R[j, i]:.2f}', ha='center', va='center', fontsize=7 if nb > 5 else 8,
                    color=SURF if R[j, i] > 0.55 else INK)
    for s in seps:
        ax.axvline(s + 0.5, color=SER[1], lw=1.5); ax.axhline(s + 0.5, color=SER[1], lw=1.5)
    ax.set_xticks(range(1, nb + 1)); ax.set_yticks(range(1, nb + 1)); ax.grid(False)
    if ticks: ax.set_xticklabels(ticks, fontsize=7, rotation=45, ha='right'); ax.set_yticklabels(ticks, fontsize=7)
    ax.set_xlabel('true bin'); ax.set_ylabel('reco bin'); ax.set_title(title, fontsize=10, loc='left')

def edge_labels(e):
    return [f'{e[i]:.2g}-{e[i + 1]:.2g}' for i in range(len(e) - 1)]

def fig_response_proton():
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.6))
    for r, mode in enumerate(('fhc', 'rhc')):
        d = load(mode, '1p'); sig = (d['MC_Signal'] > 0) & (d['Selected'] > 0); w = d['w'][sig]
        specs = [('theta_p [rad], 4 bins (built)', np.arccos(np.clip(d['proton_costh_true'], -1, 1)), np.arccos(np.clip(d['proton_costh_reco'], -1, 1)), [0, 0.597, 0.942, 1.287, PI]),
                 ('cos theta_p, 4 bins', d['proton_costh_true'], d['proton_costh_reco'], [-1, 0.275, 0.575, 0.825, 1]),
                 ('theta_pi-p [rad], 3 bins (built)', d['thpipr_true'], d['thpipr_reco'], [0, 1.193, 2.010, PI])]
        for c, (t, xt, xr, e) in enumerate(specs):
            R = response(binidx(xt[sig], e, False), binidx(xr[sig], e, False), w, len(e) - 1)
            heat(axs[r, c], R, f'{mode.upper()}  {t}\nmin diagonal {np.diag(R).min():.2f}', edge_labels(e))
    fig.suptitle('Proton-tagged response matrices (column-normalised; columns sum to <1 where events reconstruct outside the range)',
                 x=0.01, ha='left', fontsize=12)
    fig.tight_layout(); save(fig, 'response_proton_angles')

def fig_response_2d():
    for mode in ('fhc', 'rhc'):
        di, dp = load(mode, 'incl'), load(mode, '1p')
        def get(d, var, reco):
            s = 'reco' if reco else 'true'
            if var == 'costhmu': return d[f'candidate_muon_costh_{s}']
            if var == 'costhpi': return d[f'candidate_pion_costh_{s}']
            if var == 'pmu': return d[f'candidate_muon_mom_{s}']
            if var == 'thmupi': return d['mu_pi_opening_angle' if reco else 'true_mu_pi_opening_angle']
            if var == 'ppi': return ppi_reco(d['candidate_pion_mom_reco']) if reco else d['candidate_pion_mom_true']
            if var == 'thetap': return np.arccos(np.clip(d[f'proton_costh_{s}'], -1, 1))
            if var == 'thpipr': return d[f'thpipr_{s}']
            if var == 'dpt': return d[f'deltaPt_{s}']
            if var == 'Wpipr': return d[f'W_pipr_{s}']
        specs = [(di, 'costhpi', [-1, 0.35, 1], False, 'costhmu', [-1, 0.65, 0.85, 1], False, 'cos theta_pi | cos theta_mu'),
                 (di, 'costhmu', [-1, 0.8, 1], False, 'pmu', [0.15, 0.55, 3.0], True, 'cos theta_mu | p_mu'),
                 (di, 'costhpi', [-1, 0.35, 1], False, 'thmupi', [0, 0.85, 1.5, 2.6], False, 'cos theta_pi | theta_mu-pi'),
                 (di, 'costhpi', [-1, 0.35, 1], False, 'ppi', [0.175, 0.205, 1.0], True, 'cos theta_pi | p_pi'),
                 (dp, 'thetap', [0, 0.8, PI], False, 'dpt', [0, 0.3, 2.5], True, 'theta_p | delta p_T'),
                 (dp, 'thpipr', [0, 1.2, PI], False, 'dpt', [0, 0.3, 2.5], True, 'theta_pi-p | delta p_T'),
                 (dp, 'thetap', [0, 0.8, PI], False, 'Wpipr', [1.08, 1.19, 2.9], True, 'theta_p | W_pi-p'),
                 (dp, 'thpipr', [0, 1.2, PI], False, 'Wpipr', [1.08, 1.19, 2.9], True, 'theta_pi-p | W_pi-p')]
        fig, axs = plt.subplots(2, 4, figsize=(18, 9.4))
        for ax, (d, xv, xe, xo, yv, ye, yo, title) in zip(axs.flat, specs):
            sig = (d['MC_Signal'] > 0) & (d['Selected'] > 0); w = d['w'][sig]; nx = len(xe) - 1
            def flat(reco):
                i = binidx(get(d, xv, reco)[sig], xe, xo); j = binidx(get(d, yv, reco)[sig], ye, yo)
                return np.where((i >= 0) & (j >= 0), j * nx + i, -1)
            nb = nx * (len(ye) - 1); R = response(flat(False), flat(True), w, nb)
            labs = [f'Y{j + 1}X{i + 1}' for j in range(len(ye) - 1) for i in range(nx)]
            heat(ax, R, f'{title}\nmin diagonal {np.diag(R).min():.2f}', labs, seps=[nx * k for k in range(1, len(ye) - 1)])
        fig.suptitle(f'2D response matrices, {mode.upper()} (flattened Y slice by Y slice; orange lines separate Y slices)',
                     x=0.01, ha='left', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97)); save(fig, f'response_2d_{mode}')

def fig_resolution():
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    d = {k: np.concatenate([load(m, '1p')[k] for m in ('fhc', 'rhc')]) for k in load('fhc', '1p')}
    sel = (d['Selected'] > 0) & (d['MC_Signal'] > 0) & (d['prmatch'] >= 0)
    cats = [(2, 'leading true proton'), (1, 'sub-leading proton'), (0, 'not a proton')]
    tot = d['w'][sel].sum()
    thp_t = np.arccos(np.clip(d['proton_costh_true'], -1, 1)); thp_r = np.arccos(np.clip(d['proton_costh_reco'], -1, 1))
    for ax, (lab, dv, rng) in zip(axs[:2], [('theta_p reco - true [rad]', thp_r - thp_t, (-1.5, 1.5)),
                                             ('theta_pi-p reco - true [rad]', d['thpipr_reco'] - d['thpipr_true'], (-2.5, 2.5))]):
        bins = np.linspace(*rng, 61)
        for k, (c, name) in enumerate(cats):
            m = sel & (d['prmatch'] == c); frac = d['w'][m].sum() / tot
            ax.hist(dv[m], bins=bins, weights=d['w'][m], histtype='step', lw=2, color=SER[k],
                    ls=['-', '--', ':'][k], label=f'{name} ({100 * frac:.0f}%)')
        ax.set_xlabel(lab); ax.set_ylabel('selected signal (data POT)'); ax.legend(fontsize=8)
    # match fraction vs true proton momentum
    ax = axs[2]
    ax.set_title('share of selected signal by candidate truth vs true theta_p', fontsize=10, loc='left')
    tb = np.linspace(0, PI, 13); ctr = 0.5 * (tb[1:] + tb[:-1]); bottom = np.zeros(len(ctr))
    idx = np.digitize(thp_t, tb) - 1
    denom = np.array([d['w'][sel & (idx == i)].sum() for i in range(len(ctr))])
    for k, (c, name) in enumerate(cats):
        f = np.array([d['w'][sel & (idx == i) & (d['prmatch'] == c)].sum() for i in range(len(ctr))])
        f = np.divide(f, denom, out=np.zeros_like(f), where=denom > 0)
        ax.bar(ctr, f, width=tb[1] - tb[0], bottom=bottom, color=SER[k], edgecolor=SURF, linewidth=2, label=name)
        bottom += f
    ax.set_xlabel('true theta_p [rad]'); ax.set_ylabel('fraction'); ax.set_ylim(0, 1.05); ax.legend(fontsize=8, loc='lower right')
    fig.suptitle('Proton-tagged selected SIGNAL events, FHC+RHC: angular resolution split by what the reconstructed proton candidate is (backtracked)',
                 x=0.01, ha='left', fontsize=12)
    fig.tight_layout(); save(fig, 'resolution_proton_angles')

# ------------------------------------------------------------------------------------------------
# 4. generator predictions
GENS = [('genie', 'GENIE'), ('gibuu', 'GiBUU'), ('neut', 'NEUT'), ('nuwro', 'NuWro')]

def fte(g, fam, mode, obs):
    import uproot
    h = uproot.open(os.path.join(GEN, f'{g}_{fam}_{mode}_fte.root'))[obs + '_fte']
    v, err = h.values(), h.errors()
    return v, err

def fig_gen_proton():
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.6), sharey='row')
    for r, (obs, e, lab) in enumerate([('thetap', [0, 0.597, 0.942, 1.287, PI], 'theta_p [rad]'),
                                         ('thpipr', [0, 1.193, 2.010, PI], 'theta_pi-p [rad]')]):
        e = np.array(e); wdt = np.diff(e)
        for c, mode in enumerate(('fhc', 'rhc', 'comb')):
            ax = axs[r, c]
            for k, (g, G) in enumerate(GENS):
                v, err = fte(g, '1p_ext', mode, obs)
                y = v / wdt
                ax.stairs(y, e, color=SER[k], lw=2, ls=['-', '--', '-.', ':'][k], label=f'{G} (sigma {v.sum():.2f})')
            ax.set_title(f'{mode.upper()}', fontsize=10, loc='left'); ax.set_xlabel(lab)
            if c == 0: ax.set_ylabel('dsigma/dx [1e-38 cm2 / Ar / rad]')
            ax.legend(fontsize=7)
    fig.suptitle('Generator predictions, proton-tagged signal, analysis binning (sigma in 1e-38 cm2/Ar)', x=0.01, ha='left', fontsize=12)
    fig.tight_layout(); save(fig, 'gen_proton_angles')

def fig_gen_2d():
    specs = [('ext', 'costhpi_costhmu', [-1, 0.35, 1], [-1, 0.65, 0.85, 1], 'cos theta_pi', 'cos theta_mu'),
             ('ext', 'costhmu_pmu', [-1, 0.8, 1], [0.15, 0.55, 3.0], 'cos theta_mu', 'p_mu [GeV/c]'),
             ('ext', 'costhpi_thmupi', [-1, 0.35, 1], [0, 0.85, 1.5, 2.6], 'cos theta_pi', 'theta_mu-pi'),
             ('ext', 'costhpi_ppi', [-1, 0.35, 1], [0.175, 0.205, 1.0], 'cos theta_pi', 'p_pi [GeV/c]'),
             ('1p_ext', 'thetap_dpt', [0, 0.8, PI], [0, 0.3, 2.5], 'theta_p', 'delta p_T [GeV/c]'),
             ('1p_ext', 'thpipr_dpt', [0, 1.2, PI], [0, 0.3, 2.5], 'theta_pi-p', 'delta p_T [GeV/c]'),
             ('1p_ext', 'thetap_Wpipr', [0, 0.8, PI], [1.08, 1.19, 2.9], 'theta_p', 'W_pi-p [GeV/c2]'),
             ('1p_ext', 'thpipr_Wpipr', [0, 1.2, PI], [1.08, 1.19, 2.9], 'theta_pi-p', 'W_pi-p [GeV/c2]')]
    for mode in ('fhc', 'rhc', 'comb'):
        fig, axs = plt.subplots(2, 4, figsize=(18, 8))
        for ax, (fam, obs, xe, ye, xl, yl) in zip(axs.flat, specs):
            nx, ny = len(xe) - 1, len(ye) - 1
            area = np.array([(ye[j + 1] - ye[j]) * (xe[i + 1] - xe[i]) for j in range(ny) for i in range(nx)])
            xs = np.arange(1, nx * ny + 1)
            for k, (g, G) in enumerate(GENS):
                v, err = fte(g, fam, mode, obs)
                ax.errorbar(xs + (k - 1.5) * 0.12, v / area, yerr=err / area, fmt=MARK[k], ms=6, color=SER[k], lw=1.5, label=G)
            for s in range(1, ny): ax.axvline(nx * s + 0.5, color=INK2, lw=1)
            ax.set_xticks(xs)
            ax.set_xticklabels([f'{xe[i]:.3g} to {xe[i + 1]:.3g}' for j in range(ny) for i in range(nx)], fontsize=8)
            ax.set_xlabel(f'{xl} bin, grouped by {yl} slice')
            top = ax.get_ylim()[1]
            for j in range(ny):
                ax.text(nx * j + (nx + 1) / 2, top, f'{ye[j]:.3g} to {ye[j + 1]:.3g}',
                        ha='center', va='bottom', fontsize=8, color=INK2)
            ax.set_ylabel('d2sigma/dXdY [1e-38 cm2/Ar/unit area]'); ax.set_title(f'{xl} | {yl}', fontsize=10, loc='left', pad=16)
        axs.flat[0].legend(fontsize=8)
        fig.suptitle(f'Generator predictions for the 2D pairs, {mode.upper()} (open top Y slices use the nominal width)',
                     x=0.01, ha='left', fontsize=12)
        fig.tight_layout(); save(fig, f'gen_2d_{mode}')

if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['diagonals', '2d_leak', 'gen_proton', 'gen_2d', 'response_proton', 'response_2d', 'resolution']
    for w in which: globals()['fig_' + w]()
