"""binreview_figs.py -- figures of the binning review (0.68 versus 0.50 migration-diagonal criteria):
the degrees of freedom for signal (trace of A_C, the sum of the Wiener filter factors) against the
number of bins, for every binning evaluated in logs/binreview/eval.tsv, in FHC, RHC and combined.
Marker fill gives the criterion the binning meets (worst true-bin diagonal over FHC and RHC); the
released binning is drawn larger. Writes figures/binreview_dfs_incl.pdf and binreview_dfs_1p.pdf.
    python3 report/tools/binreview_figs.py
"""
import os, sys, csv, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import xsec_figs_style as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S.apply()
REP = os.path.dirname(HERE)
EVAL = os.path.join(os.path.dirname(REP), 'logs', 'binreview', 'eval.tsv')
CFG = [('fhc5', 'FHC', '#0072B2', '-'), ('rhcfull', 'RHC', '#D55E00', '-'), ('comb', 'combined', '#666666', (0, (4, 2)))]
LAB = {'pmu': r'$p_\mu$', 'ppi': r'$p_\pi$', 'costhmu': r'$\cos\theta_\mu$', 'thetamu': r'$\theta_\mu$',
       'costhpi': r'$\cos\theta_\pi$', 'thmupi': r'$\theta_{\mu\pi}$', 'Whad': r'$W_{\rm had}$', 'Wpipr': r'$W_{\pi p}$',
       'dpt': r'$\delta p_T$', 'dalphat': r'$\delta\alpha_T$', 'dphit': r'$\delta\phi_T$', 'pn': r'$p_n$',
       'thetap': r'$\theta_p$', 'thpipr': r'$\theta_{\pi p}$', 'pp': r'$p_p$'}
PANELS = {'incl': ['pmu', 'ppi', 'costhmu', 'thetamu', 'costhpi', 'thmupi'],
          '1p': ['Whad', 'Wpipr', 'dpt', 'dalphat', 'dphit', 'pn', 'thetap', 'thpipr', 'pp']}

rows = list(csv.DictReader(open(EVAL), delimiter='\t'))
pts = collections.defaultdict(lambda: collections.defaultdict(dict))   # (fam, obs) -> binning -> cfg -> row
for r in rows: pts[(r['family'], r['observable'])][r['binning']][r['config']] = r

def status(d):
    w = min(float(d[c]['min_diag']) for c in ('fhc5', 'rhcfull') if c in d)
    return 2 if w >= 0.68 else (1 if w >= 0.50 else 0)

def panel(ax, fam, obs):
    bb = pts[(fam, obs)]
    kmax, kmin, ymax = 1, 99, 0.0
    for cfg, name, col, ls in CFG:
        # one line through the maximin candidates (the edge-optimised sequence); released and the
        # theta-grid partitions of cos theta_mu are separate points
        seq = sorted(((int(d[cfg]['bins']), float(d[cfg]['dfs']), status(d)) for b, d in bb.items()
                      if cfg in d and b != 'released' and 'T' not in b.replace(obs, '')), key=lambda t: t[0])
        # at a bin count with two candidates (W_had K3 at a lower event floor) keep the one listed last
        seq = list({k: (k, y, s) for k, y, s in seq}.values())
        if seq:
            ax.plot([k for k, _, _ in seq], [y for _, y, _ in seq], color=col, ls=ls, lw=1.2, zorder=2)
            for k, y, s in seq:
                ax.plot(k, y, marker='o' if s else 'x', ms=5.5, mec=col, mfc=col if s == 2 else 'white', mew=1.2, color=col, zorder=3)
                kmax, kmin, ymax = max(kmax, k), min(kmin, k), max(ymax, y)
        for b, d in bb.items():
            if cfg not in d: continue
            k, y, s = int(d[cfg]['bins']), float(d[cfg]['dfs']), status(d)
            kmax, kmin, ymax = max(kmax, k), min(kmin, k), max(ymax, y)
            if b == 'released':      # drawn just left of its bin count so it does not hide a candidate there
                ax.plot(k - 0.22, y, marker='*' if s else 'X', ms=12 if s else 8, mec=col, mfc=col if s == 2 else 'white', mew=1.2, zorder=4)
            elif 'T' in b.replace(obs, ''):
                ax.plot(k + 0.22, y, marker='s', ms=6, mec=col, mfc=col if s == 2 else 'white', mew=1.2, zorder=4)
    ax.set_xlim(kmin - 0.7, kmax + 0.7); ax.set_ylim(0, max(2.0, ymax * 1.18))
    ax.set_xticks(range(kmin, kmax + 1))
    ax.text(0.04, 0.95, LAB[obs], transform=ax.transAxes, va='top', ha='left', fontsize=12)

def figure(fam, ncol, fn):
    obs = PANELS[fam]; nrow = (len(obs) + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.6 * nrow + 0.7), squeeze=False)
    for ax, o in zip(axs.flat, obs): panel(ax, fam, o)
    for ax in list(axs.flat)[len(obs):]: ax.set_visible(False)
    for ax in axs[-1]: ax.set_xlabel('bins')
    for ax in axs[:, 0]: ax.set_ylabel(r'trace $A_C$')
    h = [Line2D([], [], color=c, ls=ls, lw=1.2, label=n) for _, n, c, ls in CFG]
    h += [Line2D([], [], ls='', marker='o', ms=6, mfc='black', mec='black', label=r'passes 0.68'),
          Line2D([], [], ls='', marker='o', ms=6, mfc='white', mec='black', label=r'passes 0.50 only'),
          Line2D([], [], ls='', marker='x', ms=6, mec='black', label='fails 0.50'),
          Line2D([], [], ls='', marker='*', ms=11, mfc='black', mec='black', label='released')]
    if fam == 'incl':
        h.append(Line2D([], [], ls='', marker='s', ms=6, mfc='black', mec='black', label=r'$\theta_\mu$ partition in $\cos\theta_\mu$'))
    fig.legend(handles=h, loc='lower center', ncol=len(h) if len(h) <= 7 else 4, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07 if len(h) <= 7 else 0.1, 1, 1))
    out = os.path.join(REP, 'figures', fn); fig.savefig(out); plt.close(fig); print('wrote', out)

figure('incl', 3, 'binreview_dfs_incl.pdf')
figure('1p', 4, 'binreview_dfs_1p.pdf')
