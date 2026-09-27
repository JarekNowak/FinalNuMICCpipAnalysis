"""xsec2d_figs.py -- the two-dimensional cross sections, as slices in the outer variable and as maps.

The extraction writes the flattened measurement (sigma per analysis bin, ordered outer slice by outer
slice) to closure_hists_all_xsec_*.root; dividing by the bin area gives the double-differential cross
section. Two figures per pair and configuration:
  <pair>_<cfg>_slices : one panel per outer bin, d2sigma/dXdY against the inner variable, with the
                        four generators and the fake-data truth;
  <pair>_<cfg>_map    : the same numbers as a map, beside the ratio to the fake-data truth, in the
                        kBird palette the released note uses for its matrix figures.
The data are the per-run Poisson FAKE DATA: the signal region is blind.
2026-09-27 (review 6.2, 6.3): the double-differential axis and colour-bar labels name the two variables
and their units instead of "dX dY ... /unit", bin ranges are written as intervals (with pi and minus
signs), and an open-ended outer bin is labelled "> edge" and drawn with a dashed frame (slices) or a
dashed outline (maps), with a note that its height is the integral over the plotted width. Numbers
unchanged; larger text; review terminology (CV pseudo-data, pseudo-data truth).
   python3 report/tools/xsec2d_figs.py [pair ...]      (default: every pair with a sidecar)
"""
import os, sys
import numpy as np, uproot
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xsec_figs_style import apply, GEN, TRUTH, DATA, BIRD, ink

apply()
plt.rcParams.update({'font.size': 12, 'axes.titlesize': 12, 'axes.labelsize': 12,
                     'xtick.labelsize': 11, 'ytick.labelsize': 11})
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
RB = '/data/uboone/processed/rebuild_2d'
OUT = os.path.join(REPO, 'report', 'figures', 'bkgfit')
sys.path.insert(0, os.path.join(REPO, 'xsec_analyzer', 'scripts'))
from newobs_candidates import CANDS

LAB = {'costhpi': r'$\cos\theta_\pi$', 'costhmu': r'$\cos\theta_\mu$', 'pmu': r'$p_\mu$ [GeV/c]',
       'thmupi': r'$\theta_{\mu\pi}$ [rad]', 'ppi': r'$p_\pi$ [GeV/c]', 'thetap': r'$\theta_p$ [rad]',
       'thpipr': r'$\theta_{\pi p}$ [rad]', 'dpt': r'$\delta p_T$ [GeV/c]', 'Wpipr': r'$W_{\pi p}$ [GeV/c$^2$]'}
UNIT = {'costhpi': '', 'costhmu': '', 'pmu': ' GeV/c', 'thmupi': ' rad', 'ppi': ' GeV/c',
        'thetap': ' rad', 'thpipr': ' rad', 'dpt': ' GeV/c', 'Wpipr': ' GeV/c$^2$'}
SYM = {'costhpi': r'\cos\theta_\pi', 'costhmu': r'\cos\theta_\mu', 'pmu': r'p_\mu', 'thmupi': r'\theta_{\mu\pi}',
       'ppi': r'p_\pi', 'thetap': r'\theta_p', 'thpipr': r'\theta_{\pi p}', 'dpt': r'\delta p_T', 'Wpipr': r'W_{\pi p}'}
CFG = [('FHC5', 'FHC'), ('RHCFULL', 'RHC'), ('COMB', 'FHC+RHC combined')]


def edge(x):
    """A bin edge in mathtext: pi by name, a minus sign rather than a hyphen."""
    return r'\pi' if abs(x - np.pi) < 2e-3 else f'{x:g}'


def d2label(c):
    """d2sigma/dX dY with the two variables named and their units, over two lines."""
    u = [UNIT[v].strip() for v in (c['xvar'], c['yvar']) if UNIT[v].strip()]
    per = '' if not u else ('/' + u[0] if len(u) == 1 and '/' not in u[0] else '/(' + ' '.join(u) + ')')
    return (rf"$\mathrm{{d}}^2\sigma/\mathrm{{d}}{SYM[c['xvar']]}\,\mathrm{{d}}{SYM[c['yvar']]}$" + '\n'
            + rf'[$10^{{-38}}$ cm$^2${per}/Ar]')


def is_top_open(c, ax, e, i):
    return i == len(e) - 2 and (c.get(ax + 'open', False) or not np.isfinite(e[i + 1]))


def is_low_open(c, ax, i):
    return i == 0 and c.get(ax + 'low_open', False)


def slice_title(c, ye, j):
    v = c['yvar']; s = SYM[v]
    if is_top_open(c, 'y', ye, j): return rf'${s} > {edge(ye[j])}${UNIT[v]}'
    if is_low_open(c, 'y', j): return rf'${s} < {edge(ye[j + 1])}${UNIT[v]}'
    return rf'${edge(ye[j])} < {s} < {edge(ye[j + 1])}${UNIT[v]}'


def tick_label(c, ax, e, i):
    if is_top_open(c, ax, e, i): return rf'$> {edge(e[i])}$'
    if is_low_open(c, ax, i): return rf'$< {edge(e[i + 1])}$'
    return rf'$[{edge(e[i])},\ {edge(e[i + 1])}]$'


def open_note(c, ye, j, mark):
    """What the height of an open-ended outer bin means, with the plotted range it is divided by."""
    v = c['yvar']
    return (f'{mark}: open-ended ${SYM[v]}$ bin, height = integral / plotted width '
            rf'(${SYM[v]}$ drawn over ${edge(ye[j])}$–${edge(ye[j + 1])}${UNIT[v]})')


def pairs():
    out = {}
    for c in CANDS:
        if c.get('yvar'): out[c['name']] = c
    return out


def load(pfx, cfg, name, xe, ye):
    p = f'{RB}/closure_hists_all_xsec_{pfx}_{cfg}_{name}.root'
    if not os.path.exists(p): return None
    f = uproot.open(p); nx, ny = len(xe) - 1, len(ye) - 1
    area = np.array([[(xe[i + 1] - xe[i]) * (ye[j + 1] - ye[j]) for i in range(nx)] for j in range(ny)])
    get = lambda k: (f[k].values().reshape(ny, nx) / area) if k in [kk.split(';')[0] for kk in f.keys()] else None
    d = dict(data=get('h_unfolded_nuwro'), truth=get('h_fakedata_truth'),
             err=f['h_unfolded_nuwro'].errors().reshape(ny, nx) / area)
    for k, lab, col, ls in GEN: d[lab] = get(k)
    return d


def fig_slices(name, c, cfg, lab, d):
    xe, ye = np.array(c['xedges']), np.array(c['yedges']); ny = len(ye) - 1
    fig, axs = plt.subplots(1, ny, figsize=(3.4 * ny + 0.8, 3.6), squeeze=False)   # own scale per slice:
    # the outer bins differ by an order of magnitude and a shared axis hides the structure of all but the first
    for j, ax in enumerate(axs[0]):
        ax.stairs(d['truth'][j], xe, color=TRUTH, lw=2.2, label='Pseudo-data truth ($A_C$-smeared)' if j == 0 else None)
        for k, gl, col, ls in GEN:
            if d[gl] is None: continue
            ax.stairs(d[gl][j], xe, color=col, ls=ls, lw=1.8, label=gl if j == 0 else None)
        ctr = 0.5 * (xe[1:] + xe[:-1])
        ax.errorbar(ctr, d['data'][j], yerr=d['err'][j], xerr=np.diff(xe) / 2, fmt='o', color=DATA,
                    ms=4, lw=1.2, capsize=0, label='Unfolded CV pseudo-data' if j == 0 else None)
        ax.set_xlabel(LAB[c['xvar']]); ax.set_xlim(xe[0], xe[-1]); ax.set_ylim(bottom=0)
        ax.set_ylabel(d2label(c) if j == 0 else None)
        ax.set_title(slice_title(c, ye, j), fontsize=12)
        if is_top_open(c, 'y', ye, j):   # open-ended outer bin: dashed frame (explained under the figure)
            for sp in ax.spines.values(): sp.set_linestyle((0, (4, 2.5))); sp.set_edgecolor('0.35')
    fig.suptitle(f"{LAB[c['xvar']].split(' [')[0]} $\\times$ {LAB[c['yvar']].split(' [')[0]}, {lab}"
                 "  (CV pseudo-data; signal region blind)", fontsize=12, y=1.02)
    fig.tight_layout()
    # one legend for the figure, under the panels, so that it never covers a curve
    h, l = axs[0][0].get_legend_handles_labels()
    order = [l.index('Unfolded CV pseudo-data')] + [i for i in range(len(l)) if l[i] != 'Unfolded CV pseudo-data']
    ncol = 6 if ny >= 3 else 3
    fig.legend([h[i] for i in order], [l[i] for i in order], loc='upper center', bbox_to_anchor=(0.5, 0.0),
               ncol=ncol, fontsize=10.5, frameon=False, columnspacing=1.6, handlelength=2.2)
    if is_top_open(c, 'y', ye, ny - 1):
        rows = -(-len(order) // ncol)
        fig.text(0.5, -0.03 - 0.075 * rows, open_note(c, ye, ny - 1, 'Dashed frame'), ha='center', va='top',
                 fontsize=10.5, color='0.25')
    for ext in ('pdf', 'png'): fig.savefig(f'{OUT}/xsec2d_{name}_{cfg}_slices.{ext}', bbox_inches='tight')
    plt.close(fig)


def fig_map(name, c, cfg, lab, d):
    """Equal-size cells (index space) labelled by their ranges: the physical edges differ by an order of
    magnitude, which leaves the interesting cells too small to read."""
    xe, ye = np.array(c['xedges'], float), np.array(c['yedges'], float)
    nx, ny = len(xe) - 1, len(ye) - 1
    fig, axs = plt.subplots(1, 2, figsize=(5.2 + 2.0 * nx, 1.6 + 1.5 * ny))
    r = np.divide(d['data'], d['truth'], out=np.full_like(d['data'], np.nan), where=d['truth'] > 0)
    lim = max(0.35, float(np.nanmax(np.abs(r - 1))))
    for ax, (Z, cmap, vmin, vmax, title, cbl, fmt) in zip(axs, [
            (d['data'], BIRD, 0, None, 'Unfolded CV pseudo-data', d2label(c), '.3f'),
            (r, BIRD, 1 - lim, 1 + lim, 'Closure: unfolded / pseudo-data truth', 'ratio', '.2f')]):
        m = ax.pcolormesh(np.arange(nx + 1), np.arange(ny + 1), Z, cmap=cmap, vmin=vmin, vmax=vmax)
        for sp in ax.spines.values(): sp.set_zorder(3)   # contiguous cells, as COLZ draws them
        fig.colorbar(m, ax=ax, label=cbl, fraction=0.046, pad=0.03)
        for j in range(ny):
            for i in range(nx):
                if not np.isfinite(Z[j, i]): continue
                e = d['err'][j, i] if fmt == '.3f' else d['err'][j, i] / d['truth'][j, i]
                ax.text(i + 0.5, j + 0.5, f'{Z[j, i]:{fmt}}\n$\pm${e:{fmt}}', ha='center', va='center',
                        fontsize=11, color=ink(m.cmap(m.norm(Z[j, i]))))
        ax.set_xticks(np.arange(nx) + 0.5); ax.set_xticklabels([tick_label(c, 'x', xe, i) for i in range(nx)])
        ax.set_yticks(np.arange(ny) + 0.5); ax.set_yticklabels([tick_label(c, 'y', ye, j) for j in range(ny)])
        ax.set_xlabel(LAB[c['xvar']]); ax.set_ylabel(LAB[c['yvar']]); ax.set_title(title, fontsize=12)
        ax.tick_params(length=0)
        if is_top_open(c, 'y', ye, ny - 1):   # open-ended outer bin: dashed outline of its row
            ax.add_patch(plt.Rectangle((0.03, ny - 1 + 0.03), nx - 0.06, 0.94, fill=False, ls=(0, (4, 2.5)),
                                       lw=1.6, ec='white', zorder=4))
    fig.suptitle(f"{LAB[c['xvar']].split(' [')[0]} $\\times$ {LAB[c['yvar']].split(' [')[0]}, {lab}"
                 "  (CV pseudo-data; signal region blind)", fontsize=12, y=1.03)
    fig.tight_layout()
    if is_top_open(c, 'y', ye, ny - 1):
        fig.text(0.5, -0.02, open_note(c, ye, ny - 1, 'Dashed outline'), ha='center', va='top',
                 fontsize=10.5, color='0.25')
    for ext in ('pdf', 'png'): fig.savefig(f'{OUT}/xsec2d_{name}_{cfg}_map.{ext}', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    want = sys.argv[1:]
    P = pairs(); n = 0
    for name, c in P.items():
        if want and name not in want: continue
        for cfg, lab in CFG:
            ye = np.array(c['yedges'], float)
            d = load(c['pfx'], cfg, name, np.array(c['xedges'], float), np.where(np.isinf(ye), 1e9, ye))
            if d is None or d['data'] is None: continue
            fig_slices(name, c, cfg, lab, d); fig_map(name, c, cfg, lab, d); n += 1
            print(f'wrote xsec2d_{name}_{cfg}_{{slices,map}}')
    print(f'{n} configurations plotted')
