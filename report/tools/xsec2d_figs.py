"""xsec2d_figs.py -- the two-dimensional cross sections, as slices in the outer variable and as maps.

The extraction writes the flattened measurement (sigma per analysis bin, ordered outer slice by outer
slice) to closure_hists_all_xsec_*.root; dividing by the bin area gives the double-differential cross
section. Two figures per pair and configuration:
  <pair>_<cfg>_slices : one panel per outer bin, d2sigma/dXdY against the inner variable, with the
                        four generators and the fake-data truth;
  <pair>_<cfg>_map    : the same numbers as a map, beside the ratio to the fake-data truth, in the
                        kBird palette the released note uses for its matrix figures.
The data are the per-run Poisson FAKE DATA: the signal region is blind.
   python3 report/tools/xsec2d_figs.py [pair ...]      (default: every pair with a sidecar)
"""
import os, sys
import numpy as np, uproot
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xsec_figs_style import apply, GEN, TRUTH, DATA, BIRD, ink

apply()
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
CFG = [('FHC5', 'FHC'), ('RHCFULL', 'RHC'), ('COMB', 'FHC+RHC combined')]


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
    fig, axs = plt.subplots(1, ny, figsize=(3.6 * ny + 0.6, 3.6), squeeze=False)   # own scale per slice:
    # the outer bins differ by an order of magnitude and a shared axis hides the structure of all but the first
    for j, ax in enumerate(axs[0]):
        ax.stairs(d['truth'][j], xe, color=TRUTH, lw=2.2, label='fake-data truth' if j == 0 else None)
        for k, gl, col, ls in GEN:
            if d[gl] is None: continue
            ax.stairs(d[gl][j], xe, color=col, ls=ls, lw=1.8, label=gl if j == 0 else None)
        ctr = 0.5 * (xe[1:] + xe[:-1])
        ax.errorbar(ctr, d['data'][j], yerr=d['err'][j], xerr=np.diff(xe) / 2, fmt='o', color=DATA,
                    ms=4, lw=1.2, capsize=0, label='unfolded fake data' if j == 0 else None)
        ax.set_xlabel(LAB[c['xvar']]); ax.set_xlim(xe[0], xe[-1]); ax.set_ylim(bottom=0)
        ax.set_ylabel(r'$\mathrm{d}^2\sigma/\mathrm{d}X\,\mathrm{d}Y$  [$10^{-38}$cm$^2$/Ar/unit]' if j == 0 else None)
        hi = ye[j + 1] if np.isfinite(ye[j + 1]) else np.inf
        ax.set_title(f"{LAB[c['yvar']].split(' [')[0]} {ye[j]:g}--{hi:g}{UNIT[c['yvar']]}", fontsize=9)
        if j == 0: ax.legend(fontsize=7.5, loc='best')
    fig.suptitle(f"{LAB[c['xvar']].split(' [')[0]} $\\times$ {LAB[c['yvar']].split(' [')[0]}, {lab}"
                 "  (fake data; signal region blind)", fontsize=10, y=1.02)
    fig.tight_layout()
    for ext in ('pdf', 'png'): fig.savefig(f'{OUT}/xsec2d_{name}_{cfg}_slices.{ext}', bbox_inches='tight')
    plt.close(fig)


def fig_map(name, c, cfg, lab, d):
    """Equal-size cells (index space) labelled by their ranges: the physical edges differ by an order of
    magnitude, which leaves the interesting cells too small to read."""
    xe, ye = np.array(c['xedges'], float), np.array(c['yedges'], float)
    nx, ny = len(xe) - 1, len(ye) - 1
    rng = lambda e, i: f'{e[i]:g}--' + (r'$\infty$' if not np.isfinite(e[i + 1]) else f'{e[i + 1]:g}')
    fig, axs = plt.subplots(1, 2, figsize=(5.2 + 2.0 * nx, 1.6 + 1.5 * ny))
    r = np.divide(d['data'], d['truth'], out=np.full_like(d['data'], np.nan), where=d['truth'] > 0)
    lim = max(0.35, float(np.nanmax(np.abs(r - 1))))
    for ax, (Z, cmap, vmin, vmax, title, cbl, fmt) in zip(axs, [
            (d['data'], BIRD, 0, None, 'unfolded fake data', r'$\mathrm{d}^2\sigma/\mathrm{d}X\,\mathrm{d}Y$ [$10^{-38}$cm$^2$/Ar/unit]', '.3f'),
            (r, BIRD, 1 - lim, 1 + lim, 'closure: unfolded / fake-data truth', 'ratio', '.2f')]):
        m = ax.pcolormesh(np.arange(nx + 1), np.arange(ny + 1), Z, cmap=cmap, vmin=vmin, vmax=vmax)
        for sp in ax.spines.values(): sp.set_zorder(3)   # contiguous cells, as COLZ draws them
        fig.colorbar(m, ax=ax, label=cbl, fraction=0.046, pad=0.03)
        for j in range(ny):
            for i in range(nx):
                if not np.isfinite(Z[j, i]): continue
                e = d['err'][j, i] if fmt == '.3f' else d['err'][j, i] / d['truth'][j, i]
                ax.text(i + 0.5, j + 0.5, f'{Z[j, i]:{fmt}}\n$\pm${e:{fmt}}', ha='center', va='center',
                        fontsize=9, color=ink(m.cmap(m.norm(Z[j, i]))))
        ax.set_xticks(np.arange(nx) + 0.5); ax.set_xticklabels([rng(xe, i) for i in range(nx)])
        ax.set_yticks(np.arange(ny) + 0.5); ax.set_yticklabels([rng(ye, j) for j in range(ny)])
        ax.set_xlabel(LAB[c['xvar']]); ax.set_ylabel(LAB[c['yvar']]); ax.set_title(title, fontsize=10)
        ax.tick_params(length=0)
    fig.suptitle(f"{LAB[c['xvar']].split(' [')[0]} $\\times$ {LAB[c['yvar']].split(' [')[0]}, {lab}"
                 "  (fake data; signal region blind)", fontsize=10, y=1.03)
    fig.tight_layout()
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
