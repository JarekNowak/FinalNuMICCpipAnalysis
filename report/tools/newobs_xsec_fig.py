# newobs_xsec_fig.py -- extracted cross sections of the observables added 2026-09-16 (theta_p, theta_pi-p),
# from the closure sidecars of the slurm_extract runs. The "data" are the per-run Poisson FAKE DATA: the signal
# region is blind. Curves are A_C-smeared, as the extraction plots them.
#   python3 report/tools/newobs_xsec_fig.py  -> report/figures/bkgfit/newobs_xsec.{pdf,png}
import os, sys
import numpy as np, uproot
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xsec_figs_style import apply, GEN, TRUTH, DATA   # the style of the released figures

apply()
P = '/data/uboone/processed/closure_hists_xsec_ccpi1p_%s_%s.root'
OBS = [('thetap', r'$\theta_p$ [rad]'), ('thpipr', r'$\theta_{\pi p}$ [rad]')]
CFG = [('FHC5', 'FHC'), ('RHCFULL', 'RHC'), ('COMB', 'FHC+RHC combined')]

fig, axs = plt.subplots(2, 3, figsize=(15, 8))
for r, (obs, xl) in enumerate(OBS):
    for c, (cfg, lab) in enumerate(CFG):
        ax = axs[r, c]; p = P % (cfg, obs)
        if not os.path.exists(p): ax.text(0.5, 0.5, 'pending', ha='center', transform=ax.transAxes); continue
        f = uproot.open(p); h = f['h_unfolded_nuwro']; e = h.axis().edges(); ctr = 0.5 * (e[1:] + e[:-1])
        v, err = h.values(), h.errors()
        for key, gl, col, ls in GEN:
            if key not in [kk.split(';')[0] for kk in f.keys()]: continue
            g = f[key].values(); ax.stairs(g, e, color=col, ls=ls, lw=1.8,
                                           label=f'{gl} ({(g * np.diff(e)).sum():.2f})')
        t = f['h_fakedata_truth'].values(); ax.stairs(t, e, color=TRUTH, lw=2.5, label=f'fake-data truth ({(t * np.diff(e)).sum():.2f})')
        ax.errorbar(ctr, v, yerr=err, xerr=np.diff(e) / 2, fmt='o', color=DATA, ms=5, capsize=0,
                    label=f'unfolded fake data ({(v * np.diff(e)).sum():.2f})')
        ax.set_xlabel(xl); ax.set_ylim(bottom=0); ax.set_xlim(e[0], e[-1])
        if c == 0: ax.set_ylabel(r'd$\sigma$/d$\theta$ [$10^{-38}$ cm$^2$ / Ar / rad]')
        ax.set_title(f'{lab}', loc='left', fontsize=10); ax.legend(fontsize=7)
fig.suptitle('Proton-tagged CC1$\\mu$1$\\pi^{\\pm}$: new observables, FAKE DATA (signal region blind); '
             'curves and integrated $\\sigma$ [$10^{-38}$ cm$^2$/Ar] are $A_C$-smeared', x=0.01, ha='left')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(f'report/figures/bkgfit/newobs_xsec.{ext}', bbox_inches='tight')
print('wrote report/figures/bkgfit/newobs_xsec.{pdf,png}')
