"""update_talk_figures.py -- figures of the update talk report/slides/update_talk_2026-10.tex.

  cr_mcscal_chi2.pdf   shape chi2 of the muon momentum in the four control regions, FHC and RHC, beam-on data, frozen
                       protocol, with the version-1.7 simulation and with the MCS-calibrated simulation (version 1.9),
                       read from report/results/cr_protocol_mcscal.md (report/tools/cr_protocol_compare.py)
  totals_update.pdf    one-bin totals at the collaboration meeting (report/slides/collab_talk_imperial_2026.tex, Draft 1.5)
                       and in version 1.9 (report/data_release/total_xsec.tsv), inclusive and proton-tagged

    python3 report/tools/update_talk_figures.py  ->  report/figures/update_talk/
"""
import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chi2dist

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(R, 'figures', 'update_talk')
plt.rcParams.update({'font.size': 11})
# the one-bin totals shown at the meeting (collab_talk_imperial_2026.tex, table of the totals)
MEETING = {'incl': {'fhc5': (1.03, 0.40), 'rhcfull': (0.93, 0.35), 'comb': (0.97, 0.37)},
           '1p': {'fhc5': (0.62, 0.26), 'rhcfull': (0.53, 0.25), 'comb': (0.57, 0.24)}}


def cr_chi2():
    txt = open(os.path.join(R, 'results', 'cr_protocol_mcscal.md')).read()
    beam = txt[txt.index('## beam-on data'):txt.index('## pseudo-data')]
    rows = re.findall(r'^\| (fhc|rhc) \| (\w+)_pmu \| ([\d.]+)/(\d+) \| [^|]* \| ([\d.]+)/(\d+) \|', beam, re.M)
    names = {'cc0pi': 'CC0$\\pi$', 'pi0': '$\\pi^0$', 'multipi': 'multi-$\\pi$', 'cosmic': 'cosmic'}
    lab = [m.upper() + ' ' + names.get(r, r) for m, r, *_ in rows]
    old = np.array([float(x[2]) for x in rows]); new = np.array([float(x[4]) for x in rows]); ndf = int(rows[0][3])
    fig, ax = plt.subplots(figsize=(8.0, 3.6))
    x = np.arange(len(rows))
    ax.bar(x - 0.2, old, 0.4, color='#c8553d', label='version 1.7 simulation')
    ax.bar(x + 0.2, new, 0.4, color='#1f5fa8', label='MCS-calibrated simulation (v1.9)')
    ax.axhline(chi2dist.ppf(0.95, ndf), color='k', ls='--', lw=1, label=f'95% quantile, {ndf} dof')
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=30, ha='right', fontsize=10)
    ax.set_ylabel(f'shape $\\chi^2$ of $p_\\mu$ ({ndf} dof)')
    ax.set_yscale('log'); ax.set_ylim(0.3, 1500)
    ax.legend(fontsize=9.5, frameon=False, loc='upper left', ncol=3)
    s = re.search(r'^\| sum \| p_mu \| ([\d.]+)/(\d+) \|\s*\| ([\d.]+)/(\d+)', beam, re.M)
    ax.text(0.99, 0.74, f'$p_\\mu$ sum: {s.group(1)}/{s.group(2)} $\\to$ {s.group(3)}/{s.group(4)}', transform=ax.transAxes,
            ha='right', fontsize=10)
    ax.text(0.0, 1.02, 'MicroBooNE NuMI data, control regions, internal', transform=ax.transAxes, fontsize=8.5, color='0.45')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'cr_mcscal_chi2.pdf')); plt.close(fig)


def totals():
    v19 = {}
    for l in open(os.path.join(R, 'data_release', 'total_xsec.tsv')):
        if l.startswith('#') or l.startswith('family'): continue
        f = l.split('\t'); v19.setdefault(f[0], {})[f[1]] = (float(f[2]), float(f[3]))
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.4), sharey=False)
    for ax, (fam, title) in zip(axs, (('incl', 'inclusive'), ('1p', 'proton-tagged'))):
        for k, (cfg, nm) in enumerate((('fhc5', 'FHC'), ('rhcfull', 'RHC'), ('comb', 'combined'))):
            for off, src, col, lab in ((-0.12, MEETING, '#c8553d', 'collaboration meeting (Draft 1.5)'),
                                       (0.12, v19, '#1f5fa8', 'version 1.9')):
                v, e = src[fam][cfg]
                ax.errorbar(k + off, v, yerr=e, fmt='o', color=col, capsize=3, ms=5, label=lab if k == 0 else None)
                ax.text(k + off, v + e + 0.02, f'{100 * e / v:.0f}%', fontsize=8.5, color=col, ha='center', va='bottom')
        ax.set_xticks(range(3)); ax.set_xticklabels(['FHC', 'RHC', 'combined']); ax.set_title(title, fontsize=11)
        ax.set_xlim(-0.5, 2.5); ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, 1.7 if fam == 'incl' else 1.0)
    axs[0].set_ylabel(r'one-bin total [$10^{-38}$ cm$^2$/Ar]')
    axs[0].legend(fontsize=8.5, frameon=False, loc='lower left')
    fig.text(0.01, 0.97, 'MicroBooNE NuMI simulation, pseudo-data (blind), internal', fontsize=8.5, color='0.45', va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(OUT, 'totals_update.pdf')); plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    cr_chi2(); totals()
    print('written', OUT)
