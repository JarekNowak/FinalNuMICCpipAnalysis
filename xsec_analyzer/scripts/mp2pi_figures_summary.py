"""mp2pi_figures_summary.py -- summary figures of the two-pion one-bin total for the slides and the note.

unc_breakdown.pdf: the uncertainty terms of the framework extraction with the MCS-calibrated classifier inputs
(SYSTDUMP of /data/uboone/processed/mp2pi_mcscal/unfold/<cfg>[_constr]_dvfrac_fluxw_mcscal.log), FHC, RHC and
combined, without and with the sideband constraint.

unc_evolution.pdf: the expected total of the combined one-bin total along the analysis (Tables wp, constraint,
framework and mcscal of the note).

    python3 scripts/mp2pi_figures_summary.py  ->  report/figures/multipion/unc_{breakdown,evolution}.pdf
"""
import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'multipion'))
U = '/data/uboone/processed/mp2pi_mcscal/unfold/'
CFGS = [('fhc5', 'FHC'), ('rhcfull', 'RHC'), ('comb', 'Combined')]
TERMS = [('DataStats', 'Data stat.'), ('flux', 'Flux'), ('xsec_total', 'Cross section'), ('reint', 'Reinteraction'),
         ('detVar_total', 'Detector'), ('SimulationStats', 'Sim. stat.')]
COL = ['#4c72b0', '#dd8452', '#55a868', '#8172b3', '#c44e52', '#937860']
plt.rcParams.update({'font.size': 12, 'axes.titlesize': 13})


def dump(cfg, mode):
    t = open(f'{U}{cfg}{mode}_dvfrac_fluxw_mcscal.log').read()
    return {k: float(x) for k, x in re.findall(r'^\[SYSTDUMP\] (\S+) ([-\d.e+]+)$', t, re.M)}


def label(ax, txt):
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=9, color='0.45', va='top')


def breakdown():
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.5), sharey=True)
    for ax, (mode, title) in zip(axs, (('', 'No constraint'), ('_constr', 'Sideband constraint'))):
        D = [dump(c, mode) for c, _ in CFGS]
        x = np.arange(len(CFGS)); w = 0.12
        for j, ((k, n), c) in enumerate(zip(TERMS, COL)):
            ax.bar(x + (j - 2.5) * w, [d[k] for d in D], w, color=c, label=n)
        for i, d in enumerate(D):
            ax.plot([i - 3.2 * w, i + 3.2 * w], [d['total']] * 2, color='k', lw=2)
            ax.text(i, d['total'] + 1.2, f"total {d['total']:.1f}%", ha='center', fontsize=11)
        ax.set_xticks(x); ax.set_xticklabels([n for _, n in CFGS]); ax.set_title(title)
        ax.set_ylim(0, 60); ax.grid(axis='y', alpha=0.3)
        label(ax, 'MicroBooNE simulation, internal')
    axs[0].set_ylabel('Uncertainty of the one-bin total [%]')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=10.5, ncol=6, loc='lower center', frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f'{OUT}/unc_breakdown.pdf'); plt.close(fig)


def evolution():
    steps = [('Cut-based\nfirst draft', 89), ('+ particle\nclassifier', 61), ('+ event\nclassifier', 47),
             ('LLR-and-\ngeometry PID', 46), ('+ sideband\nconstraint', 16.8), ('framework', 15.7),
             ('MCS-calibrated\ninputs', 15.2)]
    fig, ax = plt.subplots(figsize=(9.5, 3.9))
    x = np.arange(len(steps)); v = [s[1] for s in steps]
    c = ['#9e9e9e'] * 4 + ['#c44e52'] * 3
    ax.bar(x, v, 0.62, color=c)
    for i, y in enumerate(v): ax.text(i, y + 1.5, f'{y:g}%', ha='center', fontsize=12)
    ax.set_xticks(x); ax.set_xticklabels([s[0] for s in steps], fontsize=10.5)
    ax.set_ylabel('Total uncertainty [%]'); ax.set_ylim(0, 105); ax.grid(axis='y', alpha=0.3)
    ax.axvline(3.5, color='k', ls=':', lw=1)
    ax.text(1.5, 86, 'study, no constraint', ha='center', fontsize=10.5, color='0.3')
    ax.text(5.0, 86, 'with the sideband constraint', ha='center', fontsize=10.5, color='#c44e52')
    label(ax, 'MicroBooNE simulation, internal; combined FHC+RHC, one-bin total')
    fig.tight_layout(); fig.savefig(f'{OUT}/unc_evolution.pdf'); plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    breakdown(); evolution()
    print('written', OUT)
