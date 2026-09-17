"""bkgfit.figures -- figures of the background-fit studies from ../logs/bkgfit/*.json.
   python3 -m bkgfit.figures   -> ../report/figures/bkgfit/*.pdf,png
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import model as M

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.abspath(os.path.join(HERE, '..', '..', 'logs', 'bkgfit'))
FIG = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'bkgfit')); os.makedirs(FIG, exist_ok=True)
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
MARK = ['o', 's', '^', 'D', 'v', 'P']
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'axes.edgecolor': INK2,
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'axes.grid': True,
                     'grid.color': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'font.size': 10,
                     'legend.frameon': False, 'lines.linewidth': 2})


def load(name):
    p = os.path.join(LOG, name + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def save(fig, name):
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(FIG, f'{name}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig); print('wrote', name)


def fig_composition():
    c = load('composition')
    if not c: return
    regions = ['SR', 'CC0pi', 'pi0', 'multipi', 'cosmic']
    fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, beam in zip(axs, ('FHC', 'RHC')):
        row = c[f'full exposure {beam}']; x = np.arange(len(regions)); bottom = np.zeros(len(regions))
        comps = [(n, [row[r]['classes'][n]['frac'] for r in regions]) for n in M.CLASS_NAMES] + [('beam-off', [row[r]['ext'] for r in regions]), ('dirt', [row[r]['dirt'] for r in regions])]
        for k, (n, v) in enumerate(comps):
            ax.bar(x, v, bottom=bottom, color=SER[k % 8], edgecolor=SURF, linewidth=2, label=n); bottom += np.array(v)
        for i, r in enumerate(regions):
            ax.text(i, 1.02, f"nubar {row[r]['nubar_all'] * 100:.0f}%", ha='center', fontsize=8, color=INK2)
        ax.set_xticks(x); ax.set_xticklabels(['signal region', 'CC0pi', 'pi0', 'multi-pi', 'cosmic']); ax.set_ylim(0, 1.1)
        ax.set_title(f'{beam} (full exposure, CV prediction)', loc='left', fontsize=10); ax.set_ylabel('fraction of prediction')
    axs[1].legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.suptitle('Composition by truth class, with the antineutrino share of the neutrino MC above each bar', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'composition')


def fig_sensitivity():
    s = load('sensitivity')
    if not s: return
    F = s['fisher']
    pars = ['CC0pi', 'CC0pi+pi0', 'CC0pi nu/nubar', 'CC0pi+pi0 nu/nubar', 'all classes']
    fig, axs = plt.subplots(1, len(pars), figsize=(19, 4), sharey=True)
    for ax, pn in zip(axs, pars):
        k0 = f'full exposure | FHC+RHC | {pn}'; names = F[k0]['names']; x = np.arange(len(names))
        for j, beams in enumerate(('FHC', 'RHC', 'FHC+RHC')):
            sig = F[f'full exposure | {beams} | {pn}']['sigma']
            ax.bar(x + (j - 1) * 0.26, sig, width=0.26, color=SER[j], edgecolor=SURF, linewidth=2, label=beams)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=30, ha='right', fontsize=8); ax.set_yscale('log')
        corr = np.array(F[k0]['corr']); mc = np.max(np.abs(corr - np.eye(len(names)))) if len(names) > 1 else 0
        ax.set_title(f'{pn}\nmax |corr| joint = {mc:.2f}', loc='left', fontsize=9)
    axs[0].set_ylabel('expected sigma(theta), Asimov'); axs[0].legend(fontsize=8)
    fig.suptitle('Identifiability: expected uncertainty of each class normalisation, FHC alone, RHC alone and joint (full exposure)', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'sensitivity_fisher')
    g = s['gain']['full exposure']
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
    for ax, beam in zip(axs, ('FHC', 'RHC')):
        row = g[beam]; keys = ['prior total', 'conditional (A)'] + [k for k in row if k.startswith('fit ')]
        v = [100 * row[k] for k in keys]
        ax.barh(np.arange(len(keys)), v, color=[INK2] + [SER[0]] + [SER[1]] * (len(keys) - 2), edgecolor=SURF, linewidth=2)
        for i, val in enumerate(v): ax.text(val, i, f' {val:.1f}%', va='center', fontsize=8)
        ax.set_yticks(np.arange(len(keys))); ax.set_yticklabels(keys, fontsize=8); ax.invert_yaxis()
        ax.set_xlabel('relative uncertainty on the one-bin total [%]')
        ax.set_title(f"{beam}: prior flux term {100 * row['prior flux']:.1f}%, B/S = {row['B/S']:.2f}", loc='left', fontsize=10)
    fig.suptitle('Expected effect of the control-region constraint on the one-bin total (Asimov, full exposure)', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'sensitivity_gain')


def fig_detvar():
    import uproot
    z = np.load('/data/uboone/processed/bkgfit/detvar_matched.npz')
    T = M.Templates()
    fig, axs = plt.subplots(2, 2, figsize=(14, 7.5))
    for r, mode in enumerate(('run4fhc', 'run4rhc')):
        ax = axs[r, 0]; f = z[f'{mode}_f']; C = z[f'{mode}_noise']; raw = T.dv[mode]
        x = np.arange(len(M.KNOBS)); tot = [(M.TOT_CC0PI, 'CC0pi total'), (M.TOT_PI0, 'pi0 total'), (M.SR, 'signal region')]
        for j, (b, lab) in enumerate(tot):
            ax.scatter(x + (j - 1) * 0.25, 100 * raw[:, b], marker='x', color=SER[j], s=40, label=f'{lab}: unscaled ratio' if r == 0 else None)
            ax.errorbar(x + (j - 1) * 0.25, 100 * f[:, b], yerr=100 * np.sqrt(C[:, b, b]), fmt='o', color=SER[j], ms=5, capsize=2,
                        label=f'{lab}: POT-normalised, matched (+- bootstrap noise)' if r == 0 else None)
        ax.axhline(0, color=INK2, lw=1); ax.set_xticks(x); ax.set_xticklabels(M.KNOBS, rotation=30, fontsize=8)
        ax.set_ylabel('shift [%]'); ax.set_title(f'{mode}: knob shifts on the totals', loc='left', fontsize=10)
        ax = axs[r, 1]; bins = M.CC0PI + M.PI0
        q = np.sqrt((f[:, bins] ** 2).sum(0)); n = np.sqrt(np.array([C[:, b, b].sum() for b in bins])); c = np.sqrt(np.clip(q ** 2 - n ** 2, 0, None))
        qr = np.sqrt((raw[:, bins] ** 2).sum(0))
        xb = np.arange(len(bins))
        ax.plot(xb, 100 * qr, drawstyle='steps-mid', color=INK2, lw=1.5, ls=':', label='unscaled ratios')
        ax.plot(xb, 100 * q, drawstyle='steps-mid', color=SER[0], label='POT-normalised shift')
        ax.plot(xb, 100 * n, drawstyle='steps-mid', color=SER[1], ls='--', label='bootstrap statistical noise')
        ax.plot(xb, 100 * c, drawstyle='steps-mid', color=SER[2], label='noise subtracted')
        ax.axvline(20.5, color=GRID, lw=1.5); ax.set_ylim(0, 35)
        ax.set_xlabel('control-region bin (CC0pi | pi0; cos theta_mu x p_mu)'); ax.set_ylabel('quadrature over 8 knobs [%]')
        ax.set_title(f'{mode}: per-bin detector term', loc='left', fontsize=10)
    axs[0, 0].legend(fontsize=6, loc='lower left'); axs[0, 1].legend(fontsize=8)
    fig.suptitle('Detector variations: POT normalisation and removal of the statistical part (event-matched Poisson bootstrap)', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'detvar_matched')


def fig_robustness():
    r = load('robustness')
    if not r: return
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.2))
    ax = axs[0]; labels = list(r['variants'])
    for j, beam in enumerate(('FHC', 'RHC')):
        pri = [100 * r['variants'][l][beam]['prior total'] for l in labels]; A = [100 * r['variants'][l][beam]['conditional (A)'] for l in labels]
        y = np.arange(len(labels)) + (j - 0.5) * 0.38
        ax.barh(y, pri, height=0.36, color=GRID, edgecolor=SER[j], linewidth=1.5)
        ax.barh(y, A, height=0.36, color=SER[j], edgecolor=SURF, linewidth=1, label=f'{beam}: constrained (bar outline: prior)')
    ax.set_yticks(np.arange(len(labels))); ax.set_yticklabels([l.replace('A, ', '').replace('A: ', '') for l in labels], fontsize=7); ax.invert_yaxis()
    ax.set_xlabel('one-bin total, relative uncertainty [%]'); ax.legend(fontsize=7, loc='lower right'); ax.set_title('which correlations the constraint uses', loc='left', fontsize=10)
    ax = axs[1]; fams = list(r['closure'])
    for j, beam in enumerate(('FHC', 'RHC')):
        for k, (key, lab) in enumerate((('prior_rel', 'prior'), ('post_rel', 'constrained'))):
            v = [100 * r['closure'][f]['summary'][beam][key]['rms'] for f in fams]
            ax.bar(np.arange(len(fams)) + (2 * j + k - 1.5) * 0.2, v, width=0.2, color=SER[j] if k else GRID, edgecolor=SER[j], linewidth=1.5,
                   label=f'{beam} {lab}')
        for i, f in enumerate(fams):
            ax.text(i + (2 * j - 0.5) * 0.2, 1, f"pull\n{r['closure'][f]['summary'][beam]['post_pull']['rms']:.2f}", ha='center', fontsize=6, color=SURF)
    ax.set_xticks(np.arange(len(fams))); ax.set_xticklabels(['flux universes (30)', 'GENIE universes (30)', 'detector knobs (8)'])
    ax.set_ylabel('RMS bias of the one-bin total [%]'); ax.legend(fontsize=7); ax.set_title('leave-one-out closure (numbers: RMS pull, constrained)', loc='left', fontsize=10)
    fig.tight_layout(); save(fig, 'robustness')


def fig_validate():
    v = load('validate')
    if not v: return
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    for row, bs in enumerate(('totals', '2D')):
        if bs not in v: continue
        vb = v[bs]
        ax = axs[row, 0]; bins = np.linspace(-4, 4, 33); xx = np.linspace(-4, 4, 200)
        for k, (pn, n) in enumerate([('CC0pi+pi0', 'CC0pi'), ('CC0pi+pi0', 'pi0'), ('CC0pi', 'CC0pi')]):
            t = vb['toys'][pn]; names = list(vb['asimov'][pn]['values']); j = names.index(n); p = np.array(t['pulls'])[:, j]
            ax.hist(p, bins=bins, histtype='step', lw=2, color=SER[k], density=True, label=f'{pn}: {n}  mean {t["pull_mean"][j]:+.2f} width {t["pull_width"][j]:.2f} cov {100 * t["coverage"][j]:.0f}%')
        ax.plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), color=INK2, lw=1, ls=':')
        ax.set_xlabel('pull'); ax.legend(fontsize=6); ax.set_title(f'{bs}: toy pulls ({vb["toys"]["CC0pi"]["n"]} toys)', loc='left', fontsize=10)
        ax = axs[row, 1]; inj = vb['injections']; labels = list(inj)
        for k, n in enumerate(('CC0pi', 'pi0')):
            ax.barh(np.arange(len(labels)) + (k - 0.5) * 0.4, [inj[l]['bias_over_sigma'][n] for l in labels], height=0.4, color=SER[k], edgecolor=SURF, label=n)
        ax.set_yticks(np.arange(len(labels))); ax.set_yticklabels(labels, fontsize=7); ax.invert_yaxis(); ax.axvline(0, color=INK2, lw=1)
        ax.set_xlabel('(fitted - true) / sigma'); ax.legend(fontsize=8); ax.set_title(f'{bs}: Asimov injections (two-parameter model)', loc='left', fontsize=10)
        ax = axs[row, 2]; yl = []
        for j, lab in enumerate(vb['mcmc']):
            mc = vb['mcmc'][lab]; mi = mc['minuit']
            for k, n in enumerate(mc['mean']):
                lo, hi = (mi['minos'][n] if mi.get('minos') and n in mi['minos'] else (-mi['errors'][n], mi['errors'][n]))
                y = len(yl); yl.append(f'{lab}: {n}')
                ax.errorbar(mi['values'][n], y + 0.15, xerr=[[-lo], [hi]], fmt='o', color=SER[0], capsize=3, label='MINUIT (MINOS)' if y == 0 else None)
                ax.errorbar(mc['mean'][n], y - 0.15, xerr=[[mc['mean'][n] - mc['q16'][n]], [mc['q84'][n] - mc['mean'][n]]], fmt='s', color=SER[1], capsize=3, label='MCMC 16-84%' if y == 0 else None)
        ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl, fontsize=7); ax.invert_yaxis(); ax.legend(fontsize=8); ax.set_xlabel('theta')
        ax.set_title(f'{bs}: profile likelihood vs MCMC', loc='left', fontsize=10)
    fig.tight_layout(); save(fig, 'validation')


def fig_data():
    d = load('data')
    if not d: return
    r = d['residuals']; labels = r['labels']; data = np.array(r['data']); pre = np.array(r['prefit']); post = np.array(r['postfit']); sig = np.array(r['sigma'])
    nb = len(labels) // 2
    fig, axs = plt.subplots(2, 2, figsize=(15, 7.5))
    for gi, beam in enumerate(('FHC (runs 4, 5)', 'RHC (runs 1, 3, 4)')):
        for col, (reg, rng) in enumerate((('CC0pi', range(0, 21)), ('pi0', range(21, 42)))):
            ax = axs[gi, col]; idx = np.arange(gi * nb, (gi + 1) * nb)[list(rng)]; x = np.arange(len(idx))
            ax.errorbar(x, data[idx] / pre[idx], yerr=np.sqrt(data[idx]) / pre[idx], fmt='o', color=INK, ms=4, label='data / pre-fit')
            ax.plot(x, post[idx] / pre[idx], drawstyle='steps-mid', color=SER[1], label='post-fit / pre-fit (CC0pi, pi0 norms from the totals fit)')
            ax.fill_between(x, 1 - sig[idx] / pre[idx], 1 + sig[idx] / pre[idx], step='mid', color=SER[0], alpha=0.18, label='total uncertainty (corrected detector)')
            for k in range(1, 7): ax.axvline(3 * k - 0.5, color=GRID, lw=1)
            ax.axhline(1, color=INK2, lw=1)
            ax.set_xticks([3 * k + 1 for k in range(7)]); ax.set_xticklabels([f'{M.CTH[k]:g}..{M.CTH[k + 1]:g}' for k in range(7)], fontsize=8)
            ax.set_xlabel('reco cos theta_mu group (within each: p_mu < 0.35, 0.35-0.75, > 0.75 GeV/c)')
            ax.set_title(f'{beam}  {reg}', loc='left', fontsize=10); ax.set_ylim(0.4, 2.0); ax.set_ylabel('ratio to pre-fit prediction')
    axs[0, 0].legend(fontsize=7, loc='upper left')
    fig.suptitle('Beam-on control regions (FHC Run 1 excluded): data and post-fit prediction relative to pre-fit', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'data_residuals')
    # f_FHC_R1 profiles
    fig, ax = plt.subplots(figsize=(6, 4))
    for k, key in enumerate([x for x in d if 'per period + f_FHC_R1' in x and 'profile' in d[x]]):
        pr = d[key]['profile']; g = np.array(pr['grid']); c = np.array(pr['chi2']); ax.plot(g, c - c.min(), color=SER[k], ls=['-', '--', '-.', ':'][k], label=key.replace(' | per period + f_FHC_R1 | ', ': '))
    ax.axhline(1, color=INK2, lw=1, ls=':'); ax.set_ylim(0, 25); ax.set_xlabel('f_FHC_R1 (exposure factor on FHC Run 1)'); ax.set_ylabel('Delta chi2'); ax.legend(fontsize=8)
    ax.set_title('FHC Run-1 exposure hypothesis', loc='left', fontsize=10)
    fig.tight_layout(); save(fig, 'data_fR1_profile')


def fig_ensembles():
    e = load('ensembles')
    if not e: return
    variants = ['standard', 'poisson', 'gauss', 'pearson']
    vlab = {'standard': 'syst. throw + Poisson, CNP', 'poisson': 'Poisson only, CNP', 'gauss': 'syst. throw + Gaussian, CNP', 'pearson': 'syst. throw + Poisson, Pearson'}
    panels = [(bs, pn, k) for bs in ('totals', '2D') for pn in ('CC0pi', 'CC0pi+pi0') for k in range(len(e[f'{bs} | {pn} | standard']['names']))]
    fig, axs = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    for ax, (bs, pn, k) in zip(axs.flat, panels):
        for j, v in enumerate(variants):
            key = f'{bs} | {pn} | {v}'
            if key not in e: continue
            run = e[key]['running']; n = np.array([r['n'] for r in run]); mm = np.array([r['mean'][k] for r in run]); ee = np.array([r['err'][k] for r in run])
            ax.plot(n, mm, color=SER[j], ls=['-', '--', '-.', ':'][j], label=f"{vlab[v]}: {e[key]['pull_mean'][k]:+.3f} +- {e[key]['pull_mean_err'][k]:.3f}")
            if v == 'standard': ax.fill_between(n, mm - ee, mm + ee, color=SER[j], alpha=0.18)
        ax.axhline(0, color=INK2, lw=1); ax.set_xscale('log'); ax.set_ylim(-0.3, 0.3)
        name = e[f'{bs} | {pn} | standard']['names'][k]
        ax.set_title(f'{bs}: {pn} model, {name}', loc='left', fontsize=10); ax.legend(fontsize=6.5, loc='upper right')
        ax.set_ylabel('running pull mean')
    for ax in axs[1]: ax.set_xlabel('number of toys')
    fig.suptitle('Pull-mean convergence with the ensemble size (band: statistical error of the standard ensemble)', x=0.01, ha='left')
    fig.tight_layout(); save(fig, 'ensembles')


if __name__ == '__main__':
    fig_composition(); fig_sensitivity(); fig_detvar(); fig_robustness(); fig_validate(); fig_ensembles(); fig_data()
