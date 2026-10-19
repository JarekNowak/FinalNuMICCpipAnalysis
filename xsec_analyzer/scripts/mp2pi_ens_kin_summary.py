"""mp2pi_ens_kin_summary.py -- pseudo-experiment ensembles of the two-pion differential extractions (two-fold
production; slurm/slurm_mp2pi_ens.sbatch with KINVARS on branch multipion/kin, output
/data/uboone/processed/mp2pi_2f/ens/<ENS>_kin/m<NNN>/<v>/). Members are the throws of the one-bin ensemble of the
same name (same seeds), extracted per variable with the nominal response, combined configuration, without and with
the sideband constraint.

Per member and extraction (bin-integrated, the units of the framework's covariance tables): x, the unfolded vector,
its total and data-statistics covariance; t, the member's own truth smeared by A_C; e, the expected truth, the
A_C-smeared tune from which the nominal members are thrown (closure sidecar, per unit x, times the bin width).
  - bias per bin: mean of x/t - 1, and mean(x - t) in units of the total uncertainty (every ensemble);
  - coverage per bin (nominal ensembles): pull (x - e)/sqrt(C_stat), its width and the share inside the 68%
    interval; chi2_stat = (x - e)^T C_stat^-1 (x - e), whose mean is the number of bins for a correct C_stat;
  - the framework's closure chi2 against the truth (parsed from the member's log, the statistic of Table 16 of the
    note; it adds the Poisson variance of the true count to the unfolded covariance) and the empirical p-value of
    the nominal realisation (report/planning/multipion/mp2pi_kin2f_results.json);
  - the integral over the bins: offset x/t and pull against e with the summed data-statistics covariance.

    python3 scripts/mp2pi_ens_kin_summary.py [ENS ...]   (default nominal2)
        ->  report/planning/multipion/mp2pi_ens_kin_<ENS>.{md,json}, report/figures/multipion/ens_kin_<ENS>.pdf
"""
import glob, json, os, re, sys
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
FIG = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'multipion'))
E = '/data/uboone/processed/mp2pi_2f/ens/'
VARS = {'pmu': r'$p_\mu$', 'cmu': r'$\cos\theta_\mu$', 'cfwd': r'$\cos\theta_\pi^{\rm fwd}$', 'cbwd': r'$\cos\theta_\pi^{\rm bwd}$',
        'cpipi': r'$\cos\theta_{\pi\pi}$', 'cmupi': r'$\cos\theta_{\mu\pi}$', 'np': r'$N_p$'}
MODES = {'': 'none', '_constr': 'sidebands'}


def tab(path):
    rows = [l.split() for l in open(path) if l.strip() and l.split()[0][0].isdigit()]
    if len(rows[0]) == 2: return np.array([float(r[1]) for r in rows])
    n = max(int(r[0]) for r in rows) + 1; m = np.zeros((n, n))
    for r in rows: m[int(r[0]), int(r[1])] = float(r[2])
    return m


def member(d, v, mode):
    p = f'{d}/{v}/comb{mode}_'
    if not os.path.exists(p + 'mat_table_cov_DataStats.txt'): return None
    f = uproot.open(f'{d}/{v}/closure_hists_xsec_comb{mode}.root')
    w = np.diff(f['h_fakedata_truth'].axis().edges())
    c = re.search(r'^truth: χ² = ([-\d.e+]+)/', open(f'{d}/{v}/comb{mode}.log').read(), re.M)
    return dict(x=tab(p + 'vec_table_unfolded_signal.txt'), t=f['h_fakedata_truth'].values() * w,
                e=f['h_genie_tune'].values() * w, Ct=tab(p + 'mat_table_cov_total.txt'),
                Cs=tab(p + 'mat_table_cov_DataStats.txt'), chi2_fw=float(c.group(1)) if c else np.nan)


def summarise(ens):
    nominal = json.load(open(os.path.join(PLAN, 'mp2pi_kin2f_results.json')))
    mem = sorted(glob.glob(f'{E}{ens}_kin/m[0-9][0-9][0-9]'))
    res = {}
    for v in VARS:
        for mode, ml in MODES.items():
            M = [r for r in (member(d, v, mode) for d in mem) if r is not None]
            if len(M) < 2: continue
            n = len(M)
            x, t, e = (np.array([r[k] for r in M]) for k in ('x', 't', 'e'))
            st = np.array([np.sqrt(np.diag(r['Ct'])) for r in M]); ss = np.array([np.sqrt(np.diag(r['Cs'])) for r in M])
            off = x / t - 1; bias = (x - t) / st; pull = (x - e) / ss
            chis = np.array([(r['x'] - r['e']) @ np.linalg.solve(r['Cs'], r['x'] - r['e']) for r in M])
            chifw = np.array([r['chi2_fw'] for r in M])
            xi, ti, ei = x.sum(1), t.sum(1), e.sum(1); si = np.array([np.sqrt(r['Cs'].sum()) for r in M])
            nom = nominal.get(f'{v} comb {ml}', {}).get('chi2_truth')
            res[f'{v} {ml}'] = dict(
                members=n, bins=x.shape[1], offset=off.mean(0).tolist(), offset_se=(off.std(0, ddof=1) / np.sqrt(n)).tolist(),
                bias_sigma=bias.mean(0).tolist(), pull_mean=pull.mean(0).tolist(), pull_width=pull.std(0, ddof=1).tolist(),
                cover68=(np.abs(pull) < 1).mean(0).tolist(), chi2_stat_mean=float(chis.mean()), chi2_fw=chifw.tolist(),
                int_offset=float((xi / ti - 1).mean()), int_offset_se=float((xi / ti - 1).std(ddof=1) / np.sqrt(n)),
                int_pull_width=float(((xi - ei) / si).std(ddof=1)), chi2_nominal=nom,
                p_empirical=(float((chifw >= nom).mean()) if nom is not None else None))
    json.dump(res, open(os.path.join(PLAN, f'mp2pi_ens_kin_{ens}.json'), 'w'), indent=1)
    nominal_ens = ens.startswith('nominal')
    L = [f'# Pseudo-experiment ensembles of the two-pion differential extractions ({ens})', '',
         'Generated by `scripts/mp2pi_ens_kin_summary.py` (see its docstring). Combined configuration, two-fold production, '
         'nominal response. Offset: extracted over the member\'s A_C-smeared truth, per bin (range over the bins, and the '
         'largest in units of its standard error); bias: the same difference in units of the total uncertainty (largest '
         'over the bins). ' + ('Pulls against the expected truth with the data-statistics covariance; chi2_stat mean: '
         'expected to equal the number of bins; framework chi2: the closure statistic of Table 16, p = share of members at '
         'or above the nominal realisation.' if nominal_ens else ''), '']
    if nominal_ens:
        L += ['| variable | constraint | members | offset per bin | max |offset|/SE | pull width (range) | 68% coverage (range) | chi2_stat mean / bins | integral offset | integral pull width | framework chi2: median, nominal, p |',
              '|---|---|---|---|---|---|---|---|---|---|---|']
    else:
        L += ['| variable | constraint | members | offset per bin | max |offset|/SE | largest |bias| / total unc. | integral offset |',
              '|---|---|---|---|---|---|---|']
    for k, r in res.items():
        v, ml = k.split()
        o = np.array(r['offset']); se = np.array(r['offset_se']); b = np.array(r['bias_sigma'])
        head = f"| {v} | {ml} | {r['members']} | {100 * o.min():+.1f} to {100 * o.max():+.1f}% | {np.abs(o / se).max():.1f} | "
        if nominal_ens:
            pw = np.array(r['pull_width']); cv = np.array(r['cover68'])
            L.append(head + f"{pw.min():.2f}-{pw.max():.2f} | {100 * cv.min():.0f}-{100 * cv.max():.0f}% | "
                     f"{r['chi2_stat_mean']:.2f} / {r['bins']} | {100 * r['int_offset']:+.1f} +- {100 * r['int_offset_se']:.1f}% | "
                     f"{r['int_pull_width']:.2f} | {np.nanmedian(r['chi2_fw']):.2f}, "
                     + (f"{r['chi2_nominal']:.2f}, {r['p_empirical']:.2f} |" if r['chi2_nominal'] is not None else '|'))
        else:
            L.append(head + f"{np.abs(b).max():.2f} | {100 * r['int_offset']:+.1f} +- {100 * r['int_offset_se']:.1f}% |")
    tf = os.path.join(PLAN, 'mp2pi_softpion_term.json')
    if ens == 'stress_psoft' and os.path.exists(tf):
        term = json.load(open(tf))['differential']
        L += ['', '## Softest-pion stress against the soft-pion efficiency term (generator envelope, per true bin)', '',
              '| variable | stress offset, none | stress offset, sidebands | soft-pion term | largest |offset| / term |', '|---|---|---|---|---|']
        for v in VARS:
            if f'{v} none' not in res or v not in term: continue
            o0, o1 = np.array(res[f'{v} none']['offset']), np.array(res[f'{v} sidebands']['offset'])
            tm = np.array(term[v]['per_bin'])
            L.append(f'| {v} | ' + ', '.join(f'{100 * x:+.1f}' for x in o0) + ' | ' + ', '.join(f'{100 * x:+.1f}' for x in o1)
                     + ' | ' + ', '.join(f'{100 * x:.1f}' for x in tm) + f' | {(np.maximum(np.abs(o0), np.abs(o1)) / tm).max():.2f} |')
    open(os.path.join(PLAN, f'mp2pi_ens_kin_{ens}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    if nominal_ens: figure(ens, res)
    return res


def figure(ens, res):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, len(VARS), figsize=(2.3 * len(VARS), 5.4))
    for j, v in enumerate(VARS):
        for ml, col, dx in (('none', '#c8553d', -0.12), ('sidebands', '#1f5fa8', 0.12)):
            r = res.get(f'{v} {ml}')
            if not r: continue
            b = np.arange(r['bins'])
            axs[0, j].errorbar(b + dx, r['pull_width'], yerr=np.array(r['pull_width']) / np.sqrt(2 * (r['members'] - 1)),
                               fmt='o', ms=3, color=col, label=ml)
            axs[1, j].hist(r['chi2_fw'], bins=20, histtype='step', color=col, lw=1.3, label=ml)
            if r['chi2_nominal'] is not None: axs[1, j].axvline(r['chi2_nominal'], color=col, ls='--', lw=1)
        axs[0, j].axhline(1, color='0.5', lw=0.8); axs[0, j].set_ylim(0.4, 1.4); axs[0, j].set_title(VARS[v], fontsize=11)
        axs[0, j].set_xlabel('bin'); axs[1, j].set_xlabel(r'closure $\chi^2$')
    axs[0, 0].set_ylabel('pull width (data stat.)'); axs[1, 0].set_ylabel('members'); axs[0, 0].legend(fontsize=8)
    fig.text(0.01, 0.995, f'MicroBooNE NuMI simulation, pseudo-experiments ({ens}), combined, internal; dashed: nominal realisation',
             fontsize=8, color='0.45', va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(FIG, f'ens_kin_{ens}.pdf')); plt.close(fig)


if __name__ == '__main__':
    for ens in (sys.argv[1:] or ['nominal2']):
        summarise(ens)
