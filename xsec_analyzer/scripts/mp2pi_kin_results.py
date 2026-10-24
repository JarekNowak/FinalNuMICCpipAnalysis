"""mp2pi_kin_results.py -- the two-pion differential extractions through the framework (branch multipion/kin:
slurm/slurm_mp2pi_kin_extract.sbatch, outputs /data/uboone/processed/mp2pi_kin/unfold/<variable>/), summarised per
variable, configuration and constraint mode: the unfolded pseudo-data per bin with its total uncertainty (log lines
"Bin: i, Val, Unc" of the physical slice), the closure chi2 against the A_C-smeared truth of the pseudo-data and against
the tune (the framework's chi2 with the total covariance), the bin-averaged uncertainty (SYSTDUMP: mean of the fractional
bin uncertainties, total and data statistics) and, with the constraint, the uncertainty of the integral over the bins
from the full covariance (unfold_output of the configuration, written by the constrained run, the last of the two;
the unconstrained combined run is repeated into <v>/comb_none). Combined configuration: also with the soft-pion
efficiency term (softpion_cov_<none|sidebands>.txt of scripts/mp2pi_softpion.py cov) added to the covariance, per bin
(rel_soft, unc_soft) and on the integral (int_total_soft); the figure shows the uncertainties with it.
Generator predictions (GENIE, GiBUU, NEUT, NuWro of scripts/mp2pi_gen/mp2pi_gen.py, A_C-smeared by the framework, closure
sidecar h_gen_*): the framework's chi2 (unfolded covariance) and, combined, the chi2 with the covariance including the
soft-pion term; one-bin total: generator over pseudo-data result, in units of the total uncertainty with the soft-pion
and MCS-calibration terms. Figure: the combined configuration with the constraint, every variable, against the truth and the tune.

    python3 scripts/mp2pi_kin_results.py [--prod 2f]  ->  report/planning/multipion/mp2pi_kin_results.{md,json},
                                              report/figures/multipion/kin_results_comb.pdf
"""
import json, os, re
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
FIG = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'multipion'))
import sys
PROD = sys.argv[sys.argv.index('--prod') + 1] if '--prod' in sys.argv else ''   # 2f: the two-fold production
U = f'/data/uboone/processed/mp2pi_{PROD or "kin"}/unfold/'
VARS = {'pmu': r'$p_\mu$ [GeV/$c$]', 'cmu': r'$\cos\theta_\mu$', 'cfwd': r'$\cos\theta_\pi^{\rm fwd}$',
        'cbwd': r'$\cos\theta_\pi^{\rm bwd}$', 'cpipi': r'$\cos\theta_{\pi\pi}$', 'cmupi': r'$\cos\theta_{\mu\pi}$',
        'np': r'$N_p$ (0, $\geq1$)'}
CFGS = {'fhc5': 'FHC', 'rhcfull': 'RHC', 'comb': 'combined'}
MODES = {'': 'none', '_constr': 'sidebands'}


def read_tab(path):
    rows = [l.split() for l in open(path) if l.strip() and l.split()[0][0].isdigit()]
    if len(rows[0]) == 2: return np.array([float(r[1]) for r in rows])
    n = int(max(int(r[0]) for r in rows)) + 1; m = np.zeros((n, n))
    for r in rows: m[int(r[0]), int(r[1])] = float(r[2])
    return m


def integral(v, cfg):
    d = f'{U}{v}/{cfg}/unfold_output/'
    x = read_tab(d + 'vec_table_unfolded_signal.txt').sum()
    return [100 * np.sqrt(read_tab(d + f'mat_table_cov_{c}.txt').sum()) / x for c in ('total', 'DataStats')]


GENS = ['GENIE', 'GiBUU', 'NEUT', 'NuWro']
GCOL = {'GENIE': '#7b4fa0', 'GiBUU': '#2a9d58', 'NEUT': '#e9a23b', 'NuWro': '#1f5fa8'}


def soft(v, cfg, mode, gens=None):
    """per-bin relative uncertainty, per-bin uncertainty per unit x and integral uncertainty with the soft-pion term;
    chi2 of each generator (A_C-smeared, per unit x) with that covariance"""
    S = f'{U}{v}/softpion_cov_{MODES[mode]}.txt'
    d = f'{U}{v}/' + ('comb_none' if mode == '' else 'comb') + '/unfold_output/'
    if cfg != 'comb' or not os.path.exists(S) or not os.path.exists(d + 'mat_table_cov_total.txt'): return {}
    x = read_tab(d + 'vec_table_unfolded_signal.txt'); n = len(x)
    C = read_tab(d + 'mat_table_cov_total.txt')[:n, :n] + np.atleast_2d(np.loadtxt(S))
    w = np.diff(uproot.open(f'{U}{v}/closure_hists_xsec_{cfg}{mode}.root')['h_unfolded_nuwro'].axis().edges())
    Ci = np.linalg.inv(C)
    chi = {g: float((x - np.array(h) * w) @ Ci @ (x - np.array(h) * w)) for g, h in (gens or {}).items()}
    return dict(rel_soft=(100 * np.sqrt(np.diag(C)) / np.abs(x)).tolist(), unc_soft=(np.sqrt(np.diag(C)) / w).tolist(),
                int_total_soft=float(100 * np.sqrt(C.sum()) / x.sum()), avg_soft=float(np.mean(100 * np.sqrt(np.diag(C)) / np.abs(x))),
                chi2_gen_soft=chi, gen_over_x={g: float((np.array(h) * w).sum() / x.sum()) for g, h in (gens or {}).items()})


def parse(v, cfg, mode):
    log = f'{U}{v}/{cfg}{mode}.log'
    if not os.path.exists(log): return None
    t = open(log).read()
    sd = {k: float(x) for k, x in re.findall(r'^\[SYSTDUMP\] (\S+) ([-\d.e+]+)$', t, re.M)}
    blk = t[t.index('[closure-dump]'):]
    chi_t = re.search(r'^truth: χ² = ([\d.e+-]+)/(\d+) bins', blk, re.M)
    chi_c = re.search(r'^MicroBooNE Tune: χ² = ([\d.e+-]+)/(\d+) bins', blk, re.M)
    first = blk[:blk.index('Info in <TCanvas::Print>')]
    bins = [(float(a), float(b)) for a, b in re.findall(r'^Bin: \d+, Val = ([\d.e+-]+), Unc = ([\d.e+-]+)', first, re.M)]
    f = uproot.open(f'{U}{v}/closure_hists_xsec_{cfg}{mode}.root')
    gens = {g: f[f'h_gen_{g}'].values().tolist() for g in GENS if f'h_gen_{g}' in f}
    chig = {}
    for g in gens:
        m = re.search(rf'^{g}: χ² = ([\d.e+-]+)/', blk, re.M)
        if m: chig[g] = float(m.group(1))
    return dict(gens=gens, chi2_gen=chig, val=[b[0] for b in bins], unc=[b[1] for b in bins], truth=f['h_fakedata_truth'].values().tolist(),
                tune=f['h_genie_tune'].values().tolist(), edges=f['h_unfolded_nuwro'].axis().edges().tolist(),
                chi2_truth=float(chi_t.group(1)), ndf=int(chi_t.group(2)), chi2_tune=float(chi_c.group(1)),
                avg_total=sd['total'], avg_stat=sd['DataStats'],
                **(dict(zip(('int_total', 'int_stat'), integral(v, cfg))) if mode == '_constr' else {}), **soft(v, cfg, mode, gens))


def totals():
    """one-bin total: pseudo-data result with the soft-pion and MCS-calibration terms against the generators"""
    if PROD != '2f' or not os.path.exists(f'{U}total/comb_constr.log'): return []
    SOFT = {c: 100 * r['term'] for c, r in json.load(open(os.path.join(PLAN, 'mp2pi_softpion_term.json')))['total'].items()}
    MCS = {'fhc5': 5.85, 'rhcfull': 5.88, 'comb': 4.72}   # constrained, mp2pi_mcscal_var_2f.md (branch multipion/kin)
    L = ['', '## One-bin total against the generators (constrained; total uncertainty with the soft-pion and MCS-calibration '
         'terms; generator / pseudo-data and (generator - pseudo-data) / uncertainty)', '',
         '| configuration | pseudo-data | uncertainty | ' + ' | '.join(GENS) + ' |', '|---|---|---|' + '---|' * len(GENS)]
    for cfg, nm in CFGS.items():
        t = open(f'{U}total/{cfg}_constr.log').read()
        sd = {k: float(x) for k, x in re.findall(r'^\[SYSTDUMP\] (\S+) ([-\d.e+]+)$', t, re.M)}
        x = sd['sigma_int']; u = np.sqrt(sd['total'] ** 2 + SOFT[cfg] ** 2 + MCS[cfg] ** 2) / 100
        gm = {'fhc5': 'fhc', 'rhcfull': 'rhc', 'comb': 'comb'}[cfg]
        gv = {g: float(uproot.open(f'/data/uboone/processed/mp2pi_2f/gen/{g.lower()}_2pi_{gm}_fte.root')['total_fte'].values()[0]) for g in GENS}
        L.append(f'| {nm} | {x:.4f} | {100 * u:.1f}% | ' + ' | '.join(f'{gv[g]:.4f}: {gv[g] / x:.2f} ({(gv[g] - x) / (u * x):+.1f})' for g in GENS) + ' |')
    return L


def main():
    res = {}
    for v in VARS:
        for cfg in CFGS:
            for mode in MODES:
                r = parse(v, cfg, mode)
                if r: res[f'{v} {cfg} {MODES[mode]}'] = r
    json.dump(res, open(os.path.join(PLAN, f'mp2pi_kin{PROD}_results.json'), 'w'), indent=1)
    L = ['# Two-pion differential extractions through the framework (pseudo-data)', '',
         'Generated by `scripts/mp2pi_kin_results.py` from /data/uboone/processed/mp2pi_kin/unfold (branch multipion/kin). '
         'Wiener-SVD (filter on, second derivative), detector variations fractional and flux-weighted, MCS-calibrated '
         'classifier inputs; bins of `mp2pi_kinematics.md`. Per-bin uncertainty: total, relative to the unfolded value. '
         'Closure: chi2 of the unfolded pseudo-data against their A_C-smeared truth with the total covariance; ratio = '
         'unfolded / truth per bin. Bin average: mean of the fractional bin uncertainties (data statistics). Integral: '
         'uncertainty of the sum over the bins from the full covariance (data statistics), with the constraint.', '',
         '| variable | configuration | constraint | bins | per-bin uncertainty | bin average (data stat.) | integral (data stat.) | with soft-pion term: per bin, average, integral | closure chi2 | unfolded / truth | chi2 vs tune |',
         '|---|---|---|---|---|---|---|---|---|---|---|']
    for k, r in res.items():
        v, cfg, mode = k.split()
        rel = 100 * np.array(r['unc']) / np.abs(np.array(r['val'])); ratio = np.array(r['val']) / np.array(r['truth'])
        L.append(f"| {v} | {CFGS[cfg]} | {mode} | {r['ndf']} | {rel.min():.0f}-{rel.max():.0f}% | {r['avg_total']:.1f}% "
                 f"({r['avg_stat']:.1f}%) | " + (f"{r['int_total']:.1f}% ({r['int_stat']:.1f}%)" if 'int_total' in r else '') + ' | ' +
                 (f"{min(r['rel_soft']):.0f}-{max(r['rel_soft']):.0f}%, {r['avg_soft']:.1f}%, {r['int_total_soft']:.1f}%" if 'rel_soft' in r else '') +
                 f" | {r['chi2_truth']:.2f}/{r['ndf']} | {ratio.min():.2f}-{ratio.max():.2f} | "
                 f"{r['chi2_tune']:.2f}/{r['ndf']} |")
    if any(r.get('chi2_gen_soft') for r in res.values()):
        L += ['', '## Generator predictions (combined; chi2 with the covariance including the soft-pion term, framework chi2 in '
              'brackets; ratio: integral over the bins against the unfolded pseudo-data)', '',
              '| variable | constraint | bins | ' + ' | '.join(GENS) + ' |', '|---|---|---|' + '---|' * len(GENS)]
        for k, r in res.items():
            v, cfg, mode = k.split()
            if not r.get('chi2_gen_soft'): continue
            L.append(f"| {v} | {mode} | {r['ndf']} | " + ' | '.join(
                f"{r['chi2_gen_soft'][g]:.1f} ({r['chi2_gen'].get(g, float('nan')):.1f}), {r['gen_over_x'][g]:.2f}" for g in GENS) + ' |')
    L += totals()
    open(os.path.join(PLAN, f'mp2pi_kin{PROD}_results.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    # figure: combined, constrained, all variables
    fig, axs = plt.subplots(2, 4, figsize=(13.5, 6.2))
    for ax, v in zip(axs.flat, VARS):
        r = res.get(f'{v} comb sidebands')
        if not r: ax.axis('off'); continue
        e = np.array(r['edges']); c = 0.5 * (e[1:] + e[:-1]); w = np.diff(e)
        if v == 'np': e = np.array([-0.5, 0.5, 1.5]); c = np.array([0., 1.]); w = np.ones(2)
        ax.stairs(r['truth'], e, color='#2a9d58', lw=1.6, label='truth of the pseudo-data ($A_C$)')
        ax.stairs(r['tune'], e, color='#c8553d', lw=1.3, ls='--', label='MicroBooNE tune ($A_C$)')
        for g in GENS:
            if g in r['gens']: ax.stairs(r['gens'][g], e, color=GCOL[g], lw=1.0, label=f'{g} ($A_C$)')
        ax.errorbar(c, r['val'], yerr=r.get('unc_soft', r['unc']), xerr=w / 2, fmt='o', color='k', ms=3.5, lw=1, label='unfolded pseudo-data')
        ax.set_xlabel(VARS[v]); ax.set_ylim(0, 1.45 * max(max(r['val']) + max(r.get('unc_soft', r['unc'])), max(r['tune'])))
        ax.text(0.03, 0.95, f"$\\chi^2$ truth {r['chi2_truth']:.1f}/{r['ndf']}\nintegral {r.get('int_total_soft', r['int_total']):.0f}%",
                transform=ax.transAxes, fontsize=8.5, va='top')
        if v == 'np': ax.set_xticks([0, 1]); ax.set_xticklabels(['0', r'$\geq1$'])
    axs[0, 0].set_ylabel(r'$d\sigma/dx$ [$10^{-38}$ cm$^2$/Ar per unit $x$]')
    axs[1, 0].set_ylabel(r'$d\sigma/dx$ [$10^{-38}$ cm$^2$/Ar per unit $x$]')
    h, l = axs[0, 0].get_legend_handles_labels(); axs[1, 3].axis('off'); axs[1, 3].legend(h, l, loc='center', frameon=False)
    fig.text(0.01, 0.995, 'MicroBooNE NuMI simulation, pseudo-data, combined FHC+RHC, sideband constraint, internal',
             fontsize=9, color='0.45', va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(FIG, f'kin{PROD}_results_comb.pdf')); plt.close(fig)


if __name__ == '__main__':
    main()
