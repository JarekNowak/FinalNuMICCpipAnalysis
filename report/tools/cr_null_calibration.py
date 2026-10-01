"""cr_null_calibration.py -- null calibration of the global control-region chi2 over the eight totals.

Pseudo-experiments: prediction ~ N(pred, C_pred), with C_pred built from the correlation matrix and the
sigma(full)/N column of a sb_protocol.C log, minus the Poisson term; data ~ Poisson(prediction); chi2
with the nominal V = C_pred + diag(pred). Writes the summary next to the log and the figure
report/figures/cr_null_calibration.pdf (first produced inline on 2026-09-07; scripted 2026-10-01).

  python3 report/tools/cr_null_calibration.py logs/sidebands_perrun/protocol_dirtfix.log [N]
"""
import os, re, sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figures', 'cr_null_calibration.pdf')


def parse(path):
    rows, corr, chi2 = [], [], None
    lines = open(path).read().splitlines()
    for i, l in enumerate(lines):
        m = re.match(r'^(fhc|rhc)\s+(\S+)\s+([\d.]+)\s+([\d.]+)\s+.*\|\s*([-\d.]+)\s+([-\d.]+)\s+([\d.]+)%\s*\|', l)
        if m: rows.append((m.group(1) + '-' + m.group(2), float(m.group(4)), float(m.group(7)) / 100.))
        m = re.match(r'^GLOBAL normalisation chi2 = ([\d.]+) / 8', l)
        if m: chi2 = float(m.group(1))
        if l.startswith('prediction correlation'):
            corr = [[float(x) for x in lines[i + 1 + k].split()[-8:]] for k in range(8)]
    if len(rows) != 8 or len(corr) != 8 or chi2 is None: sys.exit(f'cannot parse {path}')
    return rows, np.array(corr), chi2


def main():
    log = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    rows, rho, chi2_obs = parse(log)
    pred = np.array([r[1] for r in rows]); frac = np.array([r[2] for r in rows])
    sig_pred = np.sqrt(np.maximum((frac * pred) ** 2 - pred, 0.))       # remove the Poisson term
    C = rho * np.outer(sig_pred, sig_pred)
    Vinv = np.linalg.inv(C + np.diag(pred))
    rng = np.random.default_rng(20260907)
    mu = np.clip(rng.multivariate_normal(pred, C, size=n), 0., None)
    d = rng.poisson(mu) - pred
    chi2 = np.einsum('ij,jk,ik->i', d, Vinv, d)
    q95 = np.quantile(chi2, 0.95); c95 = stats.chi2.ppf(0.95, 8)
    p = stats.chi2.sf(chi2, 8)
    below = (chi2 < chi2_obs).mean()
    out = os.path.splitext(log)[0] + '_null_calibration.txt'
    with open(out, 'w') as f:
        f.write(f'null calibration of the global normalisation chi2 (8 totals), {n} pseudo-experiments\n')
        f.write(f'source: {log}\n')
        f.write('model: prediction ~ N(pred, C_pred) with C_pred from the log correlation and sigma(full)/N minus the '
                'Poisson term; data ~ Poisson(prediction); chi2 with the nominal V = C_pred + diag(pred)\n')
        f.write('prediction fractional uncertainties: ' + ' '.join(f'{r[0]} {100 * s / v:.1f}%' for r, s, v in zip(rows, sig_pred, pred)) + '\n')
        f.write(f'mean chi2 = {chi2.mean():.3f} (nominal 8); 95% quantile = {q95:.2f} (chi2_8 {c95:.2f}); '
                f'fraction with nominal p<0.05 = {100 * (p < 0.05).mean():.2f}% (5%); p<0.01 = {100 * (p < 0.01).mean():.2f}% (1%)\n')
        f.write(f'observed pseudo-data chi2 {chi2_obs:.2f}: {100 * below:.2f}% of null pseudo-experiments lie below it '
                f'(chi2_8: {100 * stats.chi2.cdf(chi2_obs, 8):.2f}%)\n')
        f.write('limitation: the covariance is treated as known (finite-universe and detector-unisim sampling noise not re-thrown)\n')
    print(open(out).read())

    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    b = np.linspace(0, 30, 121)
    ax.hist(chi2, bins=b, density=True, histtype='step', color='k', label='pseudo-experiments under the null')
    x = np.linspace(0, 30, 400)
    ax.plot(x, stats.chi2.pdf(x, 8), 'r-', label=r'$\chi^2_8$')
    ax.axvline(chi2_obs, color='b', ls='--', label=f'full-stack pseudo-data, {chi2_obs:.2f}')
    ax.axvline(c95, color='r', ls=':', label=f'$p=0.05$ ({c95:.1f})')
    ax.set_xlabel(r'global normalisation $\chi^2$ over the eight totals'); ax.set_ylabel('density')
    ax.set_xlim(-1, 31); ax.set_ylim(bottom=0); ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(FIG); print('wrote', os.path.abspath(FIG))


if __name__ == '__main__':
    main()
