"""mcs_cal_release.py -- the MCS term of the release from version 1.9: the calibration uncertainty of the MCS momentum of
exiting muons (xsec_analyzer/scripts/mcs_cal_compare.py: model, grid, production and beam variations of the
control-region calibration, fully correlated terms) replaces the data-side +-5% MCS momentum-scale term of earlier
versions (report/tools/mcs_eval.py).

For every result in data_release/cov/: cov_MCScal.txt is the calibration term (the results that use the muon momentum,
from /data/uboone/processed/mcscal/nominal_out/<tag>/cov_MCScal.txt; identically zero for the others),
cov_total_plusMCS.txt = cov_total + cov_MCScal (the official covariance, report/tools/official_cov.py), and
cov_MCSscale.txt is removed. The calibrated unfolded result of the release must equal the one the term was computed
around. Summary: data_release/mcs_cal.tsv (bin-averaged term and its four components, percent of the result).

    python3 report/tools/mcs_cal_release.py
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.join(os.path.dirname(HERE), 'data_release')
CAL = '/data/uboone/processed/mcscal/'
CMP = os.path.join(os.path.dirname(HERE), 'planning', 'multipion', 'mcs_cal_compare.json')


def vec(p):
    return np.array([float(a[1]) for a in (l.split() for l in open(p)) if len(a) == 2 and a[0] != 'numXbins'])


def cov(p):
    r = [l.split() for l in open(p)]; n = int(r[0][1]); m = np.zeros((n, n))
    for a in r[3:]:
        if len(a) == 3: m[int(a[0]), int(a[1])] = float(a[2])
    return m


def write(p, C):
    n = len(C)
    with open(p, 'w') as o:
        o.write(f'numXbins {n}\nnumYbins {n}\nxbin  ybin  z\n')
        for i in range(n):
            for j in range(n): o.write(f'{i}  {j}  {C[i, j]:.17e}\n')


def main():
    terms = json.load(open(CMP)) if os.path.exists(CMP) else {}
    rows, ncal = [], 0
    for tag in sorted(os.listdir(os.path.join(REL, 'cov'))):
        d = os.path.join(REL, 'cov', tag)
        x = vec(os.path.join(d, 'unfolded_signal.txt')); T = cov(os.path.join(d, 'cov_total.txt'))
        src = CAL + 'nominal_out/' + tag + '/'
        if os.path.exists(src + 'cov_MCScal.txt'):
            xc = vec(src + 'unfolded_signal.txt')
            if not np.allclose(x, xc, rtol=1e-9, atol=0):
                raise SystemExit(f'{tag}: the release result is not the calibrated one the term was computed around')
            C = cov(src + 'cov_MCScal.txt'); ncal += 1
        else:
            C = np.zeros_like(T)
        write(os.path.join(d, 'cov_MCScal.txt'), C)
        write(os.path.join(d, 'cov_total_plusMCS.txt'), T + C)
        old = os.path.join(d, 'cov_MCSscale.txt')
        if os.path.exists(old): os.remove(old)
        s = np.sqrt(np.diag(T)); f = lambda v: 100 * float(np.mean(np.abs(v) / np.abs(x)))
        comp = terms.get(tag, {}).get('terms', {})
        rows.append([tag, f(np.sqrt(np.diag(C))), f(np.sqrt(np.diag(T + C))) - f(s)] +
                    [f(np.array(comp[k])) if k in comp else 0. for k in ('model', 'grid', 'production', 'beam')])
    with open(os.path.join(REL, 'mcs_cal.tsv'), 'w') as o:
        o.write('# MCS calibration term (report/tools/mcs_cal_release.py), bin-averaged percent of the result: the term, the\n'
                '# increase of the official total over the framework total, and the four components (shifts)\n')
        o.write('tag\tMCScal_pct\ttotal_increase_pct\tmodel_pct\tgrid_pct\tproduction_pct\tbeam_pct\n')
        for r in rows: o.write(r[0] + '\t' + '\t'.join(f'{v:.2f}' for v in r[1:]) + '\n')
    print(f'{len(rows)} results, {ncal} with a calibration term; data_release/mcs_cal.tsv written')


if __name__ == '__main__':
    main()
