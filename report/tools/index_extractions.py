# index_extractions.py -- regenerate data_release/index_extractions.tsv from the current release. One row
# per released differential extraction (every directory of data_release/cov/), with the status of
# report/status.tsv in the four-label vocabulary (decision of 2026-09-27), so the file index cannot disagree
# with the notes. Uncertainties are those of the official covariance (framework + MCS term, official_cov.py).
# The one-bin totals are indexed separately in total_xsec.tsv.
import csv, os, glob, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import official_cov as oc
import numpy as np

R = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; D = R + 'data_release/'
LOGS = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
LABEL = {'primary': 'Primary result', 'secondary': 'Secondary result',
         'validation': 'Validation-only product', 'notreported': 'Not reported'}
CF = {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}


def status_map():
    lines = [l for l in open(R + 'status.tsv') if l.strip() and not l.startswith('#')]
    m = {}
    for r in csv.DictReader(lines, delimiter='\t'):
        for key in r['observables'].split(','):
            if key not in ('-', 'total'):
                m[(r['family'], key)] = (r, LABEL[r['status']])
    return m


def closure_from_log(tag, nbins):
    """closure chi2/p of an extraction outside the systematics-dump set, from its harvest log (the line that
    covers every bin; two-dimensional logs also carry one line per outer slice)."""
    best = None
    for line in open(LOGS + f'harvest_proposed_{tag}.log', errors='replace'):
        m = re.match(r'^truth: .* = ([0-9.e+-]+)/(\d+) bins?, p-value = ([0-9.e+-]+)', line)
        if m and int(m.group(2)) == nbins:
            best = (float(m.group(1)), float(m.group(3)))
    return best


def main():
    sm = status_map()
    res = {}
    for fn in ('current_results.tsv', 'current_results_1p_new.tsv'):
        for r in csv.DictReader(open(R + fn), delimiter='\t'):
            res[(r['family'], CF[r['config']], r['observable'])] = r
    rows = []
    for cdir in sorted(glob.glob(D + 'cov/*')):
        tag = os.path.basename(cdir); fam, cfg, obs = tag.split('_', 2)
        cur = D + f'curves_{tag}.tsv'; ac = D + f'A_C_{tag}.tsv'
        assert os.path.exists(cur) and os.path.exists(ac), tag
        prod, st = sm[(fam, obs)]
        if prod[{'FHC5': 'FHC', 'RHCFULL': 'RHC', 'COMB': 'COMB'}[cfg]] != 'Y' and st != LABEL['notreported']:
            st = LABEL['notreported']            # e.g. FHC and RHC extractions of a combined-only result
        nb = sum(1 for l in open(cur) if l[0].isdigit()); nc = len(glob.glob(cdir + '/*.txt'))
        x = oc.vec(cdir + '/unfolded_signal.txt')
        if (fam, cfg, obs) in res:
            r = res[(fam, cfg, obs)]
            sig, tot, det = float(r['sigma_int']), r['PredTotal_pct'], r['detVar_pct']
            chi2, p = float(r['chi2_truth']), float(r['p_truth'])
        else:
            M = oc.cov(cdir + '/cov_MCSscale.txt') if os.path.exists(cdir + '/cov_MCSscale.txt') else 0
            P = oc.cov(cdir + '/cov_PredTotal.txt') + M
            sig = float(x.sum())
            tot = f'{oc._binavg(P, x):.1f}'; det = f"{oc._binavg(oc.cov(cdir + '/cov_detVar_total.txt'), x):.1f}"
            chi2, p = closure_from_log(tag, nb)
        rows.append([tag, fam, cfg, obs, str(nb), f'{sig:.4e}', tot, det, f'{chi2:.3f}', f'{p:.3f}',
                     os.path.basename(ac), os.path.basename(cur), str(nc), st])
    with open(D + 'index_extractions.tsv', 'w') as o:
        o.write('# Index of every released differential extraction (binnings of the 0.50 migration criterion).\n'
                '# One row per extraction; A_C / curves / cov files are present for each. total_pct is the\n'
                '# prediction-total fractional uncertainty (bin-averaged, excluding data statistics) of the official\n'
                '# covariance, framework + MCS term (cov_total_plusMCS.txt); chi2/p are the closure against the\n'
                '# realised pseudo-data truth. status: report/status.tsv. The one-bin totals are in total_xsec.tsv.\n')
        o.write('tag\tfamily\tconfig\tobservable\tbins\tsigma_int\ttotal_pct\tdetVar_pct\tchi2\tp\tA_C\tcurves\tcov_components\tstatus\n')
        for r in rows:
            o.write('\t'.join(r) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', dict(Counter(r[-1] for r in rows)))


if __name__ == '__main__':
    main()
