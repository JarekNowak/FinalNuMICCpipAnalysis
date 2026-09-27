"""official_cov.py -- the official covariance of every extraction and the numbers derived from it.

Decision of 2026-09-27: the official total covariance of a result is the framework covariance plus the
data-side MCS muon-momentum-scale term, i.e. data_release/cov/<family>_<CFG>_<obs>/cov_total_plusMCS.txt
(cov_total.txt keeps the framework covariance alone; cov_MCSscale.txt is the term, written by
report/tools/mcs_eval.py). Every summary table takes its uncertainties from here, so the tables and the
release cannot disagree about whether the term is included.

The bin-averaged fractional uncertainties of the systematics dumps (logs/systdump/*.dump) are
reproduced from the release covariances to better than 0.03 percentage points. Where the MCS term is
non-zero, the official value is the dump value plus the exact increment computed from the release
(framework+MCS minus framework, both bin-averaged over sqrt(C_ii)/x_i), so extractions without the
term keep their dump values unchanged.

    import official_cov as oc
    d = oc.summary('incl', 'FHC5', 'pmu')   # dump keys, with PredTotal/total official and 'MCS' added
    oc.chi2('incl', 'FHC5', 'pmu', 'GENIE_smeared')   # model chi2 with the official covariance
"""
import os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REP = os.path.dirname(HERE)
REL = os.path.join(REP, 'data_release')
DUMP = os.path.join(os.path.dirname(REP), 'logs', 'systdump')
CF = {'FHC5': 'fhc5', 'RHCFULL': 'rhcfull', 'COMB': 'comb'}


def cov(path):
    C = None
    for line in open(path):
        s = line.split()
        if s[0] == 'numXbins':
            n = int(s[1]); C = np.zeros((n, n))
        elif s[0] not in ('numYbins', 'xbin') and len(s) == 3:
            C[int(s[0]), int(s[1])] = float(s[2])
    return C


def vec(path):
    return np.array([float(l.split()[1]) for l in open(path) if not l.startswith('numXbins') and l.strip()])


def reldir(fam, cfg, obs):
    return os.path.join(REL, 'cov', f'{fam}_{cfg}_{obs}')


def dumpfile(fam, cfg, obs):
    if fam == 'incl':
        return os.path.join(DUMP, f"{CF[cfg]}_{'ppi' if obs == 'ppi3bin' else obs}.dump")
    return os.path.join(DUMP, f'1p_{CF[cfg]}_{obs}.dump')


def load_dump(fam, cfg, obs):
    d = {}
    for line in open(dumpfile(fam, cfg, obs)):
        if line.startswith('[SYSTDUMP]'):
            _, k, v = line.split(); d[k] = float(v)
    return d


def _binavg(C, x):
    return 100 * float(np.mean(np.sqrt(np.clip(np.diag(C), 0, None)) / np.abs(x)))


def official_cov_path(fam, cfg, obs):
    d = reldir(fam, cfg, obs)
    p = os.path.join(d, 'cov_total_plusMCS.txt')
    return p if os.path.exists(p) else os.path.join(d, 'cov_total.txt')


def summary(fam, cfg, obs):
    """Dump values with PredTotal and total replaced by their official (framework + MCS) values.
    Adds 'MCS' (bin-averaged %, 0 when the term is absent or zero) and keeps the framework values as
    'PredTotal_fw' and 'total_fw'."""
    d = dict(load_dump(fam, cfg, obs))
    d['PredTotal_fw'], d['total_fw'] = d['PredTotal'], d['total']
    rd = reldir(fam, cfg, obs); m = os.path.join(rd, 'cov_MCSscale.txt')
    d['MCS'] = 0.0
    if os.path.exists(m):
        M = cov(m)
        if np.any(M != 0):
            x = vec(os.path.join(rd, 'unfolded_signal.txt'))
            T = cov(os.path.join(rd, 'cov_total.txt')); P = cov(os.path.join(rd, 'cov_PredTotal.txt'))
            d['MCS'] = _binavg(M, x)
            d['PredTotal'] = d['PredTotal_fw'] + _binavg(P + M, x) - _binavg(P, x)
            d['total'] = d['total_fw'] + _binavg(T + M, x) - _binavg(T, x)
    return d


def curves(fam, cfg, obs):
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(REL, f'curves_{fam}_{cfg}_{obs}.tsv')) if not l.startswith('#')]
    h = rows[0]
    return {k: np.array([float(r[i]) for r in rows[1:]]) for i, k in enumerate(h)}


def chi2(fam, cfg, obs, column, official=True):
    """chi2 of the unfolded data against an A_C-smeared model column of curves_*.tsv, with the official
    covariance (or the framework covariance with official=False). With the framework covariance this
    reproduces the unfolder's model chi2 exactly; the unfolder's 'truth' chi2 also contains the
    covariance of the fake-data truth, which is not released, and is not reproduced here."""
    c = curves(fam, cfg, obs)
    C = cov(official_cov_path(fam, cfg, obs) if official else os.path.join(reldir(fam, cfg, obs), 'cov_total.txt'))
    dlt = (c['unfolded_data'] - c[column]) * c['width']
    return float(dlt @ np.linalg.solve(C, dlt)), len(dlt)
