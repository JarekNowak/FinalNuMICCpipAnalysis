"""mp2pi_llr_recomb.py -- size of the recombination variation that the LLR mismatch of the control regions calls for,
and its effect on the two-pion results (two-fold selection).

1. Strength alpha of the Recomb2 shape in the opened control regions: LLR score of the classifier's pool tracks
   (the 'LLR score' test of scripts/cr_mp2pi_syst_check.py --prod 2f, its cache syst_check_hist.pkl), prediction
   p + alpha * (Recomb2 shape deviation), shape chi2 with the statistics and every other term of that test (flux,
   cross section, reinteraction, the seven other detector variations); linear in alpha. Cross-check: the LLR of the
   long tracks alone (pool tracks >= 40 cm, cache of scripts/mp2pi_figures_cr.py), statistics only.
2. Effect: the framework carries Recomb2 as one detector variation applied fractionally, so a variation alpha times as
   large moves each result by alpha times its Recomb2 term (cov_detVarRecomb2 of the extraction). Quoted against the
   detector term and the total uncertainty, with the constraint: the one-bin total per configuration and, combined,
   the integral and the largest bin of each differential variable. Selected yields: Recomb2 over its CV per region
   (Run-4 detector-variation samples in the universe file of the one-bin total).

    python3 scripts/mp2pi_llr_recomb.py  ->  report/planning/multipion/mp2pi_llr_recomb.{md,json}
"""
import json, os, pickle, re, sys
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
if '--prod' not in sys.argv: sys.argv += ['--prod', '2f']
import cr_mp2pi_syst_check as C

PLAN = C.OUTDIR
U = '/data/uboone/processed/mp2pi_2f/unfold/'
VARS = ['total', 'pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'np']
CFG = {'fhc5': 'FHC', 'rhcfull': 'RHC', 'comb': 'combined'}


def fit_alpha(d, p, Cother, dv):
    """chi2(alpha) = x^T V^-1 x, x = d - r (p + alpha dv), V = r^2 Cother + diag(d), last bin dropped (as C.test)"""
    r = d.sum() / p.sum(); V = r * r * Cother + np.diag(d)
    k = np.flatnonzero((np.diag(V) > 0) | (d > 0))[:-1]
    Vi = np.linalg.pinv(V[np.ix_(k, k)]); x0 = (d - r * p)[k]; g = (r * dv)[k]
    a = float(g @ Vi @ x0 / (g @ Vi @ g)); s = float(1 / np.sqrt(g @ Vi @ g))
    chi = lambda al: float((x0 - al * g) @ Vi @ (x0 - al * g))
    return dict(alpha=a, err=s, chi2_fit=chi(a), chi2_0=chi(0.), chi2_1=chi(1.), ndf=len(k) - 1)


def control_region():
    HP, R = pickle.load(open(C.SK + 'syst_check_hist.pkl', 'rb'))
    q = 'LLR score'
    d = sum(HP[h][('data', q, 'cv')] for h in HP)
    p_h = {h: HP[h][('mc', q, 'cv')] + HP[h].get(('ext', q, 'cv'), 0.) + HP[h].get(('dirt', q, 'cv'), 0.) for h in HP}
    p = sum(p_h.values())
    Cs = np.diag(sum(HP[h][('mc', q, 'w2')] + HP[h].get(('ext', q, 'w2'), 0.) + HP[h].get(('dirt', q, 'w2'), 0.) for h in HP))
    Co = Cs.copy()
    for term in C.UNIV:
        Pu = sum(HP[h][('mc', q, term)] + (p_h[h] - HP[h][('mc', q, 'cv')])[:, None] for h in HP)
        dev = np.stack([C.shape_dev(Pu[:, u], p) for u in range(Pu.shape[1])], axis=1); Co += dev @ dev.T / dev.shape[1]
    dvs = {}
    for v in C.KNOBS:
        pv = sum(p_h[h] + HP[h][('mc', q, 'cv')] * (R[(h.lower(), v)][q] - 1.) for h in HP); dvs[v] = C.shape_dev(pv, p)
    for v in C.KNOBS:
        if v != 'Recomb2': Co += np.outer(dvs[v], dvs[v])
    res = fit_alpha(d, p, Co, dvs['Recomb2'])
    # cross-check: long tracks, statistics only (cache of scripts/mp2pi_figures_cr.py)
    H = pickle.load(open('/data/uboone/processed/mp2pi_cr_skim_mcscal/fig_llr_hist.pkl', 'rb'))
    dl = sum(H[h]['data'] for h in H); pl = sum(H[h]['overlay'] + H[h]['beam-off'] + H[h]['dirt'] for h in H)
    pr = sum(H[h]['overlay'] * H[h]['ratio_Recomb2'] + H[h]['beam-off'] + H[h]['dirt'] for h in H)
    long_ = fit_alpha(dl, pl, np.diag(sum(H[h]['w2'] for h in H)), C.shape_dev(pr, pl))
    return res, long_


def cov(path):
    rows = [l.split() for l in open(path) if l.strip() and l.split()[0][0].isdigit()]
    if len(rows[0]) == 2: return np.array([float(r[1]) for r in rows])
    n = max(int(r[0]) for r in rows) + 1; m = np.zeros((n, n))
    for r in rows: m[int(r[0]), int(r[1])] = float(r[2])
    return m


def results(alpha):
    out = {}
    for v in VARS:
        for cfg in (CFG if v == 'total' else ['comb']):
            dd = f'{U}{v}/{cfg}/unfold_output/'       # written by the constrained run, the last of the two
            if not os.path.exists(dd + 'mat_table_cov_detVarRecomb2.txt'): continue
            x = cov(dd + 'vec_table_unfolded_signal.txt')
            M = {k: cov(dd + f'mat_table_cov_{k}.txt') for k in ('detVarRecomb2', 'detVar_total', 'total')}
            integ = {k: np.sqrt(m.sum()) / x.sum() for k, m in M.items()}
            perbin = {k: np.sqrt(np.diag(m)) / np.abs(x) for k, m in M.items()}
            ratio_bins = alpha * perbin['detVarRecomb2'] / perbin['total']
            out[f'{v} {cfg}'] = dict(recomb_int=float(integ['detVarRecomb2']), det_int=float(integ['detVar_total']),
                                     total_int=float(integ['total']), shift_over_det=float(alpha * integ['detVarRecomb2'] / integ['detVar_total']),
                                     shift_over_total=float(alpha * integ['detVarRecomb2'] / integ['total']),
                                     max_bin_shift_over_total=float(ratio_bins.max()), nbins=len(x))
    return out


def yields():
    f = uproot.open('/data/uboone/processed/mp2pi_COMB_total_2f_fluxw_total_univmake.root')
    d = f['mp2pi_CC1mu2piBDT_total_1D']
    get = lambda horn, v: d[f'+data+uboone+processed+mp2pi_2f+detvar_fluxw+xsec-ana-detvar_run4{horn}_{v}.root/unweighted_0_reco'].values()
    return {horn.upper(): (get(horn, 'Recomb2') / get(horn, 'CV')).tolist() for horn in ('fhc', 'rhc')}


def main():
    cr, lg = control_region()
    a = cr['alpha']; res = results(a); y = yields()
    out = dict(control_region=cr, long_tracks=lg, results=res, yield_ratio_recomb2=y)
    json.dump(out, open(os.path.join(PLAN, 'mp2pi_llr_recomb.json'), 'w'), indent=1)
    L = ['# Recombination variation called for by the LLR of the control regions, and its effect on the two-pion results', '',
         'Generated by `scripts/mp2pi_llr_recomb.py` (see its docstring). alpha: strength of the Recomb2 shape (1 = the '
         'detector variation as simulated).', '',
         f"Control regions, all pool tracks, statistics and every other term: alpha = {cr['alpha']:.2f} +- {cr['err']:.2f}; "
         f"chi2 {cr['chi2_0']:.1f} (alpha 0), {cr['chi2_1']:.1f} (alpha 1), {cr['chi2_fit']:.1f} (fit) for {cr['ndf']} dof.",
         f"Long tracks (>= 40 cm), statistics only: alpha = {lg['alpha']:.2f} +- {lg['err']:.2f}; chi2 {lg['chi2_0']:.1f} / "
         f"{lg['chi2_1']:.1f} / {lg['chi2_fit']:.1f} for {lg['ndf']} dof.", '',
         'Selected yields, Recomb2 over its CV (Run 4) per region (signal region, three-pion-like, proton-like, score band, low score):', '']
    for h, r in y.items(): L.append(f'- {h}: ' + ', '.join(f'{x:.3f}' for x in r))
    L += ['', f'Effect with alpha = {a:.2f}, with the constraint (integral over the bins; largest bin against its total):', '',
          '| result | Recomb2 term | detector term | total | shift / detector | shift / total | largest bin shift / total |', '|---|---|---|---|---|---|---|']
    for k, r in res.items():
        L.append(f"| {k} | {100 * r['recomb_int']:.1f}% | {100 * r['det_int']:.1f}% | {100 * r['total_int']:.1f}% | "
                 f"{r['shift_over_det']:.2f} | {r['shift_over_total']:.2f} | {r['max_bin_shift_over_total']:.2f} |")
    open(os.path.join(PLAN, 'mp2pi_llr_recomb.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
