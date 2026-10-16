"""mcs_cal_holdout.py -- hold-out test of the MCS calibration of exiting muon candidates (scripts/mcs_cal_derive.py):
the nominal quantile maps derived on part of the control-region candidates and tested on the rest, which the
derivation did not see.

  1. Random halves: data and prediction candidates are split at random (fixed seed) into halves A and B; the maps
     derived on A are applied to the simulation of B and compared with the data of B, and the other way round.
  2. Groups of periods: maps derived on FHC are tested on RHC and the other way round, and likewise old production
     (FHC Run 1, RHC Runs 1 and 3) against new production (FHC Runs 4 and 5, RHC Run 4).

Each test sample is shown without calibration, with the maps derived on all candidates (the release, in-sample for
the test events) and with the hold-out maps. Closure quantities as in mcs_cal_derive.py: the released p_mu bins of
all candidates (data / prediction, shape, with their chi2), the fraction of exiting candidates with an MCS momentum
above 1.2 GeV/c per |dir x| and length bin (data / prediction), and the chi2 of the MCS spectrum of exiting
candidates (shape).

    python3 scripts/mcs_cal_holdout.py  ->  report/results/mcs_cal_holdout.{md,json}
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mcs_cal_derive as M

OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'results'))
OLD, NEW = ['FHC_R1', 'RHC_R1', 'RHC_R3'], ['FHC_R4', 'FHC_R5', 'RHC_R4']
FHC, RHC = ['FHC_R1', 'FHC_R4', 'FHC_R5'], ['RHC_R1', 'RHC_R3', 'RHC_R4']


def sub(X, m):
    return {k: v[m] for k, v in X.items()}


def pmu_chi2(D, P, pm):
    """shape chi2 of the released p_mu bins (data and simulation statistics)"""
    hd = np.histogram(np.clip(D['reco'][D['ok']], 0, 49.9), M.EB)[0]
    hp = np.histogram(np.clip(pm[P['ok']], 0, 49.9), M.EB, weights=P['w'][P['ok']])[0]
    hp2 = np.histogram(np.clip(pm[P['ok']], 0, 49.9), M.EB, weights=P['w'][P['ok']] ** 2)[0]
    n = hd.sum() / hp.sum(); var = hd + n * n * hp2; m = var > 0
    return float(((hd - n * hp)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)


def evaluate(D, P, sel_d, sel_p, maps_all, maps_hold):
    """closure of the test sample (sel_d, sel_p) without calibration, with the release maps, with the hold-out maps"""
    Dt, Pt = sub(D, sel_d), sub(P, sel_p)
    res = {}
    for lab, maps in (('none', None), ('release (all candidates)', maps_all), ('hold-out', maps_hold)):
        pm = Pt['reco'] if maps is None else M.apply(Pt, M.NOMINAL, maps)
        c = M.closure(Dt, Pt, pm)
        c['chi2_pmu'] = list(pmu_chi2(Dt, Pt, pm))
        res[lab] = c
    res['test_data_candidates'] = int(Dt['ok'].sum()); res['test_data_exiting'] = int(Dt['exit'].sum())
    return res


def toys(D, P, maps_all, n=10, seed=20261017):
    """Statistical expectation of the random-half test: pseudo-data drawn (with replacement, by weight) from one
    third of the prediction with the release calibration applied, so that the true map is known; the maps are derived
    on half of the pseudo-data with a second third of the prediction and tested on the other half with the last third.
    The chi2 of the test then measures the statistics of the derivation and of the test alone (the thirds of the
    prediction carry less simulation statistics than the real test, so the expectation is on the high side)."""
    rng = np.random.default_rng(seed)
    cal = M.apply(P, M.NOMINAL, maps_all)
    res = []
    for _ in range(n):
        part = rng.integers(0, 3, len(P['w']))
        src = np.flatnonzero(part == 2)
        pr = P['w'][src] / P['w'][src].sum()
        pick = src[rng.choice(len(src), size=len(D['w']), p=pr)]
        T = {k: P[k][pick].copy() for k in ('mcs', 'adx', 'dz', 'len', 'sx', 'ok', 'exit', 'period')}
        T['reco'] = cal[pick]; T['mcs'] = cal[pick]; T['w'] = np.ones(len(pick)); T['sim'] = np.zeros(len(pick), bool)
        half = rng.random(len(pick)) < 0.5
        maps_h, _ = M.derive(T, P, 'quantile', M.NOMINAL, half, part == 0)
        Tt, Pt = sub(T, ~half), sub(P, part == 1)
        pm = M.apply(Pt, M.NOMINAL, maps_h)
        res.append((pmu_chi2(Tt, Pt, pm)[0], M.closure(Tt, Pt, pm)['chi2_mcs'][0]))
    r = np.array(res)
    return dict(n=n, chi2_pmu=[float(r[:, 0].mean()), float(r[:, 0].std())], chi2_mcs=[float(r[:, 1].mean()), float(r[:, 1].std())])


def main():
    D, P = M.candidates()
    alld, allp = np.ones(len(D['w']), bool), np.ones(len(P['w']), bool)
    maps_all, _ = M.derive(D, P, 'quantile', M.NOMINAL, alld, allp)
    rng = np.random.default_rng(20261016)
    ad, ap = rng.random(len(D['w'])) < 0.5, rng.random(len(P['w'])) < 0.5
    pid = lambda X, names: np.isin(X['period'], [M.F.PNAME.index(n) for n in names])
    tests = [('random half B, maps from half A', (ad, ap), (~ad, ~ap)),
             ('random half A, maps from half B', (~ad, ~ap), (ad, ap)),
             ('RHC, maps from FHC', (pid(D, FHC), pid(P, FHC)), (pid(D, RHC), pid(P, RHC))),
             ('FHC, maps from RHC', (pid(D, RHC), pid(P, RHC)), (pid(D, FHC), pid(P, FHC))),
             ('new production, maps from old', (pid(D, OLD), pid(P, OLD)), (pid(D, NEW), pid(P, NEW))),
             ('old production, maps from new', (pid(D, NEW), pid(P, NEW)), (pid(D, OLD), pid(P, OLD)))]
    out = {}
    fmt = lambda v: ' '.join(f'{x:.2f}' for x in v)
    L = ['# Hold-out test of the MCS calibration of exiting muon candidates', '',
         'Generated by `scripts/mcs_cal_holdout.py` from the control-region candidates of `scripts/mcs_cal_derive.py` '
         '(beam-on data against overlay, beam-off and dirt, all opened control regions). The nominal quantile maps are '
         'derived on one part of the candidates and tested on the other; "release" applies the maps derived on all '
         'candidates (in-sample for the test events). p_mu bins: data / prediction in the p_mu bins of version 1.8 (the '
         'two highest are one bin from version 1.9), shape; '
         '|dir x| and length: data / prediction of the fraction of exiting candidates with an MCS momentum above '
         '1.2 GeV/c; chi2: shape chi2 of the released p_mu bins (5 dof) and of the MCS spectrum of exiting candidates.', '']
    for name, (der_d, der_p), (tst_d, tst_p) in tests:
        maps_hold, _ = M.derive(D, P, 'quantile', M.NOMINAL, der_d, der_p)
        r = evaluate(D, P, tst_d, tst_p, maps_all, maps_hold)
        out[name] = r
        L += [f'## Test: {name}', '', f'Test sample: {r["test_data_candidates"]} data candidates, '
              f'{r["test_data_exiting"]} exiting.', '',
              '| calibration | p_mu bins | chi2 p_mu | abs(dir x) | length | chi2 MCS |', '|---|---|---|---|---|---|']
        for lab in ('none', 'release (all candidates)', 'hold-out'):
            c = r[lab]
            L.append(f'| {lab} | {fmt(c["pmu"])} | {c["chi2_pmu"][0]:.1f} / {c["chi2_pmu"][1]} | {fmt(c["abs(dir x)"])} | '
                     f'{fmt(c["length [cm]"])} | {c["chi2_mcs"][0]:.0f} / {c["chi2_mcs"][1]} |')
        L.append('')
        print('==', name, flush=True)
    t = toys(D, P, maps_all)
    out['toys'] = t
    L += [f'## Statistical expectation of the random-half test ({t["n"]} toys)', '',
          'Pseudo-data from the calibrated simulation, so that the true map is the release map; the half-split test run '
          'on them measures the statistics of the derivation and of the test alone (with thirds of the prediction, a '
          f'high-side expectation): chi2 p_mu {t["chi2_pmu"][0]:.1f} +- {t["chi2_pmu"][1]:.1f} (5 dof), chi2 MCS '
          f'{t["chi2_mcs"][0]:.0f} +- {t["chi2_mcs"][1]:.0f} (24 dof).', '']
    L += ['Bins: p_mu ' + ', '.join(f'{M.EB[b]:g}-{M.EB[b + 1]:g}' for b in range(len(M.EB) - 1)) + ' GeV/c; abs(dir x) '
          + ', '.join(f'{x:g}' for x in M.SPLIT['abs(dir x)'][1]) + '; length ' + ', '.join(f'{x:g}' for x in M.SPLIT['length [cm]'][1]) + ' cm.']
    json.dump(out, open(os.path.join(OUTDIR, 'mcs_cal_holdout.json'), 'w'), indent=1, default=float)
    open(os.path.join(OUTDIR, 'mcs_cal_holdout.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
