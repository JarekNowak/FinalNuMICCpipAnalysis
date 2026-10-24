"""pelee_compare_rates.py -- true CC rates per 1e20 POT in the overlays of report/tools/pelee_compare.py (ours, custom,
reweighted custom): nu_mu, nu_mu-bar, nu_e CC with the true vertex in the fiducial volume, unweighted, tune-weighted and
tune x ppfx_cv weighted, and the mean true energy.  ->  report/planning/pelee_compare_rates.json
    nice python3 report/tools/pelee_compare_rates.py
"""
import sys, json, numpy as np, uproot
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__))); import pelee_compare as pc
from multiprocessing import Pool
B = ['nu_pdg', 'ccnc', 'true_nu_vtx_x', 'true_nu_vtx_y', 'true_nu_vtx_z', 'ppfx_cv', 'weightSplineTimesTune', 'nu_e']
def rate(paths):
    pot = 0.; n = {}
    for p in paths:
        f = uproot.open(p); pot += float(f['nuselection/SubRun']['pot'].array(library='np').sum())
        a = f['nuselection/NeutrinoSelectionFilter'].arrays(B, library='np')
        fv = (a['true_nu_vtx_x'] > 10) & (a['true_nu_vtx_x'] < 246) & (np.abs(a['true_nu_vtx_y']) < 106) & (a['true_nu_vtx_z'] > 10) & (a['true_nu_vtx_z'] < 1026)
        t = a['weightSplineTimesTune']; t = np.where(np.isfinite(t) & (t > 0) & (t < 30), t, 1.)
        for lab, s in (('numu', a['nu_pdg'] == 14), ('numubar', a['nu_pdg'] == -14), ('nue', np.abs(a['nu_pdg']) == 12)):
            m = s & (a['ccnc'] == 0) & fv
            r = n.setdefault(lab, [0, 0., 0., 0.])
            r[0] += int(m.sum()); r[1] += float(t[m].sum()); r[2] += float((t * a['ppfx_cv'])[m].sum()); r[3] += float(a['nu_e'][m].sum())
    return dict(pot=pot, **{k: dict(n=v[0], per1e20=v[0] / pot * 1e20, tune=v[1] / pot * 1e20, tune_ppfx=v[2] / pot * 1e20, mean_E=v[3] / max(v[0], 1)) for k, v in n.items()})
def job(x):
    lab, side, paths = x; return lab, side, rate(paths)
if __name__ == '__main__':
    J = []
    for lab, o, c in pc.PAIRS[:8]: J += [(lab, 'ours', o), (lab, 'custom', c)]
    for lab, o, c in pc.PAIRS_RW: J.append((lab.replace(' (rw)', ''), 'custom_rw', c))
    with Pool(3) as pool: R = pool.map(job, J, chunksize=1)
    out = {}
    for lab, side, r in R: out.setdefault(lab, {})[side] = r
    json.dump(out, open('report/planning/pelee_compare_rates.json', 'w'), indent=1)
    for lab, d in out.items():
        for side, r in d.items():
            print(f"{lab:20s} {side:10s} pot {r['pot']:.3e} " + ' '.join(f"{k} n {r[k]['per1e20']:8.1f} tune*ppfx {r[k]['tune_ppfx']:8.1f} <E> {r[k]['mean_E']:.3f}" for k in ('numu', 'numubar')))
