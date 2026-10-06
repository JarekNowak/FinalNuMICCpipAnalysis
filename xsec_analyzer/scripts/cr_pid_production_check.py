"""cr_pid_production_check.py -- do the per-track PID variables agree between beam-on data and the
prediction equally well in the run periods of the older ntuple production (FHC Run 1, RHC Runs 1 and
3b) and of the newer one (FHC Runs 4, 5, RHC Run 4)?

The older beam-on, beam-off and dirt ntuples lack about twenty track branches that the overlays carry
(among them trk_bragg_pion_v, on which the single-pion selection cuts). This checks whether the
productions also differ in the variables both carry. Only the opened control regions are used: the
beam-on control-region skim and, for the prediction, the same regions in the simulation (sb_pi0/),
the per-run beam-off and the per-run dirt, with the scales of bkgfit/farsb.py and bkgfit/dirt.py.

Tracks: primary (generation 2), track score >= 0.5, in events of any of the four control regions;
the muon candidate and the other tracks separately. For each variable and period: the data/prediction
shape chi2 (normalisation removed) and the mean difference.

    python3 scripts/cr_pid_production_check.py  ->  report/planning/multipion/cr_pid_production_check.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun

OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
S = 'CC1mu1piXp'
OLD = ('FHC_R1', 'RHC_R1', 'RHC_R3')          # periods whose beam-on ntuples are of the older production
VARS = {   # name: (branch, bins)
    'llr_score': ('trk_llr_pid_score_v', np.linspace(-1, 1, 21)),
    'llr_u': ('trk_llr_pid_u_v', np.linspace(-40, 40, 21)),
    'llr_v': ('trk_llr_pid_v_v', np.linspace(-40, 40, 21)),
    'llr_y': ('trk_llr_pid_y_v', np.linspace(-40, 40, 21)),
    'chi2_proton': ('trk_pid_chipr_v', np.linspace(0, 300, 21)),
    'track_score': ('trk_score_v', np.linspace(0.5, 1.0, 11)),
    'length_cm': ('trk_len_v', np.linspace(0, 200, 21)),
    'mcs_over_range': (None, np.linspace(0, 3, 16)),
}
VALID = {'llr_score': (-1.0001, 1.0001), 'llr_u': (-1e3, 1e3), 'llr_v': (-1e3, 1e3), 'llr_y': (-1e3, 1e3),
         'chi2_proton': (0., 1e4), 'track_score': (0., 1.0001), 'length_cm': (0., 1.1e3), 'mcs_over_range': (0., 100.)}
EVB = [f'{S}_sb_cc0pi', f'{S}_sb_pi0', f'{S}_sb_multipi', f'{S}_sb_cosmic', f'{S}_CandidateMuonIndex',
       'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
TRB = ['pfp_generation_v', 'trk_score_v', 'trk_len_v', 'trk_llr_pid_score_v', 'trk_llr_pid_u_v', 'trk_llr_pid_v_v',
       'trk_llr_pid_y_v', 'trk_pid_chipr_v', 'trk_mcs_muon_mom_v', 'trk_range_muon_mom_v']


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def fill(path, scale, weighted, H):
    """add the tracks of the control-region events of one file to H[(role, var)] = (sum w, sum w^2)"""
    t = uproot.open(path)['stv_tree']
    keys = set(t.keys())
    for a in t.iterate([b for b in EVB + TRB if b in keys], step_size='200 MB', library='ak'):
        cr = a[f'{S}_sb_cc0pi'] | a[f'{S}_sb_pi0'] | a[f'{S}_sb_multipi'] | a[f'{S}_sb_cosmic']
        a = a[cr]
        if len(a) == 0: continue
        w = safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) if weighted \
            else np.ones(len(a))
        w = w * scale
        idx = ak.local_index(a['trk_score_v'])
        prim = (a['pfp_generation_v'] == 2) & (a['trk_score_v'] >= 0.5)
        mu = idx == a[f'{S}_CandidateMuonIndex']
        ww = ak.broadcast_arrays(w, a['trk_score_v'])[0]
        for role, m in (('muon candidate', prim & mu), ('other tracks', prim & ~mu)):
            wf = ak.to_numpy(ak.flatten(ww[m]))
            for v, (br, bins) in VARS.items():
                if br is None:
                    mc, rg = ak.flatten(a['trk_mcs_muon_mom_v'][m]), ak.flatten(a['trk_range_muon_mom_v'][m])
                    x = ak.to_numpy(ak.where(rg > 0, mc / ak.where(rg > 0, rg, 1), -1))
                else:
                    x = ak.to_numpy(ak.flatten(a[br][m]))
                lo, hi = VALID[v]
                ok = np.isfinite(x) & (x >= lo) & (x <= hi)          # physical values; the rest are sentinels
                h1, _ = np.histogram(np.clip(x[ok], bins[0], bins[-1] - 1e-9), bins, weights=wf[ok])
                h2, _ = np.histogram(np.clip(x[ok], bins[0], bins[-1] - 1e-9), bins, weights=wf[ok] ** 2)
                cur = H.get((role, v)) or [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1), 0., 0., 0.]
                H[(role, v)] = [cur[0] + h1, cur[1] + h2, cur[2] + float((x[ok] * wf[ok]).sum()), cur[3] + float(wf[ok].sum()),
                                cur[4] + float(wf[~ok].sum())]


def shape_chi2(d, d2, p, p2):
    D, P = d.sum(), p.sum()
    if D <= 0 or P <= 0: return float('nan'), 0
    r = D / P
    var = d2 + p2 * r * r
    m = var > 0
    return float(((d - p * r)[m] ** 2 / var[m]).sum()), int(m.sum() - 1)


def main():
    mc, ext, data = F.samples()
    dirts = dirt_perrun('comb')            # processed/ per-run dirt (control-region flags identical in every copy)
    out = {}
    for p, name in enumerate(F.PNAME):
        if not data[p]: continue
        Hd, Hp = {}, {}
        for f in data[p]: fill(f, 1.0, False, Hd)
        for f, s in mc[p]: fill(f, s, True, Hp)
        for f, s in ext[p]: fill(f, s, False, Hp)
        for r, f, s in dirts:
            if r == F.RUN[p]: fill(f, s, True, Hp)
        res = {}
        for key in Hd:
            d1, d2, dsx, dn, dbad = Hd[key]; p1, p2, psx, pn, pbad = Hp[key]
            chi2, ndf = shape_chi2(d1, d2, p1, p2)
            res[f'{key[0]} | {key[1]}'] = dict(chi2=chi2, ndf=ndf, data_n=float(dn), pred_n=float(pn),
                                               data_mean=dsx / dn if dn else float('nan'),
                                               pred_mean=psx / pn if pn else float('nan'),
                                               data_invalid_frac=dbad / (dn + dbad) if dn + dbad else float('nan'),
                                               pred_invalid_frac=pbad / (pn + pbad) if pn + pbad else float('nan'))
        out[name] = res
        print(f'== {name} done', flush=True)
    json.dump(out, open(os.path.join(OUTDIR, 'cr_pid_production_check.json'), 'w'), indent=1)
    L = ['# PID variables in the opened control regions: beam-on data against the prediction, per period', '',
         'Older-production periods (beam-on ntuples without the extra track branches): ' + ', '.join(OLD) + '.',
         'Shape chi2/ndf of the physical values (normalisation removed; data and prediction statistics), mean data - prediction,',
         'and the fraction of sentinel (non-physical) values in data / prediction.', '']
    roles_vars = sorted({k for r in out.values() for k in r})
    periods = list(out)
    L += ['| track | variable | ' + ' | '.join(f'{p}{" (old)" if p in OLD else ""}' for p in periods) + ' |',
          '|---|---|' + '---|' * len(periods)]
    for k in roles_vars:
        role, var = k.split(' | ')
        cells = []
        for p in periods:
            r = out[p].get(k)
            cells.append(f'{r["chi2"]:.1f}/{r["ndf"]} ({r["data_mean"] - r["pred_mean"]:+.3g}; inv {100 * r["data_invalid_frac"]:.0f}/{100 * r["pred_invalid_frac"]:.0f}%)' if r else 'n/a')
        L.append(f'| {role} | {var} | ' + ' | '.join(cells) + ' |')
    open(os.path.join(OUTDIR, 'cr_pid_production_check.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
