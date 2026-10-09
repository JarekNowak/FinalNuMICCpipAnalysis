"""cr_mcs_check.py -- MCS momentum of the muon candidate, beam-on data against the prediction in the opened control
regions of the single-pion analysis, per period. Question: is the deficit of high MCS momenta of long exiting tracks
in data (scripts/cr_mp2pi_check.py) a momentum-scale bias of the MCS estimator or a property of the sample?

(1) Contained muon candidates: the ratio MCS / range momentum per range-momentum bin, mean and median, data and
    prediction. A scale bias of the MCS estimator in data shows here directly, at the momenta of contained muons.
(2) Exiting muon candidates: the MCS momentum spectrum, data over prediction with the normalisation removed.
(3) The same exiting spectrum against the true momentum in the simulation: the MCS bias of the simulation itself.
(4) Exiting muon candidates in bins of their visible length (from the CC1mu2piBDTEvt track dump): the length
    distribution and the MCS spectrum per length bin, data over prediction with the normalisation removed. If the
    data tracks were shorter (broken tracks), the MCS precision would differ at fixed momentum.
(5) The reconstructed muon momentum (range when contained, MCS when exiting) in the bins of the released p_mu result,
    all, contained and exiting candidates, data over prediction with the normalisation removed.

Input: the control-region skim /data/uboone/processed/mp2pi_cr_skim (slurm/slurm_mp2pi_cr_skim.sbatch). Scales as
scripts/cr_mp2pi_check.py.
    python3 scripts/cr_mcs_check.py  ->  report/planning/multipion/cr_mcs_check.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun
import cr_mp2pi_check as C

S = 'CC1mu1piXp'
RB = np.array([0.15, 0.25, 0.35, 0.5, 0.7, 1.0])            # range-momentum bins of contained candidates
QB = np.linspace(0.0, 2.5, 51)                              # MCS / range
PB = np.array([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 50.])   # MCS momentum of exiting candidates
LB = np.array([0., 50., 100., 150., 200., 300., 1000.])          # visible length of exiting candidates [cm]
MB = np.array([0., 0.4, 0.8, 1.2, 1.6, 50.])                     # coarse MCS bins per length bin
B = 'CC1mu2piBDTEvt'
EB = np.array([0.15, 0.3, 0.425, 0.575, 0.825, 1.25, 50.])       # released p_mu bins (ccpi_pmu_slice_config_opt)
GROUPS = {'FHC Run 1': ['FHC_R1'], 'FHC Runs 4-5': ['FHC_R4', 'FHC_R5'], 'RHC Runs 1, 3': ['RHC_R1', 'RHC_R3'],
          'RHC Run 4': ['RHC_R4'], 'all': None}


def fill(path, scale, weighted, H, truth):
    t = uproot.open(path)['stv_tree']
    br = [f'{S}_sb_{c}' for c in ('cc0pi', 'multipi', 'pi0', 'cosmic')] + \
         [f'{S}_candidate_muon_mom_mcs', f'{S}_candidate_muon_mom_range', f'{S}_CandidateMuonTrackEndContainment',
          f'{S}_CandidateMuonIndex', f'{S}_candidate_muon_mom_reco', f'{B}_etrk_idx', f'{B}_etrk_len'] + ([f'{S}_candidate_muon_mom_true'] if truth else []) + \
         (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in t.iterate(br, step_size='200 MB', library='np'):
        w = C.safe(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']) if weighted else np.ones(len(a[br[0]]))
        mi_ = a[f'{S}_CandidateMuonIndex']
        ln = np.array([float(L[np.flatnonzero(I == m)[0]]) if m >= 0 and np.any(I == m) else -1.
                       for I, L, m in zip(a[f'{B}_etrk_idx'], a[f'{B}_etrk_len'], mi_)])
        w = w * scale
        ok = a[f'{S}_CandidateMuonIndex'] >= 0
        mcs, rng = a[f'{S}_candidate_muon_mom_mcs'], a[f'{S}_candidate_muon_mom_range']
        con = a[f'{S}_CandidateMuonTrackEndContainment'].astype(bool)
        c = ok & con & (rng > RB[0]) & (rng < RB[-1]) & np.isfinite(mcs) & (mcs > 0)
        ib = np.digitize(rng[c], RB) - 1
        q = np.clip(mcs[c] / rng[c], QB[0], QB[-1] - 1e-9)
        for b in range(len(RB) - 1):
            m = ib == b
            H[('q', b)] = H.get(('q', b), 0) + np.histogram(q[m], QB, weights=w[c][m])[0]
        e = ok & ~con & np.isfinite(mcs) & (mcs > 0)
        reco = a[f'{S}_candidate_muon_mom_reco']
        for tag, m in (('all', ok), ('contained', ok & con), ('exiting', ok & ~con)):
            H[('pmu', tag)] = H.get(('pmu', tag), 0) + np.histogram(np.clip(reco[m], 0, 49.9), EB, weights=w[m])[0]
        el = e & (ln > 0)
        H['len'] = H.get('len', 0) + np.histogram(np.clip(ln[el], 0, 999), LB, weights=w[el])[0]
        for b in range(len(LB) - 1):
            m = el & (ln >= LB[b]) & (ln < LB[b + 1])
            H[('lmcs', b)] = H.get(('lmcs', b), 0) + np.histogram(np.clip(mcs[m], 0, 49.9), MB, weights=w[m])[0]
        H['exit'] = H.get('exit', 0) + np.histogram(np.clip(mcs[e], 0, 49.9), PB, weights=w[e])[0]
        H['exit_w2'] = H.get('exit_w2', 0) + np.histogram(np.clip(mcs[e], 0, 49.9), PB, weights=w[e] ** 2)[0]
        if truth:
            tr = a[f'{S}_candidate_muon_mom_true']
            g = e & (tr > 0)
            H['exit_true'] = H.get('exit_true', 0) + np.histogram(np.clip(tr[g], 0, 49.9), PB, weights=w[g])[0]
            H['exit_mcs_t'] = H.get('exit_mcs_t', 0) + np.histogram(np.clip(mcs[g], 0, 49.9), PB, weights=w[g])[0]
            for b in range(len(PB) - 1):
                m = g & (tr >= PB[b]) & (tr < PB[b + 1])
                H[('bias', b)] = H.get(('bias', b), 0) + np.histogram(np.clip(mcs[m] / tr[m], QB[0], QB[-1] - 1e-9), QB, weights=w[m])[0]


def stats(h):
    if not np.sum(h) > 0: return float('nan'), float('nan')
    x = 0.5 * (QB[1:] + QB[:-1]); c = np.cumsum(h) / np.sum(h)
    return float((h * x).sum() / h.sum()), float(np.interp(0.5, c, QB[1:]))


def main():
    mc, ext, data = F.samples()
    re = lambda p: C.NEW + os.path.basename(p)
    dirts = dirt_perrun('comb', where=lambda run: C.NEW)
    per = {}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        Hd, Hp, Ht = {}, {}, {}
        for f in data[ip]: fill(re(f), 1.0, False, Hd, False)
        for f, s in mc[ip]: fill(re(f), s, True, Hp, True); fill(re(f), s, True, Ht, True)
        for f, s in ext[ip]: fill(re(f), s, False, Hp, False)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, True, Hp, True)
        per[name] = (Hd, Hp, Ht)
        print('==', name, flush=True)
    out = {}
    for g, names in GROUPS.items():
        use = [n for n in per if names is None or n in names]
        if not use: continue
        add = lambda i, k: sum(per[n][i].get(k, 0) for n in use)
        res = {'contained': [], 'exiting': {}, 'mc_bias': []}
        for b in range(len(RB) - 1):
            md, mdd = stats(add(0, ('q', b))); mp, mpp = stats(add(1, ('q', b)))
            res['contained'].append(dict(range_bin=[float(RB[b]), float(RB[b + 1])], n_data=float(np.sum(add(0, ('q', b)))),
                                         mean_data=md, median_data=mdd, mean_pred=mp, median_pred=mpp))
        d, p = add(0, 'exit'), add(1, 'exit'); d2, p2 = add(0, 'exit_w2'), add(1, 'exit_w2')
        n = d.sum() / p.sum()
        res['exiting'] = dict(bins=list(map(float, PB)), data=list(map(float, d)), pred=list(map(float, p)), norm=float(n),
                              shape=list(map(float, (d / np.where(p > 0, p, np.nan)) / n)),
                              shape_err=list(map(float, np.sqrt(d2 / np.where(d > 0, d, np.nan) ** 2 + p2 / np.where(p > 0, p, np.nan) ** 2)
                                                 * (d / np.where(p > 0, p, np.nan)) / n)))
        res['pmu'] = {}
        for tag in ('all', 'contained', 'exiting'):
            d_, p_ = add(0, ('pmu', tag)), add(1, ('pmu', tag))
            res['pmu'][tag] = dict(norm=float(d_.sum() / p_.sum()), shape=list(map(float, (d_ / p_) / (d_.sum() / p_.sum()))),
                                   data=list(map(float, d_)))
        dl, pl = add(0, 'len'), add(1, 'len'); nl = dl.sum() / pl.sum()
        res['length'] = dict(bins=list(map(float, LB)), shape=list(map(float, (dl / pl) / nl)), data=list(map(float, dl)),
                             mcs_bins=list(map(float, MB)), mcs_shape=[])
        for b in range(len(LB) - 1):
            d_, p_ = add(0, ('lmcs', b)), add(1, ('lmcs', b))
            res['length']['mcs_shape'].append(list(map(float, (d_ / np.where(p_ > 0, p_, np.nan)) / (d_.sum() / p_.sum()))))
        for b in range(len(PB) - 1):
            mb, medb = stats(add(2, ('bias', b)))
            res['mc_bias'].append(dict(true_bin=[float(PB[b]), float(PB[b + 1])], mean=mb, median=medb))
        out[g] = res
    json.dump(out, open(os.path.join(C.OUTDIR, 'cr_mcs_check.json'), 'w'), indent=1, default=float)
    L = ['# MCS momentum of the muon candidate in the opened control regions', '',
         'Generated by `scripts/cr_mcs_check.py` from the control-region skim (`slurm/slurm_mp2pi_cr_skim.sbatch`).', '',
         '## Contained muon candidates: MCS / range momentum, median data / prediction (mean data / prediction)', '',
         '| range momentum [GeV/c] | ' + ' | '.join(out) + ' |', '|---|' + '---|' * len(out)]
    for b in range(len(RB) - 1):
        cells = []
        for g in out:
            r = out[g]['contained'][b]
            cells.append(f'{r["median_data"]:.3f} / {r["median_pred"]:.3f} ({r["mean_data"]:.3f} / {r["mean_pred"]:.3f}; n {r["n_data"]:.0f})')
        L.append(f'| {RB[b]:.2f}-{RB[b + 1]:.2f} | ' + ' | '.join(cells) + ' |')
    L += ['', '## Exiting muon candidates: MCS spectrum, data / prediction with the normalisation removed', '',
          '| MCS momentum [GeV/c] | ' + ' | '.join(f'{g} (norm {out[g]["exiting"]["norm"]:.2f})' for g in out) + ' |',
          '|---|' + '---|' * len(out)]
    for b in range(len(PB) - 1):
        cells = [f'{out[g]["exiting"]["shape"][b]:.2f} ± {out[g]["exiting"]["shape_err"][b]:.2f} (data {out[g]["exiting"]["data"][b]:.0f})' for g in out]
        L.append(f'| {PB[b]:.1f}-{PB[b + 1]:.1f} | ' + ' | '.join(cells) + ' |')
    L += ['', '## Simulation: MCS / true momentum of exiting muon candidates, median (mean)', '',
          '| true momentum [GeV/c] | ' + ' | '.join(out) + ' |', '|---|' + '---|' * len(out)]
    for b in range(len(PB) - 1):
        L.append(f'| {PB[b]:.1f}-{PB[b + 1]:.1f} | ' + ' | '.join(f'{out[g]["mc_bias"][b]["median"]:.3f} ({out[g]["mc_bias"][b]["mean"]:.3f})' for g in out) + ' |')
    L += ['', '## Reconstructed muon momentum in the released p_mu bins (all periods): data / prediction, normalisation removed', '',
          '| muon candidates | norm | ' + ' | '.join(f'{EB[b]:.3g}-{EB[b + 1]:.3g}' if b < len(EB) - 2 else f'> {EB[b]:.3g}' for b in range(len(EB) - 1)) + ' |',
          '|---|---|' + '---|' * (len(EB) - 1)]
    for tag in ('all', 'contained', 'exiting'):
        r = out['all']['pmu'][tag]
        L.append(f'| {tag} | {r["norm"]:.3f} | ' + ' | '.join(f'{x:.2f} ({d:.0f})' for x, d in zip(r['shape'], r['data'])) + ' |')
    r = out['all']['length']
    L += ['', '## Exiting muon candidates by visible length (all periods): data / prediction with the normalisation removed', '',
          '| visible length [cm] | length distribution | MCS 0-0.4 | 0.4-0.8 | 0.8-1.2 | 1.2-1.6 | > 1.6 GeV/c |',
          '|---|---|---|---|---|---|---|']
    for b in range(len(LB) - 1):
        L.append(f'| {LB[b]:.0f}-{LB[b + 1]:.0f} | {r["shape"][b]:.2f} (data {r["data"][b]:.0f}) | '
                 + ' | '.join(f'{x:.2f}' for x in r['mcs_shape'][b]) + ' |')
    open(os.path.join(C.OUTDIR, 'cr_mcs_check.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
