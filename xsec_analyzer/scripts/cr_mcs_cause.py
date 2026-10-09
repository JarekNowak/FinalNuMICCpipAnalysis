"""cr_mcs_cause.py -- where the MCS momentum of beam-on muon candidates departs from the prediction (scripts/cr_mcs_check.py:
exiting candidates far low above about 1.2 GeV/c, contained candidates agree up to 1 GeV/c). Opened control regions of
the single-pion analysis, all periods, prediction = overlay + beam-off + dirt with the scales of scripts/cr_mp2pi_check.py.

(1) Contained candidates up to 3 GeV/c range momentum: MCS / range, median data and prediction. At high momentum the
    scattering angles approach the angular resolution, so an estimator problem shows here with the momentum known.
(2) Contained candidates: the trajectory deflection spread (trk_avg_deflection_stdev_v) times the range momentum, median
    data and prediction per range-momentum bin. Pure scattering makes the product flat; a worse angular resolution in
    data raises it at high momentum.
(3) Exiting candidates: deflection spread, data over prediction with the normalisation removed, and its median per bin
    of visible length and of the drift-direction component. (2) and (3) use the periods whose beam-on files carry the
    deflection branch (new production), data and prediction alike.
(4) Exiting candidates: the fraction with MCS momentum above 1.2 GeV/c, data over prediction (1 = no deficit), in bins
    of the track direction, start position, visible length and the face through which the track leaves.
(5) The Run-4 detector-variation samples, unweighted as they enter the detector term: the same fraction, the MCS / true
    momentum and the deflection spread of exiting candidates, variation over CV.
(6) Tracks longer than 1 m: the angle between the start direction and the start-to-end chord and the chord / length.
    A coherent bend (a distortion) raises the angle; noise on the scale of the MCS segments does not.

    python3 scripts/cr_mcs_cause.py  ->  report/planning/multipion/cr_mcs_cause.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun
import cr_mp2pi_check as C

S = 'CC1mu1piXp'
RB = np.array([0.15, 0.35, 0.5, 0.7, 1.0, 1.25, 1.5, 2.0, 3.0])
QB = np.linspace(0., 2.5, 251)
DB = np.logspace(-6, 0, 241)                   # deflection spread [rad], positive values only
DP = np.logspace(-6, 0, 241)                   # deflection spread x momentum [rad GeV/c]
DV = np.array([1e-6, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1.])   # coarse deflection-spread spectrum
DVAR = ['LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']
TB = np.array([0.2, 0.6, 1.0, 1.5, 2.0, 3.0])
AB = np.linspace(0., 45., 451)                 # angle between start direction and chord [deg]
SB = np.linspace(0.8, 1.1, 301)                # chord / length
NEWPROD = ['FHC_R4', 'FHC_R5', 'RHC_R4']        # periods whose beam-on data carry trk_avg_deflection_stdev_v  # true momentum of exiting true muons in the detector variations
LB = np.array([0., 50., 100., 150., 200., 300., 1000.])
PHI = 1.2                                      # GeV/c, high-MCS threshold
SPLITS = {
    '|dir x| (drift)': ('adx', [0., 0.3, 0.6, 0.85, 1.01]),
    'dir y (vertical)': ('dy', [-1.01, -0.5, 0., 0.5, 1.01]),
    'dir z (beam axis)': ('dz', [-1.01, 0., 0.5, 0.8, 0.95, 1.01]),
    'theta_yz to the z axis [deg]': ('tyz', [0., 15., 45., 75., 105., 180.]),
    'start x [cm]': ('sx', [-10., 64., 128., 192., 270.]),
    'start y [cm]': ('sy', [-130., -60., 0., 60., 130.]),
    'start z [cm]': ('sz', [-10., 250., 500., 750., 1050.]),
    'visible length [cm]': ('len', list(LB)),
}
FACES = ['x low (anode)', 'x high (cathode)', 'y low', 'y high', 'z low', 'z high']
TPC = np.array([[0., 256.35], [-116.5, 116.5], [0., 1036.8]])


def first(arr, idx):
    """element idx of every row of a jagged array (nan when absent)"""
    sel = ak.local_index(arr) == idx
    return ak.to_numpy(ak.fill_none(ak.firsts(arr[sel]), np.nan)).astype(float)


def fill(path, scale, weighted, H):
    t = uproot.open(path)['stv_tree']
    tv = ['trk_len_v', 'trk_sce_start_x_v', 'trk_sce_start_y_v', 'trk_sce_start_z_v', 'trk_sce_end_x_v', 'trk_sce_end_y_v',
          'trk_sce_end_z_v', 'trk_dir_x_v', 'trk_dir_y_v', 'trk_dir_z_v', 'trk_avg_deflection_stdev_v', 'trk_mcs_muon_mom_v']
    br = [f'{S}_CandidateMuonIndex', f'{S}_CandidateMuonTrackEndContainment', f'{S}_candidate_muon_mom_mcs',
          f'{S}_candidate_muon_mom_range'] + tv + (['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight'] if weighted else [])
    for a in t.iterate(br, step_size='200 MB', library='ak'):
        if len(a) == 0: continue
        w = C.safe(ak.to_numpy(a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight'])) if weighted \
            else np.ones(len(a))
        w = w * scale
        mi = ak.to_numpy(a[f'{S}_CandidateMuonIndex'])
        ok = mi >= 0
        v = {k: first(a[k], ak.Array(np.where(ok, mi, 0))) for k in tv}
        mcs = ak.to_numpy(a[f'{S}_candidate_muon_mom_mcs']).astype(float)
        rng = ak.to_numpy(a[f'{S}_candidate_muon_mom_range']).astype(float)
        con = ak.to_numpy(a[f'{S}_CandidateMuonTrackEndContainment']).astype(bool)
        H['check'] = H.get('check', 0) + np.array([np.sum(ok), np.sum(ok & np.isclose(v['trk_mcs_muon_mom_v'], mcs, rtol=1e-4))])
        dfl = v['trk_avg_deflection_stdev_v']
        # (1), (2): contained
        c = ok & con & (rng > RB[0]) & (rng < RB[-1]) & np.isfinite(mcs) & (mcs > 0)
        ib = np.digitize(rng, RB) - 1
        for b in range(len(RB) - 1):
            m = c & (ib == b)
            H[('q', b)] = H.get(('q', b), 0) + np.histogram(np.clip(mcs[m] / rng[m], QB[0], QB[-1] - 1e-9), QB, weights=w[m])[0]
            H[('dz', b)] = H.get(('dz', b), 0) + np.array([w[m].sum(), w[m & ~(dfl > 0)].sum()])
            m &= dfl > 0
            H[('dp', b)] = H.get(('dp', b), 0) + np.histogram(np.clip(dfl[m] * rng[m], DP[0], DP[-1] * 0.999), DP, weights=w[m])[0]
        # (3), (4): exiting
        e = ok & ~con & np.isfinite(mcs) & (mcs > 0)
        ed = e & (dfl > 0)
        H['edz'] = H.get('edz', 0) + np.array([w[e].sum(), w[e & ~(dfl > 0)].sum()])
        H['edfl'] = H.get('edfl', 0) + np.histogram(np.clip(dfl[ed], DV[0], DV[-1] * 0.999), DV, weights=w[ed])[0]
        ln = v['trk_len_v']
        for b in range(len(LB) - 1):
            m = ed & (ln >= LB[b]) & (ln < LB[b + 1])
            H[('ldfl', b)] = H.get(('ldfl', b), 0) + np.histogram(np.clip(dfl[m], DB[0], DB[-1] * 0.999), DB, weights=w[m])[0]
        for b in range(len(SPLITS['|dir x| (drift)'][1]) - 1):
            bins = SPLITS['|dir x| (drift)'][1]
            m = ed & (np.abs(v['trk_dir_x_v']) >= bins[b]) & (np.abs(v['trk_dir_x_v']) < bins[b + 1])
            H[('xdfl', b)] = H.get(('xdfl', b), 0) + np.histogram(np.clip(dfl[m], DB[0], DB[-1] * 0.999), DB, weights=w[m])[0]
        dx, dy, dz = v['trk_dir_x_v'], v['trk_dir_y_v'], v['trk_dir_z_v']
        q = dict(adx=np.abs(dx), dy=dy, dz=dz, tyz=np.degrees(np.arccos(np.clip(dz / np.sqrt(dy ** 2 + dz ** 2 + 1e-12), -1, 1))),
                 sx=v['trk_sce_start_x_v'], sy=v['trk_sce_start_y_v'], sz=v['trk_sce_start_z_v'], len=ln)
        hi = mcs > PHI
        for name, (k, bins) in SPLITS.items():
            x = q[k]
            for b in range(len(bins) - 1):
                m = e & (x >= bins[b]) & (x < bins[b + 1])
                H[('split', name, b)] = H.get(('split', name, b), 0) + np.array([w[m].sum(), w[m & hi].sum(), (w[m] ** 2).sum(), (w[m & hi] ** 2).sum()])
        end = np.stack([v['trk_sce_end_x_v'], v['trk_sce_end_y_v'], v['trk_sce_end_z_v']])
        dist = np.concatenate([end - TPC[:, :1], TPC[:, 1:] - end])           # x lo, y lo, z lo, x hi, y hi, z hi
        face = np.array([0, 2, 4, 1, 3, 5])[np.argmin(np.where(np.isfinite(dist), dist, 1e9), axis=0)]
        for f in range(6):
            m = e & (face == f)
            H[('face', f)] = H.get(('face', f), 0) + np.array([w[m].sum(), w[m & hi].sum(), (w[m] ** 2).sum(), (w[m & hi] ** 2).sum()])
        # (6) straightness of tracks longer than 1 m: angle between the start direction and the start-to-end chord
        ch = end - np.stack([v['trk_sce_start_x_v'], v['trk_sce_start_y_v'], v['trk_sce_start_z_v']])
        cl = np.linalg.norm(ch, axis=0)
        ang = np.degrees(np.arccos(np.clip((ch / np.where(cl > 0, cl, np.nan) * np.stack([dx, dy, dz])).sum(0), -1, 1)))
        st = cl / np.where(ln > 0, ln, np.nan)
        bx = SPLITS['|dir x| (drift)'][1]
        for tag, m0 in (('exit', e), ('cont', ok & con & (rng > 0.8))):
            for b in range(len(bx) - 1):
                m = m0 & (ln > 100) & np.isfinite(ang) & (np.abs(dx) >= bx[b]) & (np.abs(dx) < bx[b + 1])
                H[('ang', tag, b)] = H.get(('ang', tag, b), 0) + np.histogram(np.clip(ang[m], 0, 44.99), AB, weights=w[m])[0]
                H[('st', tag, b)] = H.get(('st', tag, b), 0) + np.histogram(np.clip(st[m], SB[0], SB[-1] * 0.9999), SB, weights=w[m])[0]


def med(h, edges):
    h = np.asarray(h, float)
    if not h.sum() > 0: return float('nan')
    return float(np.interp(0.5, np.cumsum(h) / h.sum(), edges[1:]))


def detvar_one(path):
    br = [f'{S}_CandidateMuonIndex', f'{S}_CandidateMuonTrackEndContainment', f'{S}_candidate_muon_mom_mcs',
          f'{S}_candidate_muon_mom_true', 'trk_avg_deflection_stdev_v']
    H = {'n': np.zeros(2), 'r': [np.zeros(len(QB) - 1) for _ in range(len(TB) - 1)], 'd': np.zeros(len(DB) - 1)}
    for a in uproot.open(path)['stv_tree'].iterate(br, step_size='300 MB', library='ak'):
        mi = ak.to_numpy(a[f'{S}_CandidateMuonIndex']); ok = mi >= 0
        mcs = ak.to_numpy(a[f'{S}_candidate_muon_mom_mcs']).astype(float)
        tr = ak.to_numpy(a[f'{S}_candidate_muon_mom_true']).astype(float)
        e = ok & ~ak.to_numpy(a[f'{S}_CandidateMuonTrackEndContainment']).astype(bool) & np.isfinite(mcs) & (mcs > 0)
        H['n'] += [e.sum(), (e & (mcs > PHI)).sum()]
        for b in range(len(TB) - 1):
            m = e & (tr >= TB[b]) & (tr < TB[b + 1])
            H['r'][b] += np.histogram(np.clip(mcs[m] / tr[m], QB[0], QB[-1] - 1e-9), QB)[0]
        dfl = first(a['trk_avg_deflection_stdev_v'], ak.Array(np.where(ok, mi, 0)))
        m = e & (dfl > 0)
        H['d'] += np.histogram(np.clip(dfl[m], DB[0], DB[-1] * 0.999), DB)[0]
    return H


def detvar(out):
    """(5) the Run-4 detector-variation samples (unweighted, as they enter the detector term): exiting muon candidates,
    high-MCS fraction, MCS / true momentum of the candidate and the deflection spread, variation over its CV"""
    L = ['', '## (5) Detector variations (Run 4): exiting muon candidates, variation / CV', '',
         f'Fraction with MCS > {PHI} GeV/c, median MCS / true momentum per true-momentum bin, and median deflection spread. '
         'The control-region data need a high-MCS fraction ratio of about 0.49 and MCS / true lowered by f(p)/p '
         '(0.78 at 1.5, 0.73 at 2.0 GeV/c).', '']
    out['detvar'] = {}
    for cfg in ('fhc', 'rhc'):
        P = '/data/uboone/processed/xsec-ana-detvar_run4' + cfg + '_'
        cv = detvar_one(P + 'CV.root')
        hd = '| ' + cfg.upper() + ' | exiting | high-MCS fraction | ratio | ' + ' | '.join(
            f'MCS/true {TB[b]:g}-{TB[b + 1]:g}' for b in range(len(TB) - 1)) + ' | deflection |'
        L += [hd, '|' + '---|' * (hd.count('|') - 1)]
        fcv = cv['n'][1] / cv['n'][0]
        qcv = [med(h, QB) for h in cv['r']]; dcv = med(cv['d'], DB)
        L.append(f'| CV | {cv["n"][0]:.0f} | {fcv:.3f} | 1 | ' + ' | '.join(f'{q:.3f}' for q in qcv) + f' | {1e3 * dcv:.3f} mrad |')
        out['detvar'][cfg] = {'CV': dict(n=cv['n'].tolist(), q=qcv, d=dcv)}
        for v in DVAR:
            if not os.path.exists(P + v + '.root'): continue
            h = detvar_one(P + v + '.root')
            f = h['n'][1] / h['n'][0]; q = [med(x, QB) for x in h['r']]; d = med(h['d'], DB)
            out['detvar'][cfg][v] = dict(n=h['n'].tolist(), q=q, d=d)
            L.append(f'| {v} | {h["n"][0]:.0f} | {f:.3f} | {f / fcv:.3f} | ' + ' | '.join(f'{a / b:.3f}' for a, b in zip(q, qcv))
                     + f' | {d / dcv:.3f} |')
            print('detvar', cfg, v, flush=True)
        L.append('')
    return L


def main():
    mc, ext, data = F.samples()
    re = lambda p: C.NEW + os.path.basename(p)
    dirts = dirt_perrun('comb', where=lambda run: C.NEW)
    Hd, Hp, Nd, Np = {}, {}, {}, {}
    for ip, name in enumerate(F.PNAME):
        if not data[ip]: continue
        hd, hp = {}, {}
        for f in data[ip]: fill(re(f), 1.0, False, hd)
        for f, s in mc[ip]: fill(re(f), s, True, hp)
        for f, s in ext[ip]: fill(re(f), s, False, hp)
        for r, f, s in dirts:
            if r == F.RUN[ip]: fill(f, s, True, hp)
        for H, T in ((hd, Hd), (hp, Hp)) + (((hd, Nd), (hp, Np)) if name in NEWPROD else ()):
            for k, x in H.items(): T[k] = T.get(k, 0) + x
        print('==', name, flush=True)
    out = {'index_check_data': Hd['check'].tolist(), 'index_check_pred': Hp['check'].tolist()}
    L = ['# Where the MCS momentum of beam-on muon candidates departs from the prediction', '',
         'Generated by `scripts/cr_mcs_cause.py` from the control-region skim (`slurm/slurm_mp2pi_cr_skim.sbatch`), all '
         'periods; prediction = overlay + beam-off + dirt. Muon-candidate track variables taken from the PFP vectors at '
         f'the candidate index (MCS momentum matches the candidate branch for {Hd["check"][1]:.0f} of {Hd["check"][0]:.0f} '
         'data candidates).', '',
         '## (1) Contained candidates: MCS / range momentum, median (data / prediction), and (2) deflection spread x range momentum, median', '',
         'Deflection spread: trk_avg_deflection_stdev_v [rad], positive values only, periods FHC Run 4, FHC Run 5 and RHC '
         'Run 4 (the old-production beam-on files lack the branch); the share of candidates without a positive value is listed.', '',
         '| range momentum [GeV/c] | data events | MCS / range data | prediction | ratio | deflection x p data [mrad GeV/c] | prediction | ratio | no deflection data | prediction |',
         '|---|---|---|---|---|---|---|---|---|---|']
    out['contained'] = []
    for b in range(len(RB) - 1):
        qd, qp = med(Hd.get(('q', b), 0), QB), med(Hp.get(('q', b), 0), QB)
        dd, dp = med(Nd.get(('dp', b), 0), DP), med(Np.get(('dp', b), 0), DP)
        zd, zp = Nd[('dz', b)][1] / Nd[('dz', b)][0], Np[('dz', b)][1] / Np[('dz', b)][0]
        n = float(np.sum(Hd.get(('q', b), 0)))
        out['contained'].append(dict(bin=[float(RB[b]), float(RB[b + 1])], n=n, q_data=qd, q_pred=qp, dp_data=dd, dp_pred=dp,
                                     nodefl_data=zd, nodefl_pred=zp))
        L.append(f'| {RB[b]:.2f}-{RB[b + 1]:.2f} | {n:.0f} | {qd:.3f} | {qp:.3f} | {qd / qp:.3f} | {1e3 * dd:.3f} | {1e3 * dp:.3f} | '
                 f'{dd / dp:.3f} | {zd:.3f} | {zp:.3f} |')
    d, p = Nd['edfl'], Np['edfl']
    sh = (d / np.where(p > 0, p, np.nan)) / (d.sum() / p.sum())
    out['exiting_deflection'] = dict(bins=DV.tolist(), shape=sh.tolist())
    L += ['', '## (3) Exiting candidates: deflection spread [mrad], data / prediction with the normalisation removed', '',
          f'Without a positive deflection spread: data {Nd["edz"][1] / Nd["edz"][0]:.3f}, prediction {Np["edz"][1] / Np["edz"][0]:.3f}.', '',
          '| ' + ' | '.join(f'{1e3 * DV[i]:.3g}-{1e3 * DV[i + 1]:.3g}' for i in range(len(DV) - 1)) + ' |',
          '|' + '---|' * (len(DV) - 1), '| ' + ' | '.join(f'{x:.2f}' for x in sh) + ' |',
          '', '| visible length [cm] | median deflection data [mrad] | prediction | ratio |', '|---|---|---|---|']
    out['exiting_deflection_length'] = []
    for b in range(len(LB) - 1):
        md_, mp_ = med(Nd.get(('ldfl', b), 0), DB), med(Np.get(('ldfl', b), 0), DB)
        out['exiting_deflection_length'].append([md_, mp_])
        L.append(f'| {LB[b]:.0f}-{LB[b + 1]:.0f} | {1e3 * md_:.3f} | {1e3 * mp_:.3f} | {md_ / mp_:.3f} |')
    L += ['', '| abs(dir x) | median deflection data [mrad] | prediction | ratio |', '|---|---|---|---|']
    bx = SPLITS['|dir x| (drift)'][1]
    out['exiting_deflection_dirx'] = []
    for b in range(len(bx) - 1):
        md_, mp_ = med(Nd.get(('xdfl', b), 0), DB), med(Np.get(('xdfl', b), 0), DB)
        out['exiting_deflection_dirx'].append([md_, mp_])
        L.append(f'| {bx[b]:g} to {bx[b + 1]:g} | {1e3 * md_:.3f} | {1e3 * mp_:.3f} | {md_ / mp_:.3f} |')
    # (4)
    tot_d = sum(Hd[('face', f)] for f in range(6)); tot_p = sum(Hp[('face', f)] for f in range(6))
    L += ['', f'## (4) Exiting candidates: fraction with MCS momentum > {PHI} GeV/c, data / prediction', '',
          f'Inclusive: data {tot_d[1] / tot_d[0]:.3f}, prediction {tot_p[1] / tot_p[0]:.3f}, ratio '
          f'{(tot_d[1] / tot_d[0]) / (tot_p[1] / tot_p[0]):.3f}. Per bin: data fraction / prediction fraction = ratio (data events).', '']
    out['split'] = {}
    for name, (k, bins) in list(SPLITS.items()) + [('exit face', (None, None))]:
        L += [f'| {name} | ' + ' | '.join(['bin', 'data', 'prediction', 'ratio', 'events']) + ' |', '|---|---|---|---|---|---|']
        rows = []
        nb = 6 if k is None else len(bins) - 1
        for b in range(nb):
            dd = Hd.get(('face', b) if k is None else ('split', name, b), np.zeros(4))
            pp = Hp.get(('face', b) if k is None else ('split', name, b), np.zeros(4))
            if dd[0] <= 0 or pp[0] <= 0: continue
            fd, fp = dd[1] / dd[0], pp[1] / pp[0]
            er = (fd / fp) * np.sqrt(1 / max(dd[1], 1) + 1 / dd[0])
            lab = FACES[b] if k is None else f'{bins[b]:g} to {bins[b + 1]:g}'
            rows.append(dict(bin=lab, f_data=fd, f_pred=fp, ratio=fd / fp, err=er, n=float(dd[0])))
            L.append(f'| | {lab} | {fd:.3f} | {fp:.3f} | {fd / fp:.2f} ± {er:.2f} | {dd[0]:.0f} |')
        out['split'][name] = rows
        L.append('')
    L += ['## (6) Straightness of tracks longer than 1 m: angle between the start direction and the start-to-end chord', '',
          'Median angle [deg] and median chord / length, data / prediction. Exiting candidates, and contained candidates with '
          'range momentum above 0.8 GeV/c.', '',
          '| candidates | abs(dir x) | data events | angle data | prediction | chord / length data | prediction |',
          '|---|---|---|---|---|---|---|']
    out['straightness'] = {}
    for tag in ('exit', 'cont'):
        for b in range(len(bx) - 1):
            hd, hp = Hd.get(('ang', tag, b), np.zeros(1)), Hp.get(('ang', tag, b), np.zeros(1))
            if not np.sum(hd) > 0: continue
            r = [med(hd, AB), med(hp, AB), med(Hd[('st', tag, b)], SB), med(Hp[('st', tag, b)], SB)]
            out['straightness'][f'{tag} {b}'] = r
            L.append(f'| {"exiting" if tag == "exit" else "contained > 0.8 GeV/c"} | {bx[b]:g} to {bx[b + 1]:g} | {np.sum(hd):.0f} | '
                     f'{r[0]:.2f} | {r[1]:.2f} | {r[2]:.4f} | {r[3]:.4f} |')
    L += detvar(out)
    json.dump(out, open(os.path.join(C.OUTDIR, 'cr_mcs_cause.json'), 'w'), indent=1, default=float)
    open(os.path.join(C.OUTDIR, 'cr_mcs_cause.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
