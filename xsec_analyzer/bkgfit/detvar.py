"""bkgfit.detvar -- detector-variation shifts with the statistical part removed.

Three facts about the Run-4 detector-variation samples, measured 2026-09-16, drive this module:
  1. The knob samples do NOT have the POT of the CV sample (WMX has 0.913 of it in FHC, Recomb2 1.039).
     Unscaled count ratios therefore carry that factor as a spurious "shift": -8% in the FHC control
     regions for most knobs, against a POT-normalised shift of about +-1%. (The extraction normalises by
     summed_pot; the sideband macros sb_protocol*.C, sb_detvar*.C and the first bkgfit templates did not.)
  2. About 97% of the events of a knob sample are the SAME generated events as in the CV sample
     (matched on true neutrino energy and vertex): the samples are re-simulations, not independent.
  3. The statistical noise of a knob-minus-CV difference therefore comes only from events that migrate
     between bins (and the unmatched few per cent), 3-10x less than the independent-sample estimate.

Method, per beam mode and knob k, on the template bins (plus the virtual region totals):
  * shift  d_k,i = N_k,i * POT_CV / POT_k - N_CV,i , fractional f_k,i = d_k,i / N_CV,i;
  * noise  C_k = covariance of f_k over B Poisson-bootstrap replicas in which every GENERATED event gets one
           Poisson(1) weight shared by its CV and knob copies (unmatched events get their own);
  * because E[f f^T] = phi phi^T + C, the unbiased estimate of the true knob covariance is f f^T - C;
    summed over knobs and projected onto the positive semi-definite cone (eigenvalue clipping).
Usage:  python3 -m bkgfit.detvar   -> /data/uboone/processed/bkgfit/detvar_matched.npz
"""
import os, sys
import numpy as np
import uproot

KNOBS = ['LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']
NBV = 48                          # 46 template bins + 2 virtual totals (see bkgfit.model)
CTH = np.array([-1, 0.0, 0.45, 0.65, 0.8, 0.9, 0.95, 1.0]); PMU = np.array([0.15, 0.35, 0.75])
OUT = '/data/uboone/processed/bkgfit/detvar_matched.npz'
S = 'CC1mu1piXp'


def load(path):
    f = uproot.open(path); t = f['stv_tree']
    br = ['mc_nu_energy', 'mc_nu_vtx_x', 'mc_nu_vtx_y', 'mc_nu_vtx_z', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight',
          f'{S}_Selected', f'{S}_MC_Signal', f'{S}_sb_cc0pi', f'{S}_sb_pi0', f'{S}_sb_multipi', f'{S}_sb_cosmic',
          f'{S}_candidate_muon_costh_reco', f'{S}_candidate_muon_mom_reco', f'{S}_EventCategory']
    a = t.arrays(br, library='np')
    pot = float(f['summed_pot'].member('fVal'))
    w = a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']
    w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.0)
    key = np.core.records.fromarrays([np.round(a['mc_nu_energy'], 5), np.round(a['mc_nu_vtx_x'], 3),
                                      np.round(a['mc_nu_vtx_y'], 3), np.round(a['mc_nu_vtx_z'], 3)])
    c, p = a[f'{S}_candidate_muon_costh_reco'], a[f'{S}_candidate_muon_mom_reco']
    ic = np.searchsorted(CTH, c, side='right') - 1; ic = np.where(c == 1.0, 6, ic)
    ip = np.searchsorted(PMU, p, side='right') - 1
    kin = np.where((ic >= 0) & (ic < 7) & (ip >= 0), ic * 3 + ip, -1)
    pairs_e, pairs_b = [], []
    def add(mask, b):
        idx = np.nonzero(mask)[0]; pairs_e.append(idx); pairs_b.append(b if np.ndim(b) else np.full(len(idx), b))
    add(a[f'{S}_Selected'].astype(bool), 0)
    cc0 = a[f'{S}_sb_cc0pi'].astype(bool) & (kin >= 0); pi0 = a[f'{S}_sb_pi0'].astype(bool) & (kin >= 0)
    add(cc0, 2 + kin[cc0]); add(pi0, 23 + kin[pi0])
    add(a[f'{S}_sb_multipi'].astype(bool), 44); add(a[f'{S}_sb_cosmic'].astype(bool), 45)
    add(cc0, 46); add(pi0, 47)
    return dict(key=key, w=w, pot=pot, ev=np.concatenate(pairs_e), bin=np.concatenate(pairs_b).astype(int), n=len(w))


def matched_ids(kcv, kal):
    allk = np.concatenate([kcv, kal])
    _, inv = np.unique(allk, return_inverse=True)
    return inv[:len(kcv)], inv[len(kcv):], inv.max() + 1


def knob_shift(cv, al, B=300, seed=1):
    idc, ida, nid = matched_ids(cv['key'], al['key'])
    s = cv['pot'] / al['pot']
    ncv = np.bincount(cv['bin'], weights=cv['w'][cv['ev']], minlength=NBV)
    nal = np.bincount(al['bin'], weights=al['w'][al['ev']], minlength=NBV) * s
    f = np.divide(nal - ncv, ncv, out=np.zeros(NBV), where=ncv > 0)
    rng = np.random.default_rng(seed); reps = np.zeros((B, NBV))
    ec, ea = idc[cv['ev']], ida[al['ev']]
    for b in range(B):
        r = rng.poisson(1.0, nid).astype(float)
        c_b = np.bincount(cv['bin'], weights=cv['w'][cv['ev']] * r[ec], minlength=NBV)
        a_b = np.bincount(al['bin'], weights=al['w'][al['ev']] * r[ea], minlength=NBV) * s
        reps[b] = np.divide(a_b - c_b, c_b, out=np.zeros(NBV), where=c_b > 0)
    C = np.cov(reps.T)
    matched = len(np.intersect1d(idc, ida)) / al['n']
    return f, C, ncv, dict(pot_ratio=1 / s, matched=matched)


def main():
    out = {}
    for mode in ('run4fhc', 'run4rhc'):
        cv = load(f'/data/uboone/processed/sb_pi0/xsec-ana-detvar_{mode}_CV.root')
        F, C, info = [], [], []
        for k in KNOBS:
            al = load(f'/data/uboone/processed/sb_pi0/xsec-ana-detvar_{mode}_{k}.root')
            f, c, ncv, inf = knob_shift(cv, al)
            F.append(f); C.append(c); info.append(inf)
            tot = [46, 47, 0]
            print(f"{mode} {k:10s} POT ratio {inf['pot_ratio']:.4f} matched {inf['matched']:.3f} | "
                  + '  '.join(f"{n}: {f[b]:+.4f} +- {np.sqrt(c[b, b]):.4f}" for n, b in zip(('CC0pi', 'pi0', 'SR'), tot)))
        out[f'{mode}_f'] = np.array(F); out[f'{mode}_noise'] = np.array(C); out[f'{mode}_cv'] = ncv
        out[f'{mode}_pot_ratio'] = np.array([i['pot_ratio'] for i in info]); out[f'{mode}_matched'] = np.array([i['matched'] for i in info])
    os.makedirs(os.path.dirname(OUT), exist_ok=True); np.savez(OUT, **out); print('wrote', OUT)


if __name__ == '__main__':
    main()
