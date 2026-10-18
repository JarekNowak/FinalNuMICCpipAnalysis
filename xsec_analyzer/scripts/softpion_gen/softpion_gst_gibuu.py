"""softpion_gst_gibuu.py -- softest charged-pion momentum of the two-pion signal in standalone GENIE (gst) and GiBUU
(FinalEvents.dat), the samples of generator_predictions/newg4 (NuMI flux, numu and numubar, argon). Signal as
CC1mu2piBDT, kinematic part: CC, one muon above 0.15 GeV/c, exactly two charged pions (at any momentum), no pi0, kaon
or heavier meson; signal with both pions above 0.10 GeV/c, sub-threshold with the softer below. Weighted counts per
bin of scripts/mp2pi_softpion.py (GEDGES) into /data/uboone/processed/mp2pi_2f/softpion/gen_<generator>_<flavour>.json.

    python3 scripts/softpion_gen/softpion_gst_gibuu.py
"""
import json, os
import numpy as np
import uproot, awkward as ak

G = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/'
OUT = '/data/uboone/processed/mp2pi_2f/softpion/'
GEDGES = [0.0, 0.05, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.60, 1e9]


def is_meson(pdg):
    a = np.abs(pdg)
    return (a < 9900000) & ((a // 1000) % 10 == 0) & ((a // 100) % 10 != 0) & ~((a >= 901) & (a <= 930)) & \
        ~np.isin(a, [110, 990, 998, 999, 100])


def save(name, psoft, w, n_events):
    h = np.histogram(np.clip(psoft, 0, 1e8), GEDGES, weights=w)[0]
    json.dump(dict(edges=GEDGES, counts=h.tolist(), n_events=int(n_events)), open(f'{OUT}gen_{name}.json', 'w'))
    print(name, 'events', n_events, 'two-pion (any p_soft)', round(float(w.sum()), 3), 'share 0.10-0.175 of signal',
          round(float(h[2:5].sum() / h[2:].sum()), 4), flush=True)


def genie(fl):
    t = uproot.open(f'{G}g_{fl}.gst.root')['gst']
    a = t.arrays(['cc', 'fspl', 'pxl', 'pyl', 'pzl', 'pdgf', 'pxf', 'pyf', 'pzf'])
    pdg = a['pdgf']; p = np.sqrt(a['pxf'] ** 2 + a['pyf'] ** 2 + a['pzf'] ** 2)
    pmu = np.sqrt(a['pxl'] ** 2 + a['pyl'] ** 2 + a['pzl'] ** 2)
    ab = abs(pdg)
    npipm = ak.sum(ab == 211, axis=1); npi0 = ak.sum(pdg == 111, axis=1)
    nk = ak.sum((ab == 321) | (ab == 311) | (pdg == 310) | (pdg == 130), axis=1)
    nheavy = ak.sum(ak.unflatten(is_meson(ak.to_numpy(ak.flatten(pdg))) & ~np.isin(ak.to_numpy(ak.flatten(ab)), [111, 211, 321, 311, 310, 130]),
                            ak.num(pdg)), axis=1)
    sel = ak.to_numpy(a['cc']).astype(bool) & (np.abs(ak.to_numpy(a['fspl'])) == 13) & (ak.to_numpy(pmu) > 0.15) & \
        (ak.to_numpy(npipm) == 2) & (ak.to_numpy(npi0) == 0) & (ak.to_numpy(nk) == 0) & (ak.to_numpy(nheavy) == 0)
    ps = ak.to_numpy(ak.min(ak.where(ab == 211, p, 1e9), axis=1))[sel]
    save(f'genie_{fl}', ps, np.ones(len(ps)), t.num_entries)


def gibuu(fl):
    d = np.loadtxt(f'{G}gibuu/run_{fl}/FinalEvents.dat', comments='#', usecols=(0, 1, 2, 3, 4, 9, 10, 11))
    run, evt, pid, ch, pw = d[:, 0].astype(np.int64), d[:, 1].astype(np.int64), d[:, 2].astype(int), d[:, 3].astype(int), d[:, 4]
    p = np.sqrt(d[:, 5] ** 2 + d[:, 6] ** 2 + d[:, 7] ** 2)
    key = run * 10 ** 9 + evt
    u, inv = np.unique(key, return_inverse=True); n = len(u)
    cnt = lambda m: np.bincount(inv, weights=m.astype(float), minlength=n)
    is_mu = pid == 902
    nmu = cnt(is_mu); npipm = cnt((pid == 101) & (ch != 0)); npi0 = cnt((pid == 101) & (ch == 0))
    nk = cnt((pid == 110) | (pid == 111)); nheavy = cnt(((pid >= 102) & (pid <= 109)) | ((pid >= 112) & (pid <= 122)))
    pmu = np.zeros(n); w = np.zeros(n); pmu[inv[is_mu]] = p[is_mu]; w[inv[is_mu]] = pw[is_mu]
    psoft = np.full(n, 1e9); pi = (pid == 101) & (ch != 0)
    np.minimum.at(psoft, inv[pi], p[pi])
    sel = (nmu == 1) & (pmu > 0.15) & (npipm == 2) & (npi0 == 0) & (nk == 0) & (nheavy == 0)
    save(f'gibuu_{fl}', psoft[sel], w[sel], n)


if __name__ == '__main__':
    for fl in ('numu', 'numubar'):
        genie(fl); gibuu(fl)
