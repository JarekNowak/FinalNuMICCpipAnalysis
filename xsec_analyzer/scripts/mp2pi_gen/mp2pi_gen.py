"""mp2pi_gen.py -- generator predictions of the two-pion signal (CC1mu2piBDT signal definition, kinematic part) for the
differential variables and the one-bin total, in the true bins of the extractions, as input to the framework
(Prediction ... file <fte file> <variable>_fte).

Samples of generator_predictions/newg4 (NuMI flux, numu and numubar, argon): GENIE (gst), GiBUU (FinalEvents.dat),
NuWro and NEUT (signal events dumped by nuwro_2pi.cc and neut_2pi.C, run_nuwro_neut.sh). Signal: CC, exactly one muon
with p > 0.15 GeV/c, exactly two charged pions (counted at any momentum) both above 0.10 GeV/c, no pi0, kaon or heavier
meson, any number of nucleons, opening angle of the muon and the leading pion below 2.6 rad. Variables as
CC1mu1piXp::fill_mp_true_kinematics, angles about the neutrino direction (+z in these samples); N_p counts protons above
0.3 GeV/c. Bins: the edges of report/planning/multipion/mp2pi_kinematics_2f.json, first bin open below, last open above.

Normalisation per nucleon in cm2, as the single-pion readers of generator_predictions/gen2d: GENIE sigma_CC per Ar / 40
times the signal fraction; GiBUU sum of perweight / 20 runs x 1e-38; NuWro the event weight and NEUT the harmonic mean of
Totcrs (1e-38) times the signal fraction. Per generator and flavour a file with one unit-width histogram per variable
(bin-integrated cross section per nucleon); generator_predictions/gen2d/combine_ext.C then mixes numu and numubar by
flux, multiplies by 40 nucleons and converts to 1e-38 cm2/Ar (FHC, RHC), and combine_comb_ext the combined exposure.

Horn mode: the samples were generated with the FHC fluxes of each flavour (newg4/flux_newg4_clean.root = the *_cv_fhc
histograms of flux_newg4.root). For RHC every event is weighted by r(E) = (phi_RHC(E)/Phi_RHC) / (phi_FHC(E)/Phi_FHC)
of its flavour; since the events follow phi_FHC sigma, sigma_FHC sum(r)/N is the RHC flux-averaged cross section.
(The released single-pion predictions use the FHC-flux samples for RHC unchanged.)

    python3 scripts/mp2pi_gen/mp2pi_gen.py   ->  /data/uboone/processed/mp2pi_2f/gen/<generator>_2pi_{numu,numubar}_{fhc,rhc}.root,
                                                 <generator>_2pi_{fhc,rhc,comb}_fte.root
"""
import json, os, subprocess
import numpy as np
import uproot, awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'report', 'planning', 'multipion'))
G = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/newg4/'
GEN2D = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/generator_predictions/gen2d/'
OUT = '/data/uboone/processed/mp2pi_2f/gen/'
VARS = ['pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'np']
GENIE_SIGMA_AR = {'numu': 1.390129e-37, 'numubar': 4.208102e-38}   # as generator_predictions/gen2d/regen_ext.sh


def flux_ratio(fl):
    """r(E) of the flavour: RHC flux shape over FHC flux shape (above the 0.06 GeV generation threshold)"""
    f = uproot.open(G + 'flux_newg4.root')
    hf, hr = f[f'{fl}_cv_fhc'], f[f'{fl}_cv_rhc']; e = hf.axis().edges(); c = 0.5 * (e[1:] + e[:-1]); m = c > 0.06
    vf, vr = hf.values() * m, hr.values() * m
    r = np.where(vf > 0, (vr / vr.sum()) / np.where(vf > 0, vf / vf.sum(), 1.), 0.)
    return lambda E: r[np.clip(np.digitize(E, e) - 1, 0, len(r) - 1)]


def is_meson(pdg):
    a = np.abs(pdg)
    return (a < 9900000) & ((a // 1000) % 10 == 0) & ((a // 100) % 10 != 0) & ~((a >= 901) & (a <= 930)) & \
        ~np.isin(a, [110, 990, 998, 999, 100])


def genie(fl):
    t = uproot.open(f'{G}g_{fl}.gst.root')['gst']
    a = t.arrays(['cc', 'fspl', 'pxl', 'pyl', 'pzl', 'pdgf', 'pxf', 'pyf', 'pzf', 'Ev'])
    pdg = a['pdgf']; ab = abs(pdg)
    p = np.sqrt(a['pxf'] ** 2 + a['pyf'] ** 2 + a['pzf'] ** 2)
    flat = ak.to_numpy(ak.flatten(pdg))
    heavy = ak.unflatten(is_meson(flat) & ~np.isin(np.abs(flat), [111, 211, 321, 311, 310, 130]), ak.num(pdg))
    sel = ak.to_numpy(a['cc']).astype(bool) & (np.abs(ak.to_numpy(a['fspl'])) == 13) & \
        (ak.to_numpy(ak.sum(ab == 211, axis=1)) == 2) & (ak.to_numpy(ak.sum(pdg == 111, axis=1)) == 0) & \
        (ak.to_numpy(ak.sum((ab == 321) | (ab == 311) | (pdg == 310) | (pdg == 130), axis=1)) == 0) & \
        (ak.to_numpy(ak.sum(heavy, axis=1)) == 0)
    s = a[sel]; pis = ab[sel] == 211
    pv = np.stack([ak.to_numpy(s[k][pis]) for k in ('pxf', 'pyf', 'pzf')], axis=-1)        # (n, 2, 3)
    mu = np.stack([ak.to_numpy(s[k]) for k in ('pxl', 'pyl', 'pzl')], axis=1)
    npr = ak.to_numpy(ak.sum((pdg[sel] == 2212) & (p[sel] > 0.3), axis=1))
    w = np.full(len(mu), GENIE_SIGMA_AR[fl] / 40. / t.num_entries)
    return mu, pv[:, 0], pv[:, 1], npr, w, ak.to_numpy(s['Ev'])


def gibuu(fl):
    d = np.loadtxt(f'{G}gibuu/run_{fl}/FinalEvents.dat', comments='#', usecols=(0, 1, 2, 3, 4, 9, 10, 11, 14))
    key = d[:, 0].astype(np.int64) * 10 ** 9 + d[:, 1].astype(np.int64)
    pid, ch, pw, v, enu = d[:, 2].astype(int), d[:, 3].astype(int), d[:, 4], d[:, 5:8], d[:, 8]
    u, inv = np.unique(key, return_inverse=True); n = len(u)
    cnt = lambda m: np.bincount(inv, weights=m.astype(float), minlength=n)
    ismu, ispi = pid == 902, (pid == 101) & (ch != 0)
    sel = (cnt(ismu) == 1) & (cnt(ispi) == 2) & (cnt((pid == 101) & (ch == 0)) == 0) & (cnt((pid == 110) | (pid == 111)) == 0) & \
        (cnt(((pid >= 102) & (pid <= 109)) | ((pid >= 112) & (pid <= 122))) == 0)
    mu = np.zeros((n, 3)); w = np.zeros(n); E = np.zeros(n); mu[inv[ismu]] = v[ismu]; w[inv[ismu]] = pw[ismu]; E[inv[ismu]] = enu[ismu]
    pirows = np.flatnonzero(ispi); ev = inv[pirows]
    order = np.argsort(ev, kind='stable'); pirows, ev = pirows[order], ev[order]
    first = np.full(n, -1); second = np.full(n, -1)
    uev, idx, c = np.unique(ev, return_index=True, return_counts=True)
    first[uev] = pirows[idx]; two = c >= 2; second[uev[two]] = pirows[idx[two] + 1]
    npr = cnt((pid == 1) & (ch == 1) & (np.linalg.norm(v, axis=1) > 0.3))
    s = np.flatnonzero(sel)
    return mu[s], v[first[s]], v[second[s]], npr[s].astype(int), w[s] / 20. * 1e-38, E[s]


def dumped(g, fl):
    path = f'{OUT}gen2pi_{g}_{fl}.txt'
    head = {l.split()[1]: float(l.split()[2]) for l in open(path) if l.startswith('#')}
    a = np.loadtxt(path, comments='#', ndmin=2)
    return a[:, 0:3], a[:, 3:6], a[:, 6:9], a[:, 9].astype(int), np.full(len(a), head['sigma'] / head['N']), a[:, 10]


def variables(mu, p0, p1, npr):
    m0, m1 = np.linalg.norm(p0, axis=1), np.linalg.norm(p1, axis=1)
    lead = np.where((m0 >= m1)[:, None], p0, p1); pmu = np.linalg.norm(mu, axis=1)
    um, u0, u1 = mu / pmu[:, None], p0 / m0[:, None], p1 / m1[:, None]
    ok = (pmu > 0.15) & (np.minimum(m0, m1) > 0.10) & \
        (np.sum(um * lead / np.linalg.norm(lead, axis=1)[:, None], axis=1) > np.cos(2.6))
    c0, c1 = u0[:, 2], u1[:, 2]
    x = dict(pmu=pmu, cmu=um[:, 2], cfwd=np.maximum(c0, c1), cbwd=np.minimum(c0, c1), cpipi=np.sum(u0 * u1, axis=1),
             cmupi=np.maximum(np.sum(um * u0, axis=1), np.sum(um * u1, axis=1)), np=npr.astype(float))
    return ok, x


def write(g, fl, mode, data, edges):
    import ROOT
    mu, p0, p1, npr, w, E = data
    if mode == 'rhc': w = w * flux_ratio(fl)(E)
    ok, x = variables(mu, p0, p1, npr)
    f = ROOT.TFile(f'{OUT}{g}_2pi_{fl}_{mode}.root', 'recreate'); keep = []
    out = {'total': float(w[ok].sum())}
    h = ROOT.TH1D('total', '', 1, 0, 1); h.SetBinContent(1, w[ok].sum()); h.SetBinError(1, np.sqrt((w[ok] ** 2).sum())); keep.append(h)
    for v in VARS:
        e = edges[v]; n = len(e) - 1
        b = np.clip(np.digitize(x[v], e) - 1, 0, n - 1)
        s = np.bincount(b[ok], weights=w[ok], minlength=n); s2 = np.bincount(b[ok], weights=w[ok] ** 2, minlength=n)
        h = ROOT.TH1D(v, '', n, 0, n)
        for i in range(n): h.SetBinContent(i + 1, s[i]); h.SetBinError(i + 1, np.sqrt(s2[i]))
        keep.append(h); out[v] = (s / s.sum()).tolist()
    f.Write(); f.Close()
    print(f'{g:6s} {fl:8s} {mode} signal events {ok.sum():6d}, sigma {1e38 * w[ok].sum():.4f} 1e-38 cm2/nucleon', flush=True)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    edges = {v: json.load(open(os.path.join(PLAN, 'mp2pi_kinematics_2f.json')))[v]['edges'] for v in VARS}
    res = {}
    for g in ('genie', 'gibuu', 'nuwro', 'neut'):
        for fl in ('numu', 'numubar'):
            data = genie(fl) if g == 'genie' else gibuu(fl) if g == 'gibuu' else dumped(g, fl)
            for mode in ('fhc', 'rhc'): res[f'{g} {fl} {mode}'] = write(g, fl, mode, data, edges)
        obs = ','.join(['total'] + VARS)
        macro = (f'.L {GEN2D}combine_ext.C+', f'combine_ext("{OUT}{g}_2pi_numu_fhc.root","{OUT}{g}_2pi_numubar_fhc.root","{OUT}{g}_2pi_fhc_fte.root","fhc","{obs}"); '
                 f'combine_ext("{OUT}{g}_2pi_numu_rhc.root","{OUT}{g}_2pi_numubar_rhc.root","{OUT}{g}_2pi_rhc_fte.root","rhc","{obs}"); '
                 f'combine_comb_ext("{OUT}{g}_2pi_fhc_fte.root","{OUT}{g}_2pi_rhc_fte.root","{OUT}{g}_2pi_comb_fte.root","{obs}")')
        r = subprocess.run(['root', '-l', '-b', '-q', '-e', macro[0], '-e', macro[1]], capture_output=True, text=True)
        print('\n'.join(l for l in r.stdout.splitlines() if 'total' in l or 'missing' in l or 'wrote' in l))
    json.dump(res, open(f'{OUT}gen2pi_shapes.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
