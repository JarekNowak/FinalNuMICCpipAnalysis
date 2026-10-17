"""mp2pi_kinematics.py -- feasibility of differential measurements with the two-pion selection (CC1mu2piBDT, MCS-calibrated
production /data/uboone/processed/mp2pi_mcscal): for each candidate variable, the predicted signal-region distribution by
true class at the full exposure (combined FHC+RHC), the response of the signal (true bin -> reco bin), the efficiency per
true bin and a binning under the 0.50 criterion (every diagonal element P(reco bin i | true bin i) >= 0.50) with at least
NMIN predicted signal events per bin.

Reco variables are built from the tracks the assignment chose (evt_mu_idx, evt_pi0_idx, evt_pi1_idx index the PFParticle
arrays). Muon momentum: range if the track end is inside the containment volume of the single-pion analysis, else MCS,
calibrated for simulated tracks with the v1.9 map (/data/uboone/processed/mcscal/maps/nominal.txt, the interpolation of
CC1mu1piXp::mp_mcs_mom). Angles about the NuMI beam: reco about the target-to-reco-vertex direction, truth about the true
neutrino direction (true_nu_p[xyz] of the raw ntuple, matched entry by entry), as include/XSecAnalyzer/NuMIBeamFrame.hh.
Pion momentum: range momentum of a contained track with the muon-to-pion mass correction (equal range at equal T/m).
The two pions are ordered by the variable itself (forward/backward pion, higher/lower momentum), the same way at truth
and reco level, so the convention is symmetric in the pions. Proton multiplicity: reco, an unassigned track with
P_p > 0.5 (event feature rest_nPp05); truth, a proton above 0.3 GeV/c (mc_n_proton_sub).

Scales: overlays and beam-off as bkgfit/farsb.py (per period), dirt per run (bkgfit/dirt.py). No beam-off sample exists
for Run 2: the Run-1 beam-off of the same horn mode stands in, scaled by the POT ratio (beam-off is 0.4% of the signal
region).

    python3 scripts/mp2pi_kinematics.py build    ->  /data/uboone/processed/mp2pi_mcscal/kin/events.npz
        (or per sample: build1 <k> for k = 0..N-1, in parallel, then merge)
    python3 scripts/mp2pi_kinematics.py report   ->  report/planning/multipion/mp2pi_kinematics.{json,md},
                                                     report/figures/multipion/kin_<variable>.pdf
"""
import json, os, re, sys
import numpy as np
import uproot, awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from bkgfit import farsb as F
from bkgfit.dirt import dirt_perrun

D = '/data/uboone/processed/mp2pi_mcscal/'
CACHE = D + 'kin/events.npz'
MAP = '/data/uboone/processed/mcscal/maps/nominal.txt'
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
FIGDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'figures', 'multipion'))
E = 'CC1mu2piBDT_'
CONT = (2.0, 254.35, -113.53, 107.47, 2.1, 1034.9)
MMU, MPI = 0.105658, 0.139570
# NuMI target in detector coordinates (NuMIBeamFrame.hh)
_R = np.array([[0.92103853804025681562, 4.6254001262154668408e-05, -0.38947144863934973769],
               [0.022713504803924120662, 0.99829162468141474651, 0.053832413938664107345],
               [0.38880857519374290021, -0.058427989452906302359, 0.91946400794392302291]])
TARGET = -(_R.T @ np.array([5502., 7259., 67270.]))
TOPO = ['signal', 'CC0pi', 'CC1pi', 'CC2pi other', 'CC3pi+', 'CCpi0', 'CC other', 'NC', 'nue CC', 'out of FV', 'beam-off', 'dirt']
NMIN = 40.
FEAT = ['dL_up', 'dL_down', 'dL_runnerup', 'L_best', 'mu_Pmu', 'mu_Ppi', 'mu_Pp', 'mu_Pother', 'mu_is_cand', 'rest_n',
        'rest_maxPp', 'rest_maxPpi', 'rest_nPp05']
NPP = FEAT.index('rest_nPp05')


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def load_map():
    cells = []
    for l in open(MAP):
        if not l.strip() or l.startswith('#'): continue
        v = [float(x) for x in l.split()]
        cells.append((v[0], v[1], v[2], v[3], np.array(v[4:])))
    return cells


def mcs_cal(mcs, adx, ln, cells):
    out = mcs.copy()
    for x0, x1, l0, l1, g in cells:
        m = (adx >= x0) & (adx < x1) & (ln >= l0) & (ln < l1) & (mcs > 0)
        if not m.any(): continue
        pmax = 0.01 * (len(g) - 1); p = mcs[m]
        u = np.minimum(p / 0.01, len(g) - 1 - 1e-9); k = u.astype(int); t = u - k
        out[m] = np.where(p >= pmax, p * g[-1] / pmax, g[k] + t * (g[np.minimum(k + 1, len(g) - 1)] - g[k]))
    return out


def inside(x, y, z):
    return (x > CONT[0]) & (x < CONT[1]) & (y > CONT[2]) & (y < CONT[3]) & (z > CONT[4]) & (z < CONT[5])


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n > 0, n, 1.)


def pi_range_mom(pmu):
    """range momentum of a pion from the muon-hypothesis range momentum (equal T/m at equal range)"""
    t = (np.sqrt(pmu ** 2 + MMU ** 2) - MMU) * MPI / MMU
    return np.sqrt((t + MPI) ** 2 - MPI ** 2)


def raw_path(out):
    log = D + os.path.basename(out)[len('xsec-ana-'):-len('.root')] + '.log'
    if not os.path.exists(log): return None
    m = re.search(r'input_filename: (\S+)', open(log).read())
    return m.group(1) if m else None


BR = ['region', 'MC_Signal', 'MC_Signal_1p', 'mc_topology', 'evt_mu_idx', 'evt_pi0_idx', 'evt_pi1_idx', 'evt_features',
      'mc_n_proton_sub']
TRK = ['trk_dir_x_v', 'trk_dir_y_v', 'trk_dir_z_v', 'trk_len_v', 'trk_sce_end_x_v', 'trk_sce_end_y_v', 'trk_sce_end_z_v',
       'trk_range_muon_mom_v', 'trk_mcs_muon_mom_v', 'backtracked_pdg']
VTX = ['reco_nu_vtx_sce_x', 'reco_nu_vtx_sce_y', 'reco_nu_vtx_sce_z']
W = ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
MC = ['mc_pdg', 'mc_px', 'mc_py', 'mc_pz']


def read(path, scale, kind, cells, out):
    is_mc = kind in ('mc', 'dirt')
    t = uproot.open(path)['stv_tree']
    br = [E + b if b != 'region' else E + 'evt_region' for b in BR] + TRK + VTX + (W + MC if is_mc else [])
    raw = None
    if kind == 'mc':
        rp = raw_path(path)
        raw = uproot.open(rp)['nuselection/NeutrinoSelectionFilter'] if rp else None
        assert raw is not None and raw.num_entries == t.num_entries, f'raw ntuple of {path}'
    start = 0
    for a in t.iterate(br, step_size=200000, library='ak'):
        n = len(a); reg = ak.to_numpy(a[E + 'evt_region'])
        sig = ak.to_numpy(a[E + 'MC_Signal']).astype(bool) if kind == 'mc' else np.zeros(n, bool)
        keep = (reg == 0) | sig
        idx = np.flatnonzero(keep)
        if len(idx) == 0: start += n; continue
        b = a[keep]; m = len(idx)
        w = (safe(ak.to_numpy(b['tuned_cv_weight']) * ak.to_numpy(b['ppfx_cv_weight']) * ak.to_numpy(b['normalisation_weight']))
             if is_mc else np.ones(m)) * scale
        rec = {k: np.full(m, np.nan) for k in ('pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'ppi_hi', 'ppi_lo',
                                               'pi_cont0', 'pi_cont1', 'mu_cont', 'np')}
        sel = ak.to_numpy(b[E + 'evt_region']) == 0
        if sel.any():
            bs = b[sel]
            vtx = np.stack([ak.to_numpy(bs[v]) for v in VTX], 1); nud = unit(vtx - TARGET)
            def trk(name, i):   # element i of each event's PFParticle array
                arr = bs[name]; n = ak.to_numpy(ak.num(arr)); st = np.concatenate([[0], np.cumsum(n)[:-1]])
                flat = ak.to_numpy(ak.flatten(arr)).astype(float); o = np.full(len(i), np.nan); ok = i < n
                o[ok] = flat[st[ok] + i[ok]]; return o
            im, i0, i1 = (np.maximum(ak.to_numpy(bs[E + x]), 0) for x in ('evt_mu_idx', 'evt_pi0_idx', 'evt_pi1_idx'))
            d = {}
            for tag, ii in (('mu', im), ('p0', i0), ('p1', i1)):
                d[tag] = unit(np.stack([trk(f'trk_dir_{c}_v', ii) for c in 'xyz'], 1))
                d[tag + 'cont'] = inside(trk('trk_sce_end_x_v', ii), trk('trk_sce_end_y_v', ii), trk('trk_sce_end_z_v', ii))
                d[tag + 'rng'] = trk('trk_range_muon_mom_v', ii)
                mcs = trk('trk_mcs_muon_mom_v', ii).astype(float)
                if is_mc:
                    sim = np.nan_to_num(trk('backtracked_pdg', ii)) != 0
                    cal = mcs_cal(mcs, np.abs(d[tag][:, 0]), trk('trk_len_v', ii), cells)
                    mcs = np.where(sim & ~d[tag + 'cont'], cal, mcs)
                d[tag + 'mcs'] = mcs
            r = {}
            r['pmu'] = np.where(d['mucont'], d['murng'], d['mumcs'])
            r['cmu'] = np.sum(d['mu'] * nud, 1)
            c0, c1 = np.sum(d['p0'] * nud, 1), np.sum(d['p1'] * nud, 1)
            r['cfwd'], r['cbwd'] = np.maximum(c0, c1), np.minimum(c0, c1)
            r['cpipi'] = np.sum(d['p0'] * d['p1'], 1)
            r['cmupi'] = np.maximum(np.sum(d['mu'] * d['p0'], 1), np.sum(d['mu'] * d['p1'], 1))
            q0 = np.where(d['p0cont'], pi_range_mom(d['p0rng']), np.nan); q1 = np.where(d['p1cont'], pi_range_mom(d['p1rng']), np.nan)
            r['ppi_hi'], r['ppi_lo'] = np.fmax(q0, q1), np.fmin(q0, q1)
            r['pi_cont0'], r['pi_cont1'], r['mu_cont'] = d['p0cont'].astype(float), d['p1cont'].astype(float), d['mucont'].astype(float)
            r['np'] = np.array([f[NPP] if len(f) > NPP else np.nan for f in ak.to_list(bs[E + 'evt_features'])], float)
            for k, v in r.items(): rec[k][sel] = v
        s = ak.to_numpy(b[E + 'MC_Signal']).astype(bool) if kind == 'mc' else np.zeros(m, bool)
        rows = dict(w=w, kind=np.full(m, {'mc': 0, 'ext': 1, 'dirt': 2}[kind]), region=ak.to_numpy(b[E + 'evt_region']), sig=s,
                    topo=ak.to_numpy(b[E + 'mc_topology']) if kind == 'mc' else np.full(m, 10 if kind == 'ext' else 11))
        for k, v in rec.items(): rows['r_' + k] = v
        # truth (signal rows only): final-state momenta, the true neutrino direction from the raw ntuple
        tv = {k: np.full(m, np.nan) for k in ('pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'ppi_hi', 'ppi_lo', 'np')}
        if s.any():
            bt = b[s]; ent = start + idx[s]
            pdg = bt['mc_pdg']; px, py, pz = bt['mc_px'], bt['mc_py'], bt['mc_pz']
            pm = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
            def vec(j):
                return np.stack([ak.to_numpy(ak.fill_none(ak.firsts(c[j]), 0.)) for c in (px, py, pz)], 1).astype(float)
            vmu = vec(ak.argmax(ak.where(np.abs(pdg) == 13, pm, -1.), axis=1, keepdims=True))
            order = ak.argsort(ak.where((np.abs(pdg) == 211) & (pm > 0), pm, -1.), axis=1, ascending=False)
            v0, v1 = vec(order[:, 0:1]), vec(order[:, 1:2])
            lo, hi = int(ent.min()), int(ent.max()) + 1
            nud = unit(np.stack([raw[f'true_nu_p{c}'].array(library='np', entry_start=lo, entry_stop=hi)[ent - lo]
                                 for c in 'xyz'], 1).astype(float))
            umu, u0, u1 = unit(vmu), unit(v0), unit(v1)
            tv['pmu'][s] = np.linalg.norm(vmu, axis=1)
            tv['cmu'][s] = np.sum(umu * nud, 1)
            c0, c1 = np.sum(u0 * nud, 1), np.sum(u1 * nud, 1)
            tv['cfwd'][s], tv['cbwd'][s] = np.maximum(c0, c1), np.minimum(c0, c1)
            tv['cpipi'][s] = np.sum(u0 * u1, 1)
            tv['cmupi'][s] = np.maximum(np.sum(umu * u0, 1), np.sum(umu * u1, 1))
            q0, q1 = np.linalg.norm(v0, axis=1), np.linalg.norm(v1, axis=1)
            tv['ppi_hi'][s], tv['ppi_lo'][s] = np.maximum(q0, q1), np.minimum(q0, q1)
            tv['np'][s] = ak.to_numpy(bt[E + 'mc_n_proton_sub'])
        for k, v in tv.items(): rows['t_' + k] = v
        for k, v in rows.items(): out.setdefault(k, []).append(v)
        start += n


def samples():
    mc, ext, _ = F.samples()
    re_ = lambda q: D + os.path.basename(q)
    out = []
    for ip, name in enumerate(F.PNAME):
        out += [(re_(f), s, 'mc') for f, s in mc[ip]]
        if ext[ip]: out += [(re_(f), s, 'ext') for f, s in ext[ip]]
        else:   # no Run-2 beam-off: Run 1 of the same horn mode, scaled by the POT ratio
            j = F.PNAME.index(name.replace('R2', 'R1'))
            out += [(re_(f), s * F.POT[ip] / F.POT[j], 'ext') for f, s in ext[j]]
    out += [(p, s, 'dirt') for _, p, s in dirt_perrun('comb', where=lambda run: D)]
    return out


def build():
    cells = load_map(); out = {}
    for p, s, kind in samples():
        print(f'  {kind:4s} {s:.5f} {os.path.basename(p)}', flush=True)
        read(p, s, kind, cells, out)
    ev = {k: np.concatenate(v) for k, v in out.items() if not k.startswith('_')}
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, **ev)
    print('written', CACHE, len(ev['w']), 'rows')


# ------------------------------------------------------------------------------------------------------ report
VARS = {
    'pmu':    dict(label=r'$p_\mu$ [GeV/$c$]', lo=0.15, hi=2.5, fine=0.05, over=True),
    'cmu':    dict(label=r'$\cos\theta_\mu$', lo=-1., hi=1., fine=0.05),
    'cfwd':   dict(label=r'$\cos\theta_\pi$, forward pion', lo=-1., hi=1., fine=0.05),
    'cbwd':   dict(label=r'$\cos\theta_\pi$, backward pion', lo=-1., hi=1., fine=0.05),
    'cpipi':  dict(label=r'$\cos\theta_{\pi\pi}$', lo=-1., hi=1., fine=0.05),
    'cmupi':  dict(label=r'$\cos\theta_{\mu\pi}$, closer pion', lo=-1., hi=1., fine=0.05),
    'ppi_lo': dict(label=r'$p_\pi$, lower (both contained) [GeV/$c$]', lo=0.10, hi=0.60, fine=0.025, over=True),
    'ppi_hi': dict(label=r'$p_\pi$, higher (both contained) [GeV/$c$]', lo=0.10, hi=0.80, fine=0.025, over=True),
    'np':     dict(label='protons above 0.3 GeV/$c$ (reco: unassigned track with $P_p>0.5$)', edges=[-0.5, 0.5, 20.5]),
}


def binning(t, r, w, fine_edges, nmin=NMIN, crit=0.50):
    """contiguous merge of fine bins, left to right, until the bin holds >= nmin signal and its diagonal >= crit"""
    edges = [fine_edges[0]]; i = 0; nf = len(fine_edges) - 1
    while i < nf:
        j = i + 1
        while True:
            lo, hi = edges[-1], fine_edges[j]
            tin = (t >= lo) & (t < hi); both = tin & (r >= lo) & (r < hi)
            if (w[tin].sum() >= nmin and w[both].sum() >= crit * w[tin].sum()) or j == nf: break
            j += 1
        edges.append(fine_edges[j]); i = j
    if len(edges) > 2:   # the last bin may fail: merge it into its neighbour
        lo, hi = edges[-2], edges[-1]; tin = (t >= lo) & (t < hi); both = tin & (r >= lo) & (r < hi)
        if w[tin].sum() < nmin or w[both].sum() < crit * w[tin].sum(): edges.pop(-2)
    return np.array(edges)


def report():
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11})
    ev = dict(np.load(CACHE))
    w, sr, sig = ev['w'], ev['region'] == 0, ev['sig']
    res, L = {}, ['# Differential feasibility of the two-pion selection', '',
                  'Generated by `scripts/mp2pi_kinematics.py` (see its docstring). Combined FHC+RHC, full exposure, CC1mu2piBDT '
                  'signal region (MCS-calibrated production). Binning: contiguous, left to right, each bin with at least '
                  f'{NMIN:.0f} predicted signal events in the signal region and P(reco bin i | true bin i) >= 0.50. '
                  'Stat.: sqrt(N)/S of the bin (the data statistics of a background-subtracted count, before unfolding).', '',
                  f'Signal region: {w[sr].sum():.1f} events, {w[sr & sig].sum():.1f} signal; signal generated {w[sig & (ev["kind"] == 0)].sum():.1f}.', '']
    L += ['| variable | reco available in SR signal | bins | edges | diagonal | efficiency | signal per bin | stat. per bin |',
          '|---|---|---|---|---|---|---|---|']
    gen = sig & (ev['kind'] == 0)
    for v, c in VARS.items():
        t, r = ev['t_' + v], ev['r_' + v]
        if v == 'np': fine = np.array(c['edges'])
        else:
            fine = np.round(np.arange(c['lo'], c['hi'] + 1e-9, c['fine']), 4)
            if c.get('over'): fine = np.append(fine, 50.)
        sel = sr & sig & np.isfinite(t) & np.isfinite(r)
        tt, rr = t.copy(), r.copy()
        if c.get('over'): tt = np.minimum(tt, 49.9); rr = np.minimum(rr, 49.9)
        avail = w[sel].sum() / w[sr & sig].sum()
        edges = fine if v == 'np' else binning(tt[sel], rr[sel], w[sel], fine)
        nb = len(edges) - 1
        M = np.histogram2d(tt[sel], rr[sel], [edges, edges], weights=w[sel])[0]
        rown = M.sum(1)
        diag = np.diag(M) / np.where(rown > 0, rown, 1)
        # efficiency per true bin: selected (any reco value) over generated
        gsel = gen & np.isfinite(t); tg = np.minimum(t, 49.9) if c.get('over') else t
        G = np.histogram(tg[gsel], edges, weights=w[gsel])[0]
        Ssel = np.histogram(tg[gsel & sr & np.isfinite(r)], edges, weights=w[gsel & sr & np.isfinite(r)])[0]
        eff = Ssel / np.where(G > 0, G, 1)
        okr = sr & np.isfinite(r); rr2 = np.minimum(r, 49.9) if c.get('over') else r
        N = np.histogram(rr2[okr], edges, weights=w[okr])[0]; S = np.histogram(rr2[okr & sig], edges, weights=w[okr & sig])[0]
        stat = np.sqrt(N) / np.where(S > 0, S, 1)
        fmt = lambda a: ', '.join(f'{x:g}' for x in a)
        res[v] = dict(edges=edges.tolist(), diag=diag.tolist(), eff=eff.tolist(), signal=S.tolist(), events=N.tolist(),
                      stat=stat.tolist(), reco_available=avail)
        L.append(f"| {v} | {100 * avail:.0f}% | {nb} | {fmt(edges)} | {fmt(np.round(diag, 2))} | {fmt(np.round(100 * eff, 1))}% | "
                 f"{fmt(np.round(S))} | {fmt(np.round(100 * stat))}% |")
        # figure: predicted events per measurement bin by true class, response of the signal, efficiency and sqrt(N)/S
        fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3), gridspec_kw=dict(width_ratios=[1.15, 1, 0.9]))
        lab = [f'{a:g}$-${b:g}' if not (c.get('over') and k == nb - 1) else f'$>${a:g}'
               for k, (a, b) in enumerate(zip(edges[:-1], edges[1:]))]
        if v == 'np': lab = ['0', r'$\geq1$']
        cls = [('signal', sig), ('CC1$\\pi$', ev['topo'] == 2), ('CC$\\geq3\\pi$', ev['topo'] == 4),
               ('CC$\\pi^0$', ev['topo'] == 5), ('CC0$\\pi$', ev['topo'] == 1), ('NC', ev['topo'] == 7),
               ('other $\\nu$', ~sig & np.isin(ev['topo'], [3, 6, 8, 9])), ('dirt, beam-off', np.isin(ev['topo'], [10, 11]))]
        bot = np.zeros(nb); xb = np.arange(nb)
        for (n, m), col in zip(cls, plt.cm.tab10.colors):
            h = np.histogram(rr2[okr & m], edges, weights=w[okr & m])[0]
            ax[0].bar(xb, h, 0.85, bottom=bot, color=col, label=n, edgecolor='k', lw=0.3); bot += h
        ax[0].set_xticks(xb); ax[0].set_xticklabels(lab, rotation=40, fontsize=8.5, ha='right')
        ax[0].set_xlabel('reco bin of ' + c['label'], fontsize=10); ax[0].set_ylabel('predicted events, signal region')
        ax[0].set_ylim(0, 1.32 * bot.max())
        ax[0].legend(fontsize=7.5, ncol=4, loc='upper left', frameon=False, columnspacing=0.8, handlelength=1.2)
        Mn = M / np.where(rown[:, None] > 0, rown[:, None], 1)
        ax[1].imshow(Mn.T, origin='lower', cmap='Blues', vmin=0, vmax=1, aspect='auto')
        for i in range(nb):
            for j in range(nb):
                ax[1].text(i, j, f'{Mn[i, j]:.2f}', ha='center', va='center', fontsize=7.5 if nb > 6 else 9,
                           color='w' if Mn[i, j] > 0.6 else 'k')
        ax[1].set_xticks(xb); ax[1].set_yticks(xb)
        ax[1].set_xticklabels(lab, rotation=40, fontsize=8, ha='right'); ax[1].set_yticklabels(lab, fontsize=8)
        ax[1].set_xlabel('true bin'); ax[1].set_ylabel('reco bin'); ax[1].set_title('P(reco bin | true bin), signal', fontsize=10)
        ax[2].bar(xb, 100 * eff, 0.6, color='#55a868', label='efficiency [%]')
        ax[2].plot(xb, 100 * stat, 'o-', color='#c44e52', label=r'$\sqrt{N}/S$ [%]')
        ax[2].set_xticks(xb); ax[2].set_xticklabels(lab, rotation=40, fontsize=8, ha='right')
        ax[2].set_ylim(0, max(30, 1.3 * max(100 * stat.max(), 100 * eff.max()))); ax[2].legend(fontsize=8.5, loc='upper left')
        ax[2].grid(alpha=0.3); ax[2].set_xlabel('true bin (efficiency), reco bin ($\\sqrt{N}/S$)', fontsize=9)
        fig.text(0.005, 0.995, 'MicroBooNE simulation, internal; combined FHC+RHC, full exposure', fontsize=8, color='0.45', va='top')
        fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(f'{FIGDIR}/kin_{v}.pdf'); plt.close(fig)
    L += ['', 'Reco available: share of the signal-region signal for which the reco variable exists (pion momenta: both '
          'assigned pions contained).']
    json.dump(res, open(os.path.join(OUTDIR, 'mp2pi_kinematics.json'), 'w'), indent=1)
    open(os.path.join(OUTDIR, 'mp2pi_kinematics.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


def build1(k):
    """one sample (index k of samples()) into kin/part_<k>.npz, for parallel builds; merge with 'merge'"""
    cells = load_map(); out = {}; p, s, kind = samples()[int(k)]
    read(p, s, kind, cells, out)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(os.path.join(os.path.dirname(CACHE), f'part_{int(k):02d}.npz'),
                        **{a: np.concatenate(v) for a, v in out.items()})
    print('done', k, os.path.basename(p), flush=True)


def merge():
    n = len(samples()); parts = [os.path.join(os.path.dirname(CACHE), f'part_{k:02d}.npz') for k in range(n)]
    missing = [q for q in parts if not os.path.exists(q)]
    if missing: sys.exit(f'missing parts: {missing}')
    P = [x for x in (dict(np.load(q)) for q in parts) if "w" in x]   # a sample with no kept event writes no rows
    np.savez_compressed(CACHE, **{k: np.concatenate([x[k] for x in P]) for k in P[0]})
    print('written', CACHE)


if __name__ == '__main__':
    if sys.argv[1] == 'build1': build1(sys.argv[2])
    else: {'build': build, 'report': report, 'merge': merge}[sys.argv[1]]()
