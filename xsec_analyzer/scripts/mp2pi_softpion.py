"""mp2pi_softpion.py -- does the interaction-model covariance cover the softest-pion dependence of the two-pion signal
efficiency? Two-fold production (/data/uboone/processed/mp2pi_2f), overlays at full exposure (scales of
scripts/mp2pi_kinematics.py), true signal (CC1mu2piBDT_MC_Signal) in bins of the softest charged-pion momentum
p_soft (mc_pionpm_min_mom), selected = signal region (evt_region 0). Also the sub-threshold class: two charged pions
with the softer below 0.10 GeV/c (mc_topology 3, fewer than two above 0.10), which migrates into the signal region
as background if the soft-pion spectrum changes.

Weight columns, as UniverseMaker: CV (tuned_cv * ppfx_cv * normalisation); the softest-pion stress model of the
ensembles (CV * (1 + 0.5 cos(pi min(p_soft, 0.5) / 0.5)), signal only, and "extended" also on the sub-threshold
class); the 600 GENIE multisim universes (weight_All_UBGenie * ppfx * normalisation, replacing the tune); the 1000
reinteraction universes (weight_reint_all * CV); the unisim knobs (weight_*_UBGenie * ppfx * normalisation) and the
xsr_scc universes (* CV).

    sbatch --array=0-<n-1> slurm/slurm_mp2pi_softpion.sbatch     (python3 scripts/mp2pi_softpion.py build1 <k>)
    python3 scripts/mp2pi_softpion.py report   ->  report/planning/multipion/mp2pi_softpion.{md,json},
                                                   report/figures/multipion/softpion.pdf

Generators: the softest-pion spectrum of the two-pion signal in standalone GENIE, GiBUU, NuWro and NEUT
(scripts/softpion_gen/, gen_<generator>_<flavour>.json), numu and numubar combined with the antineutrino share of the
CV signal (F_BAR), folded with the CV efficiency per bin: the efficiency each spectrum implies, against the CV and
against standalone GENIE (the same model as the CV, so the difference is the model's), and the sub-threshold class
relative to the signal, scaled to the CV prediction in the signal region.

Soft-pion term (generator envelope): the CV signal reweighted per p_soft bin by the spectrum of each generator over that
of standalone GENIE (numu and numubar mixed with the antineutrino share of the configuration's CV signal); the term is
the largest efficiency shift over GiBUU, NuWro and NEUT, for the one-bin total per configuration (CV efficiency per bin
of that horn mode) and per true bin of each differential variable (combined, kin/events.npz of the production).

    python3 scripts/mp2pi_softpion.py term     ->  report/planning/multipion/mp2pi_softpion_term.{md,json}

Covariance of the differential results: in the extraction a model that changes the efficiency per true bin by delta
moves the unfolded result by A_C (delta * t), t the tune truth per true bin (bin-integrated; the event counts of the
CV signal per true bin of kin/events.npz, normalised so that A_C t reproduces the tune of the extraction,
vec_table_uBTune) and A_C the regularisation matrix of that extraction (mat_table_add_smear). The term is the
envelope of the generators, in every bin the largest of their shifts (GiBUU and NEUT shift every bin the same way),
fully correlated over the bins: cov = dx dx^T, for the combined extractions without (unfold/<v>/comb_none) and with
(unfold/<v>/comb) the constraint. Shifts are quoted in % of the A_C-smeared tune.

    python3 scripts/mp2pi_softpion.py cov      ->  /data/uboone/processed/mp2pi_2f/unfold/<v>/softpion_cov_<none|sidebands>.txt,
                                                   report/planning/multipion/mp2pi_softpion_cov.{md,json}
"""
import json, os, sys
import numpy as np
if '--prod' not in sys.argv: sys.argv += ['--prod', '2f']
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp2pi_kinematics as K
import uproot, awkward as ak

OUT = '/data/uboone/processed/mp2pi_2f/softpion/'
E = 'CC1mu2piBDT_'
EDGES = np.array([0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.60, 1e9])
MEDGES = np.array([0.0, 0.05, 0.10])
UNISIM = ['Theta_Delta2Npi', 'ThetaDelta2NRad', 'DecayAngMEC', 'NormCCCOH', 'RPA_CCQE', 'AxFFCCQEshape',
          'VecFFCCQEshape', 'XSecShape_CCMEC']
SCC = ['xsr_scc_Fa3_SCC', 'xsr_scc_Fv3_SCC']
NMS, NRE = 600, 1000
F_BAR = 0.125      # antineutrino share of the CV two-pion signal, combined, two-fold production (2026-10-18)
F_BAR_CFG = {'fhc5': 0.132, 'rhcfull': 0.121, 'comb': 0.125}   # per configuration, same production
GENS = ['genie', 'gibuu', 'nuwro', 'neut']


def stress(p):
    return 1 + 0.5 * np.cos(np.pi * np.minimum(p, 0.5) / 0.5)


def columns():
    c = ['CV', 'stress', 'stress_ext'] + [f'ms{i}' for i in range(NMS)] + [f're{i}' for i in range(NRE)]
    for k in UNISIM: c += [f'{k}_0', f'{k}_1']
    for k in SCC: c += [f'{k}_{i}' for i in range(10)]
    return c


def mc_samples():
    return [x for x in K.samples() if x[2] == 'mc']


def build1(k):
    p, s, _ = mc_samples()[int(k)]
    ncol = len(columns()); nb, nm = len(EDGES) - 1, len(MEDGES) - 1
    acc = dict(S_all=np.zeros((ncol, nb)), S_sel=np.zeros((ncol, nb)), M_all=np.zeros((ncol, nm)),
               M_sel=np.zeros((ncol, nm)), pool3=np.zeros(nb), S_raw=np.zeros(nb), S_sel_raw=np.zeros(nb))
    base = [E + x for x in ('MC_Signal', 'mc_topology', 'mc_pionpm_min_mom', 'mc_n_pipm_thr010', 'evt_region',
                            'n_track_pool')] + K.W
    wbr = ['weight_All_UBGenie', 'weight_reint_all'] + [f'weight_{u}_UBGenie' for u in UNISIM] + [f'weight_{u}' for u in SCC]
    t = uproot.open(p)['stv_tree']; test = os.environ.get('SOFTPION_TEST')   # test: the first 100000 entries only
    for a in t.iterate(base + wbr, step_size=100000, library='ak', entry_stop=100000 if test else None):
        sig = ak.to_numpy(a[E + 'MC_Signal']).astype(bool)
        mig = (ak.to_numpy(a[E + 'mc_topology']) == 3) & (ak.to_numpy(a[E + 'mc_n_pipm_thr010']) < 2)
        keep = sig | mig
        if not keep.any(): continue
        a = a[keep]; sig = sig[keep]; mig = mig[keep]
        ps = ak.to_numpy(a[E + 'mc_pionpm_min_mom']).astype(float)
        sel = ak.to_numpy(a[E + 'evt_region']) == 0
        pn = K.safe(ak.to_numpy(a['ppfx_cv_weight']) * ak.to_numpy(a['normalisation_weight']))
        cv = K.safe(ak.to_numpy(a['tuned_cv_weight']) * ak.to_numpy(a['ppfx_cv_weight']) * ak.to_numpy(a['normalisation_weight']))
        # universe weights as UniverseMaker: safe_weight() of the raw weight, then the CV corrections
        arr = lambda br, n: K.safe(ak.to_numpy(ak.fill_none(ak.pad_none(a[br], n, clip=True), 1.0)).astype(float))
        cols = [cv, cv * np.where(sig, stress(ps), 1.0), cv * stress(ps)]
        W = [np.column_stack(cols), arr('weight_All_UBGenie', NMS) * pn[:, None], arr('weight_reint_all', NRE) * cv[:, None]]
        for u in UNISIM: W.append(arr(f'weight_{u}_UBGenie', 2) * pn[:, None])
        for u in SCC: W.append(arr(f'weight_{u}', 10) * cv[:, None])
        W = np.hstack(W) * s
        W[~np.isfinite(W)] = 0.
        bs = np.clip(np.digitize(ps, EDGES) - 1, 0, nb - 1); bm = np.clip(np.digitize(ps, MEDGES) - 1, 0, nm - 1)
        Hs = np.zeros((len(ps), nb)); Hs[np.arange(len(ps)), bs] = sig
        Hm = np.zeros((len(ps), nm)); Hm[np.arange(len(ps)), bm] = mig
        acc['S_all'] += W.T @ Hs; acc['S_sel'] += W.T @ (Hs * sel[:, None])
        acc['M_all'] += W.T @ Hm; acc['M_sel'] += W.T @ (Hm * sel[:, None])
        p3 = ak.to_numpy(a[E + 'n_track_pool']) >= 3
        acc['pool3'] += (Hs * (cv * s * p3)[:, None]).sum(0)
        acc['S_raw'] += Hs.sum(0); acc['S_sel_raw'] += (Hs * sel[:, None]).sum(0)
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(f'{OUT}{"test_" if test else ""}part_{int(k):02d}.npz', **acc)
    print('done', k, os.path.basename(p), 'signal', acc['S_raw'].sum(), flush=True)


def report():
    n = len(mc_samples()); P = [dict(np.load(f'{OUT}part_{k:02d}.npz')) for k in range(n)]
    A = {key: sum(x[key] for x in P) for key in P[0]}
    col = columns(); ix = {c: i for i, c in enumerate(col)}
    ms = [ix[f'ms{i}'] for i in range(NMS)]; re_ = [ix[f're{i}'] for i in range(NRE)]
    eps = A['S_sel'].sum(1) / A['S_all'].sum(1); d = eps / eps[0] - 1          # efficiency shift per column
    eb = A['S_sel'][0] / A['S_all'][0]                                          # CV efficiency per p_soft bin
    fold = (A['S_all'] * eb).sum(1) / A['S_all'].sum(1); df = fold / fold[0] - 1   # spectrum-only part
    lo = A['S_all'][:, :3].sum(1) / A['S_all'].sum(1)                           # share below 0.175 GeV/c
    # restricted phase space: both pions above 0.175 GeV/c
    epsr = A['S_sel'][:, 3:].sum(1) / A['S_all'][:, 3:].sum(1); dr = epsr / epsr[0] - 1
    mig = A['M_sel'].sum(1); sr_sig = A['S_sel'].sum(1)
    q = lambda x: (float(np.mean(x)), float(np.std(x)), float(np.percentile(x, 16)), float(np.percentile(x, 84)))
    res = dict(edges=EDGES[:-1].tolist(), eff_cv=eb.tolist(), signal_cv=A['S_all'][0].tolist(),
               signal_sel_cv=A['S_sel'][0].tolist(), pool3_cv=(A['pool3'] / A['S_all'][0]).tolist(),
               raw_signal=A['S_raw'].tolist(), raw_selected=A['S_sel_raw'].tolist(),
               eff_total_cv=float(eps[0]), stress=dict(eff_shift=float(d[1]), fold_shift=float(df[1]), share_lo=float(lo[1]),
               restricted_eff_shift=float(dr[1]), migration_sr=float(mig[1]), migration_sr_ext=float(mig[2])),
               cv=dict(share_lo=float(lo[0]), migration_sr=float(mig[0]), sr_signal=float(sr_sig[0])),
               multisim=dict(eff_shift=q(d[ms]), fold_shift=q(df[ms]), share_lo=q(lo[ms]), restricted_eff_shift=q(dr[ms]),
                             migration_sr=q(mig[ms]), sr_signal=q(sr_sig[ms])),
               reint=dict(eff_shift=q(d[re_]), fold_shift=q(df[re_]), share_lo=q(lo[re_]), restricted_eff_shift=q(dr[re_]),
                          migration_sr=q(mig[re_])),
               unisim={u: dict(eff_shift=[float(d[ix[f'{u}_0']]), float(d[ix[f'{u}_1']])],
                               share_lo=[float(lo[ix[f'{u}_0']]), float(lo[ix[f'{u}_1']])],
                               migration_sr=[float(mig[ix[f'{u}_0']]), float(mig[ix[f'{u}_1']])]) for u in UNISIM},
               scc={u: dict(eff_shift=q(d[[ix[f'{u}_{i}'] for i in range(10)]])) for u in SCC},
               ms_eff_shift_all=d[ms].tolist(), ms_share_lo_all=lo[ms].tolist(), ms_fold_all=df[ms].tolist(),
               ms_eff_bins=(A['S_sel'][ms] / A['S_all'][ms]).tolist(), ms_spec=(A['S_all'][ms] / A['S_all'][ms].sum(1, keepdims=True)).tolist(),
               re_eff_bins=(A['S_sel'][re_] / A['S_all'][re_]).tolist(),
               stress_spec=(A['S_all'][1] / A['S_all'][1].sum()).tolist(), cv_spec=(A['S_all'][0] / A['S_all'][0].sum()).tolist())
    json.dump(res, open(os.path.join(K.OUTDIR, 'mp2pi_softpion.json'), 'w'), indent=1)
    f = lambda x: f'{100 * x:+.2f}%'
    L = ['# Softest-pion dependence of the two-pion signal efficiency against the interaction-model universes', '',
         'Generated by `scripts/mp2pi_softpion.py` from the two-fold production (overlays, full exposure, combined). '
         'Efficiency: signal in the signal region over all signal. Fold: the efficiency the CV efficiency per p_soft bin '
         'gives for the spectrum of the column (spectrum-only part). Restricted: signal with both pions above 0.175 GeV/c. '
         'Sub-threshold: two charged pions, the softer below 0.10 GeV/c, predicted in the signal region. Multisim and '
         'reinteraction: mean, standard deviation, 16-84% range over the universes.', '',
         f'CV: efficiency {100 * eps[0]:.2f}%, signal-region signal {sr_sig[0]:.1f}, share of the signal below 0.175 GeV/c '
         f'{100 * lo[0]:.1f}%, sub-threshold events in the signal region {mig[0]:.1f}.', '',
         '| p_soft bin [GeV/c] | signal (CV) | efficiency (CV) | pool >= 3 tracks | multisim efficiency 16-84% | reint efficiency 16-84% | stress / CV signal |',
         '|---|---|---|---|---|---|---|']
    msb = np.array(res['ms_eff_bins']); reb = np.array(res['re_eff_bins'])
    for b in range(len(EDGES) - 1):
        L.append(f'| {EDGES[b]:g}-{EDGES[b + 1]:g} | {A["S_all"][0][b]:.1f} | {100 * eb[b]:.1f}% | {100 * A["pool3"][b] / A["S_all"][0][b]:.0f}% | '
                 f'{100 * np.percentile(msb[:, b], 16):.1f}-{100 * np.percentile(msb[:, b], 84):.1f}% | '
                 f'{100 * np.percentile(reb[:, b], 16):.1f}-{100 * np.percentile(reb[:, b], 84):.1f}% | {A["S_all"][1][b] / A["S_all"][0][b]:.2f} |')
    L += ['', '| variation | efficiency shift | spectrum-only (fold) | share below 0.175 | restricted (>0.175) efficiency shift | sub-threshold in SR |',
          '|---|---|---|---|---|---|',
          f'| stress (signal only) | {f(d[1])} | {f(df[1])} | {100 * lo[1]:.1f}% | {f(dr[1])} | {mig[1]:.1f} (extended {mig[2]:.1f}) |']
    for nm, idx in (('GENIE multisim (sd)', ms), ('reinteraction (sd)', re_)):
        L.append(f'| {nm} | {100 * np.std(d[idx]):.2f}% | {100 * np.std(df[idx]):.2f}% | {100 * np.std(lo[idx]):.1f}% | {100 * np.std(dr[idx]):.2f}% | {np.std(mig[idx]):.1f} |')
    for u in UNISIM:
        i0, i1 = ix[f'{u}_0'], ix[f'{u}_1']
        L.append(f'| {u} | {f(d[i0])} / {f(d[i1])} | {f(df[i0])} / {f(df[i1])} | {100 * lo[i0]:.1f}% / {100 * lo[i1]:.1f}% | {f(dr[i0])} / {f(dr[i1])} | {mig[i0]:.1f} / {mig[i1]:.1f} |')
    for u in SCC:
        idx = [ix[f'{u}_{i}'] for i in range(10)]
        L.append(f'| {u} (sd) | {100 * np.std(d[idx]):.2f}% | {100 * np.std(df[idx]):.2f}% | {100 * np.std(lo[idx]):.1f}% | {100 * np.std(dr[idx]):.2f}% | {np.std(mig[idx]):.1f} |')
    L += ['', f'Multisim: largest efficiency shift {f(d[ms].min())} / {f(d[ms].max())}; share below 0.175 GeV/c '
              f'{100 * lo[ms].min():.1f}-{100 * lo[ms].max():.1f}% (stress {100 * lo[1]:.1f}%). Stress shift in multisim '
              f'standard deviations: {d[1] / np.std(d[ms]):+.1f}; restricted {dr[1] / np.std(dr[ms]):+.1f}.']
    # generators
    gen = {}
    for g in GENS:
        fs = {fl: f'{OUT}gen_{g}_{fl}.json' for fl in ('numu', 'numubar')}
        if not all(os.path.exists(x) for x in fs.values()): continue
        c = {fl: np.array(json.load(open(x))['counts']) for fl, x in fs.items()}
        shp = sum(fr * c[fl][2:] / c[fl][2:].sum() for fl, fr in (('numu', 1 - F_BAR), ('numubar', F_BAR)))
        sub = sum(fr * c[fl][:2].sum() / c[fl][2:].sum() for fl, fr in (('numu', 1 - F_BAR), ('numubar', F_BAR)))
        gen[g] = dict(spec=shp.tolist(), share_lo=float(shp[:3].sum()), fold=float((eb * shp).sum()), sub_over_signal=float(sub),
                      share_lo_numu=float(c['numu'][2:5].sum() / c['numu'][2:].sum()), share_lo_numubar=float(c['numubar'][2:5].sum() / c['numubar'][2:].sum()))
    if 'genie' in gen:
        g0 = gen['genie']
        L += ['', '## Generators (spectrum folded with the CV efficiency per bin)', '',
              f'numu and numubar combined with the antineutrino share of the CV signal ({F_BAR}). Against GENIE: the model '
              'difference (standalone GENIE is the CV model). Sub-threshold: events with the softer pion below 0.10 GeV/c per '
              f'signal event, relative to GENIE, applied to the {mig[0]:.1f} predicted in the signal region (in % of its signal).', '',
              '| generator | share below 0.175 (numu, numubar) | efficiency shift vs CV | vs GENIE | sub-threshold vs GENIE | sub-threshold change, % of SR signal |',
              '|---|---|---|---|---|---|']
        for g, r in gen.items():
            dm = mig[0] * (r['sub_over_signal'] / g0['sub_over_signal'] - 1)
            r.update(shift_cv=r['fold'] / eps[0] - 1, shift_genie=r['fold'] / g0['fold'] - 1,
                     sub_ratio=r['sub_over_signal'] / g0['sub_over_signal'], sub_change_pct=100 * dm / sr_sig[0])
            L.append(f"| {g} | {100 * r['share_lo']:.1f}% ({100 * r['share_lo_numu']:.1f}, {100 * r['share_lo_numubar']:.1f}) | "
                     f"{f(r['shift_cv'])} | {f(r['shift_genie'])} | {r['sub_ratio']:.2f} | {r['sub_change_pct']:+.1f}% |")
        L.append(f"| stress | {100 * lo[1]:.1f}% | {f(df[1])} | | | |")
    res['generators'] = gen
    json.dump(res, open(os.path.join(K.OUTDIR, 'mp2pi_softpion.json'), 'w'), indent=1)
    open(os.path.join(K.OUTDIR, 'mp2pi_softpion.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    figure(res)


def gen_ratio(fbar):
    """per generator: its p_soft shape over that of standalone GENIE, numu and numubar mixed with fbar"""
    shp = {}
    for g in GENS:
        c = {fl: np.array(json.load(open(f'{OUT}gen_{g}_{fl}.json'))['counts'])[2:] for fl in ('numu', 'numubar')}
        shp[g] = (1 - fbar) * c['numu'] / c['numu'].sum() + fbar * c['numubar'] / c['numubar'].sum()
    return {g: shp[g] / shp['genie'] for g in GENS if g != 'genie'}


def term():
    smp = mc_samples(); P = [dict(np.load(f'{OUT}part_{k:02d}.npz')) for k in range(len(smp))]
    horn = ['rhc' if '_rhc' in os.path.basename(x[0]) else 'fhc' for x in smp]
    out = dict(total={}, differential={})
    L = ['# Soft-pion efficiency term: generator envelope', '',
         'Generated by `scripts/mp2pi_softpion.py term` (see its docstring). Shift of the signal efficiency when the CV '
         'signal is reweighted by the softest-pion spectrum of each generator over that of standalone GENIE; the term '
         'is the largest absolute shift.', '', '| configuration | GiBUU | NuWro | NEUT | term |', '|---|---|---|---|---|']
    for cfg, hs in (('fhc5', ['fhc']), ('rhcfull', ['rhc']), ('comb', ['fhc', 'rhc'])):
        S = sum(p['S_all'][0] for p, h in zip(P, horn) if h in hs); Ss = sum(p['S_sel'][0] for p, h in zip(P, horn) if h in hs)
        e0 = Ss.sum() / S.sum(); R = gen_ratio(F_BAR_CFG[cfg])
        sh = {g: float((Ss * r).sum() / (S * r).sum() / e0 - 1) for g, r in R.items()}
        out['total'][cfg] = dict(shifts=sh, term=max(abs(x) for x in sh.values()))
        L.append(f'| {cfg} | ' + ' | '.join(f'{100 * sh[g]:+.2f}%' for g in ('gibuu', 'nuwro', 'neut')) + f" | {100 * out['total'][cfg]['term']:.2f}% |")
    D = np.load(K.CACHE); sig = D['sig'].astype(bool); w = D['w']; sel = D['region'] == 0; ps = D['t_ppi_lo']
    b = np.clip(np.digitize(ps, EDGES) - 1, 0, len(EDGES) - 2); R = gen_ratio(F_BAR_CFG['comb'])
    edges = json.load(open(os.path.join(K.OUTDIR, 'mp2pi_kinematics_2f.json')))
    L += ['', 'Differential variables, combined, per true bin (first bin open below, last open above):', '',
          '| variable | per-bin term | integral |', '|---|---|---|']
    for v in ('pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'np'):
        e = edges[v]['edges']; t = D['t_' + v]; n = len(e) - 1
        ib = np.clip(np.digitize(t, e) - 1, 0, n - 1)
        terms = []
        for i in range(n):
            m = sig & (ib == i); e0 = w[m & sel].sum() / w[m].sum()
            terms.append(max(abs((w[m & sel] * r[b[m & sel]]).sum() / (w[m] * r[b[m]]).sum() / e0 - 1) for r in R.values()))
        integ = max(abs((w[sig & sel] * r[b[sig & sel]]).sum() / (w[sig] * r[b[sig]]).sum() / (w[sig & sel].sum() / w[sig].sum()) - 1)
                    for r in R.values())
        out['differential'][v] = dict(per_bin=[float(x) for x in terms], integral=float(integ))
        L.append(f'| {v} | ' + ', '.join(f'{100 * x:.1f}' for x in terms) + f' % | {100 * integ:.2f}% |')
    json.dump(out, open(os.path.join(K.OUTDIR, 'mp2pi_softpion_term.json'), 'w'), indent=1)
    open(os.path.join(K.OUTDIR, 'mp2pi_softpion_term.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


def read_tab(path):
    rows = [l.split() for l in open(path) if l.strip() and l.split()[0][0].isdigit()]
    if len(rows[0]) == 2: return np.array([float(r[1]) for r in rows])
    n = max(int(r[0]) for r in rows) + 1; m = np.zeros((n, n))
    for r in rows: m[int(r[0]), int(r[1])] = float(r[2])
    return m


def cov():
    U = '/data/uboone/processed/mp2pi_2f/unfold/'
    D = np.load(K.CACHE); sig = D['sig'].astype(bool); w = D['w']; sel = D['region'] == 0; ps = D['t_ppi_lo']
    b = np.clip(np.digitize(ps, EDGES) - 1, 0, len(EDGES) - 2); R = gen_ratio(F_BAR_CFG['comb'])
    edges = json.load(open(os.path.join(K.OUTDIR, 'mp2pi_kinematics_2f.json')))
    out = {}
    L = ['# Soft-pion term in the covariance of the differential results', '',
         'Generated by `scripts/mp2pi_softpion.py cov` (see its docstring). Combined configuration. Shift of the unfolded '
         'result per bin (in % of the result) for each generator; the envelope generator (largest shift of the integral) '
         'gives the covariance; "check": A_C t against the tune of the extraction, largest relative difference.', '',
         '| variable | constraint | GiBUU | NEUT | NuWro | envelope (per bin) | integral | bins set by NEUT | check |',
         '|---|---|---|---|---|---|---|---|---|']
    for v in ('pmu', 'cmu', 'cfwd', 'cbwd', 'cpipi', 'cmupi', 'np'):
        e = edges[v]['edges']; n = len(e) - 1
        ib = np.clip(np.digitize(D['t_' + v], e) - 1, 0, n - 1)
        N = np.array([w[sig & (ib == i)].sum() for i in range(n)])
        delta = {}
        for g, r in R.items():
            dl = []
            for i in range(n):
                m = sig & (ib == i); e0 = w[m & sel].sum() / w[m].sum()
                dl.append((w[m & sel] * r[b[m & sel]]).sum() / (w[m] * r[b[m]]).sum() / e0 - 1)
            delta[g] = np.array(dl)
        for mode, d in (('none', 'comb_none'), ('sidebands', 'comb')):
            uo = f'{U}{v}/{d}/unfold_output/'
            A = read_tab(uo + 'mat_table_add_smear.txt'); u = read_tab(uo + 'vec_table_uBTune.txt')[:n]
            x = read_tab(uo + 'vec_table_unfolded_signal.txt')[:n]; A = A[:n, :n]
            c = u.sum() / (A @ N).sum(); t = c * N; check = np.abs(A @ t / u - 1).max()
            dx = {g: A @ (delta[g] * t) for g in delta}
            G = np.array([dx[g] for g in dx]); arg = np.abs(G).argmax(0); names = list(dx)
            env = G[arg, np.arange(n)]
            C = np.outer(env, env)
            np.savetxt(f'{U}{v}/softpion_cov_{mode}.txt', C, header=f'soft-pion covariance (generator envelope), {v}, combined, {mode}; bin-integrated units of vec_table_unfolded_signal')
            out[f'{v} {mode}'] = dict(shift_pct={g: (100 * dx[g] / u).tolist() for g in dx}, envelope_pct=(100 * env / u).tolist(),
                                      envelope_gen=[names[a] for a in arg], integral_pct=float(100 * abs(env.sum()) / u.sum()),
                                      check=float(check))
            f = lambda y: ', '.join(f'{100 * z:+.1f}' for z in y / u)
            L.append(f"| {v} | {mode} | {f(dx['gibuu'])} | {f(dx['neut'])} | {f(dx['nuwro'])} | {f(env)} | "
                     f"{out[f'{v} {mode}']['integral_pct']:.1f}% | {sum(names[a] == 'neut' for a in arg)} | {check:.3f} |")
    json.dump(out, open(os.path.join(K.OUTDIR, 'mp2pi_softpion_cov.json'), 'w'), indent=1)
    open(os.path.join(K.OUTDIR, 'mp2pi_softpion_cov.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


def figure(res):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    e = np.array(EDGES); e[-1] = 0.8; c = 0.5 * (e[1:] + e[:-1])
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    wd = np.diff(e); ms = np.array(res['ms_spec']) / wd
    ax[0].fill_between(c, np.percentile(ms, 2.5, 0), np.percentile(ms, 97.5, 0), step='mid', color='#9ecae1', label='GENIE multisim 95%')
    ax[0].stairs(np.array(res['cv_spec']) / wd, e, color='k', lw=1.5, label='CV')
    ax[0].stairs(np.array(res['stress_spec']) / wd, e, color='#c8553d', lw=1.5, ls='--', label='stress')
    for g, col in (('gibuu', '#2a9d58'), ('nuwro', '#7b4fa0'), ('neut', '#e9a23b')):
        if g in res.get('generators', {}): ax[0].stairs(np.array(res['generators'][g]['spec']) / wd, e, color=col, lw=1.2, label={'gibuu': 'GiBUU', 'nuwro': 'NuWro', 'neut': 'NEUT'}[g])
    ax[0].set_xlabel(r'true $p_{\pi}^{\rm soft}$ [GeV/$c$]'); ax[0].set_ylabel(r'signal share per GeV/$c$'); ax[0].legend(fontsize=8)
    mb = np.array(res['ms_eff_bins']); rb = np.array(res['re_eff_bins'])
    ax[1].fill_between(c, 100 * np.percentile(mb, 16, 0), 100 * np.percentile(mb, 84, 0), step='mid', color='#9ecae1', label='GENIE multisim 68%')
    ax[1].fill_between(c, 100 * np.percentile(rb, 16, 0), 100 * np.percentile(rb, 84, 0), step='mid', color='#fdae6b', alpha=0.6, label='reinteraction 68%')
    ax[1].stairs(100 * np.array(res['eff_cv']), e, color='k', lw=1.5, label='CV')
    ax[1].set_xlabel(r'true $p_{\pi}^{\rm soft}$ [GeV/$c$]'); ax[1].set_ylabel('signal efficiency [%]'); ax[1].legend(fontsize=8)
    ax[2].hist(100 * np.array(res['ms_eff_shift_all']), bins=40, color='#9ecae1', label='GENIE multisim')
    ax[2].axvline(100 * res['stress']['eff_shift'], color='#c8553d', ls='--', label='stress')
    for g, col in (('gibuu', '#2a9d58'), ('nuwro', '#7b4fa0'), ('neut', '#e9a23b')):
        if g in res.get('generators', {}): ax[2].axvline(100 * res['generators'][g]['shift_genie'], color=col, lw=1.2,
                                                        label={'gibuu': 'GiBUU', 'nuwro': 'NuWro', 'neut': 'NEUT'}[g] + ' (vs GENIE)')
    ax[2].set_xlabel('efficiency shift [%]'); ax[2].set_ylabel('universes'); ax[2].legend(fontsize=8)
    fig.text(0.01, 0.99, 'MicroBooNE NuMI simulation, internal; combined FHC+RHC, two-fold selection', fontsize=8, color='0.45', va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(K.OUTDIR, '..', '..', 'figures', 'multipion', 'softpion.pdf')); plt.close(fig)


if __name__ == '__main__':
    {'build1': lambda: build1(sys.argv[2]), 'report': report, 'term': term, 'cov': cov, 'n': lambda: print(len(mc_samples()))}[sys.argv[1]]()
