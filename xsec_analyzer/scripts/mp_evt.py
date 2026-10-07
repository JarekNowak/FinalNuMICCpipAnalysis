#!/usr/bin/env python3
"""mp_evt.py -- Phase 2, items 1-4, of report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md: assignment of the
tracks to particle roles, event features, the two-pion event classifier and its working point.

Input: /data/uboone/processed/mp_evt (slurm/slurm_mp_evt.sbatch), every overlay, beam-off and dirt sample of
the release list processed with CC1mu2piEvt (the event-level track dump: class probabilities and kinematics of
every generation-2 track with track score >= 0.5 and valid calorimetry, in events with one neutrino slice and
the reco vertex in the fiducial volume) next to CC1mu2pi, CC1mu2piNewPID, CC1mu3pi, CC1mu3piNewPID and
CC1mu1piXp. Normalisation as the framework (scripts/mp_phase0_baseline.py): MC and dirt per run by data POT
over the distinct summed_pot, beam-off by 0.98 x triggers over gates, from the release file list.

Assignment (item 1). With the class probabilities P_mu, P_pi, P_p, P_other of each track, the hypothesis
"one muon, N pions, the rest protons or other" has the log-likelihood
    L = sum_rest log(1 - P_mu - P_pi) + log P_mu(muon) + sum_pions log P_pi(pion),
maximised over the choice of the muon and the N pions (all choices enumerated, at most 8 tracks per event).
The best and the runner-up choice are kept for N = 1, 2, 3.

Event features (item 2): the probabilities of the assigned tracks, L_2 - L_1, L_2 - L_3 and the margin over
the runner-up; class sums and counts over all tracks; track and shower counts, topological score and
CosmicIP; and, in the "full" set only, lengths, vertex distances, containment and the opening angles of the
assigned tracks.

Event classifier (item 3): XGBoost, signal = CC1mu2pi_MC_Signal, background = every other overlay event,
dirt and beam-off, with central-value weights at data exposure (classes balanced in the training). Trained on
the run periods the particle classifier was trained on (FHC Runs 1, 2, 4; RHC Runs 1, 2, 4), tested on the
periods held out of both (FHC Run 5, RHC Run 3). The dirt sample of RHC Run 3b also serves RHC Runs 1 and 2;
it is left out of the training so that no test event was seen. Domain: one slice, vertex in the fiducial
volume, software trigger, at least N + 1 dump tracks, and the reco counterpart of the truth cut of the signal
(decision D6): angle between the muon and the longest assigned pion below 2.6 rad.

Working point (item 4): the expected total uncertainty of the one-bin total, approximated as data statistics
in quadrature with the flux term 17% x (1 + B/S) of the single-pion analysis, with the test-period yields
projected to the full exposure of each horn mode (scaled by data POT).

    python3 scripts/mp_evt.py build            -> /data/uboone/processed/mp_evt/events.npz, sums.json
    python3 scripts/mp_evt.py train [--npi 3]  -> report/planning/multipion/phase2_event[_3pi].{json,md}
"""
import itertools, json, os, sys
import numpy as np
import uproot
import awkward as ak

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_phase0_baseline as B

# MP_EVT_DIR selects the processing of a particle-classifier variant (detector-robustness study)
B.NEW = os.environ.get('MP_EVT_DIR', '/data/uboone/processed/mp_evt').rstrip('/') + '/'
NEW = B.NEW
OUTDIR = B.OUTDIR
NMAX = 8                                  # tracks entering the enumeration (largest max(P_mu, P_pi) first)
FLOOR = 1e-6
TRAIN_RUNS, TEST_RUNS = (1, 2, 4, 11, 12, 14), (5, 13)
FHC_RUNS, RHC_RUNS = B.FHC_RUNS, B.RHC_RUNS
TEST_DIRT_ONLY = ('neutrinoselection_filt_run3b_dirt_overlay',)   # dirt file shared with a test period
COS_D6 = np.cos(2.6)
VTAG = ('_' + sys.argv[sys.argv.index('--tag') + 1]) if '--tag' in sys.argv else ''     # variant of the event classifier
DROP = set(sys.argv[sys.argv.index('--drop') + 1].split(',')) if '--drop' in sys.argv else set()
EAUG = float(sys.argv[sys.argv.index('--aug') + 1]) if '--aug' in sys.argv else 0.   # detector-variation events in training
ODD = '--odd' in sys.argv           # detector term on the odd event numbers only (the even ones may be in a training)
PEDGES = [0.10, 0.175, 0.25, 0.35, 0.50, 0.75, 10.]    # true pion momentum bins of the efficiency check
FLUX = 0.17
KIND = {'mc': 0, 'ext': 1, 'dirt': 2}
TOPO = ['CC2pi signal', 'CC0pi', 'CC1pi', 'CC2pi', 'CC3pi+', 'CCpi0', 'CCother', 'NC', 'nueCC', 'OOFV']
W = ['spline_weight', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
E = 'CC1mu2piEvt_'
TRK = ['idx', 'contained', 'pid_mu', 'pid_pi', 'pid_p', 'pid_other', 'len', 'ts', 'dist', 'llr', 'dirx', 'diry',
       'dirz', 'true_pdg', 'true_mom']
EVT = ['evt_presel', 'evt_ntrack', 'evt_nshower', 'evt_mu_cand_idx', 'cutflow_bits', 'MC_Signal', 'MC_Signal_1p',
       'mc_topology', 'Selected', 'mc_pi_mom']
OTHER = ['CC1mu2piNewPID_Selected', 'CC1mu3pi_Selected', 'CC1mu3pi_MC_Signal', 'CC1mu3piNewPID_Selected',
         'CC1mu1piXp_Selected', 'CC1mu1piXp_MC_Signal', 'CC1mu1piXp_sb_multipi']
# comparison selections: name -> (selected branch, signal branch)
COMPARE = {'CC1mu2pi': (E + 'Selected', E + 'MC_Signal'), 'CC1mu2piNewPID': ('CC1mu2piNewPID_Selected', E + 'MC_Signal'),
           'CC1mu3pi': ('CC1mu3pi_Selected', 'CC1mu3pi_MC_Signal'),
           'CC1mu3piNewPID': ('CC1mu3piNewPID_Selected', 'CC1mu3pi_MC_Signal')}


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


# ---------------------------------------------------------------- assignment
def combos(npi):
    """All (muon, pions) choices among NMAX slots: arrays (n_combos,) and (n_combos, npi)."""
    mu, pis = [], []
    for m in range(NMAX):
        rest = [i for i in range(NMAX) if i != m]
        for c in itertools.combinations(rest, npi):
            mu.append(m); pis.append(c)
    return np.array(mu), np.array(pis)


COMBOS = {n: combos(n) for n in (1, 2, 3)}


def pad(x, fill):
    return ak.to_numpy(ak.fill_none(ak.pad_none(x, NMAX, clip=True), fill)).astype(np.float64)


def assign(P, npi):
    """P: dict of padded (n_events, NMAX) arrays. Returns best and runner-up L and the best choice."""
    pm, pp = np.clip(P['mu'], FLOOR, 1.), np.clip(P['pi'], FLOOR, 1.)
    rest = np.log(np.clip(1. - P['mu'] - P['pi'], FLOOR, 1.))
    valid = P['valid']
    gmu = np.where(valid, np.log(pm) - rest, -np.inf)
    gpi = np.where(valid, np.log(pp) - rest, -np.inf)
    base = np.where(valid, rest, 0.).sum(axis=1)
    M, PI = COMBOS[npi]
    s = gmu[:, M] + gpi[:, PI].sum(axis=2)                  # (events, combos)
    s = np.where(np.isfinite(s), s, -np.inf)
    order = np.argsort(-s, axis=1)[:, :2]
    best = np.take_along_axis(s, order[:, :1], 1)[:, 0]
    second = np.take_along_axis(s, order[:, 1:2], 1)[:, 0]
    ib = order[:, 0]
    return base + best, base + second, M[ib], PI[ib]


def features(P, T, ev, npi):
    """Event features of the best npi-pion assignment. P, T: padded per-track arrays; ev: event arrays.
    Rows with fewer than npi + 1 tracks give -inf scores; they are masked by the caller."""
    with np.errstate(invalid='ignore'):
        return _features(P, T, ev, npi)


def _features(P, T, ev, npi):
    n = P['valid'].sum(axis=1)
    L = {}
    for k in (1, 2, 3):
        b, s, m, pis = assign(P, k)
        L[k] = (b, s, m, pis)
    lb, ls, mu, pis = L[npi]
    rows = np.arange(len(n))
    # order the assigned pions by track length, longest first (decision D6 refers to the longest)
    plen = np.take_along_axis(T['len'], pis, 1)
    pis = np.take_along_axis(pis, np.argsort(-plen, axis=1), 1)
    f = {}
    big = 50.
    f['dL_up'] = lb - L[npi - 1][0]                                   # against one pion fewer (npi >= 2)
    f['dL_down'] = (np.where((n >= npi + 2) & np.isfinite(L[npi + 1][0]), lb - L[npi + 1][0], big)
                    if npi < 3 else np.full(len(n), big))                 # against one pion more
    f['dL_runnerup'] = np.where(np.isfinite(ls), lb - ls, big)
    f['L_best'] = lb
    for cls in ('mu', 'pi', 'p', 'other'):
        f[f'mu_P{cls}'] = P[cls][rows, mu]
    f['mu_is_cand'] = (T['idx'][rows, mu] == ev['evt_mu_cand_idx']).astype(float)
    for k in range(npi):
        for cls in ('mu', 'pi', 'p'):
            f[f'pi{k}_P{cls}'] = P[cls][rows, pis[:, k]]
    # tracks not assigned: largest proton and pion probability, number of proton-like ones
    roles = np.zeros_like(P['valid'])
    roles[rows, mu] = True
    for k in range(npi):
        roles[rows, pis[:, k]] = True
    free = P['valid'] & ~roles
    f['rest_n'] = free.sum(axis=1).astype(float)
    f['rest_maxPp'] = np.where(free, P['p'], -1.).max(axis=1)
    f['rest_maxPpi'] = np.where(free, P['pi'], -1.).max(axis=1)
    f['rest_nPp05'] = (free & (P['p'] > 0.5)).sum(axis=1).astype(float)
    for cls in ('mu', 'pi', 'p', 'other'):
        f[f'sum_P{cls}'] = np.where(P['valid'], P[cls], 0.).sum(axis=1)
        f[f'n_P{cls}05'] = (P['valid'] & (P[cls] > 0.5)).sum(axis=1).astype(float)
    f['n_dump'] = n.astype(float)
    f['n_short'] = (P['valid'] & (T['len'] < 5.)).sum(axis=1).astype(float)
    f['evt_ntrack'] = ev['evt_ntrack'].astype(float)
    f['evt_nshower'] = ev['evt_nshower'].astype(float)
    f['topo'] = ev['topological_score'].astype(float)
    f['cosmicip'] = np.clip(ev['CosmicIP'].astype(float), -1., 1000.)
    # kinematic features ("full" set)
    kin = {}
    def trk(name, j): return T[name][rows, j]
    for tag, j in [('mu', mu)] + [(f'pi{k}', pis[:, k]) for k in range(npi)]:
        kin[f'{tag}_len'] = trk('len', j)
        kin[f'{tag}_dist'] = trk('dist', j)
        kin[f'{tag}_endcont'] = (trk('contained', j).astype(int) >> 1 & 1).astype(float)
        kin[f'{tag}_ts'] = trk('ts', j)
    def cosang(a, b):
        return (T['dirx'][rows, a] * T['dirx'][rows, b] + T['diry'][rows, a] * T['diry'][rows, b]
                + T['dirz'][rows, a] * T['dirz'][rows, b])
    for k in range(npi):
        kin[f'cos_mu_pi{k}'] = cosang(mu, pis[:, k])
    for a, b in itertools.combinations(range(npi), 2):
        kin[f'cos_pi{a}_pi{b}'] = cosang(pis[:, a], pis[:, b])
    # truth of the assigned tracks (diagnostics only, never a feature)
    truth = {'mu_true_pdg': trk('true_pdg', mu)}
    for k in range(npi):
        truth[f'pi{k}_true_pdg'] = trk('true_pdg', pis[:, k])
    d6 = kin['cos_mu_pi0'] > COS_D6
    return f, kin, truth, d6


def padded(a):
    """Per-track arrays padded to NMAX; tracks beyond NMAX dropped by largest max(P_mu, P_pi)."""
    key = np.maximum(a[E + 'etrk_pid_mu'], a[E + 'etrk_pid_pi'])
    order = ak.argsort(key, axis=1, ascending=False)
    P = {'valid': ak.to_numpy(ak.fill_none(ak.pad_none(ak.ones_like(key[order], dtype=bool), NMAX, clip=True), False))}
    for cls, br in (('mu', 'pid_mu'), ('pi', 'pid_pi'), ('p', 'pid_p'), ('other', 'pid_other')):
        P[cls] = pad(a[E + 'etrk_' + br][order], 0.)
    T = {x: pad(a[E + 'etrk_' + x][order], -1.) for x in TRK if x not in ('pid_mu', 'pid_pi', 'pid_p', 'pid_other')}
    return P, T


# ---------------------------------------------------------------- build
def scan_file(path, kind, run, scale, key, file_id, add, sums):
    """Rows of one processed file (preselected, software trigger, at least three dump tracks) with their
    features, and the per-sample sums (generated signal, comparison selections) under sums[...][key]."""
    tag = os.path.basename(path)[len('xsec-ana-'):-len('.root')]
    t = uproot.open(path)['stv_tree']
    br = [E + 'etrk_' + x for x in TRK] + [E + x for x in EVT] + OTHER + ['topological_score', 'CosmicIP', 'nslice']
    if kind != 'ext': br += W
    br = [b for b in br if b in t.keys()]
    print(f'  run {run:2d} {kind:4s} {tag} scale {scale:.5g} entries {t.num_entries}', flush=True)
    for a, rep in t.iterate(br, step_size=200000, library='ak', report=True):
        ne = len(a)
        w = np.ones(ne) if kind == 'ext' else safe(np.prod([ak.to_numpy(a[x]).astype(float) for x in W], axis=0))
        ww = w * scale
        sig2 = ak.to_numpy(a[E + 'MC_Signal']).astype(bool) if kind == 'mc' else np.zeros(ne, bool)
        sig3 = ak.to_numpy(a['CC1mu3pi_MC_Signal']).astype(bool) if kind == 'mc' else np.zeros(ne, bool)
        g = sums['generated'].setdefault(key, {'2pi': 0., '3pi': 0., '2pi_1p': 0.,
                                               'leading': [0.] * (len(PEDGES) - 1), 'subleading': [0.] * (len(PEDGES) - 1)})
        g['2pi'] += ww[sig2].sum(); g['3pi'] += ww[sig3].sum()
        if kind == 'mc':
            g['2pi_1p'] += ww[sig2 & ak.to_numpy(a[E + 'MC_Signal_1p']).astype(bool)].sum()
            mpm = a[E + 'mc_pi_mom'][sig2]
            for k2, x in (('leading', ak.firsts(mpm)), ('subleading', ak.firsts(mpm[:, 1:]))):
                h = np.histogram(ak.to_numpy(ak.fill_none(x, -1.)), bins=PEDGES, weights=ww[sig2])[0]
                g[k2] = [u + v for u, v in zip(g[k2], h)]
        for name, (sb, sg) in COMPARE.items():
            sel = ak.to_numpy(a[sb]).astype(bool)
            sig = ak.to_numpy(a[sg]).astype(bool) if kind == 'mc' else np.zeros(ne, bool)
            c = sums['compare'].setdefault(name, {}).setdefault(key, {'sig': 0., 'mc_bkg': 0., 'ext': 0., 'dirt': 0.})
            c['sig'] += ww[sel & sig].sum()
            c['mc_bkg' if kind == 'mc' else kind] += ww[sel & ~sig].sum() if kind == 'mc' else ww[sel].sum()
        pre = (ak.to_numpy(a[E + 'evt_presel']).astype(bool) & (ak.to_numpy(a[E + 'cutflow_bits']) & 1).astype(bool)
               & (ak.to_numpy(ak.num(a[E + 'etrk_idx'])) >= 3))
        if not pre.any(): continue
        b = a[pre]
        P, T = padded(b)
        ev = {x: ak.to_numpy(b[E + x]) for x in ('evt_ntrack', 'evt_nshower', 'evt_mu_cand_idx')}
        ev['topological_score'] = ak.to_numpy(b['topological_score']); ev['CosmicIP'] = ak.to_numpy(b['CosmicIP'])
        for npi in (2, 3):
            ok = P['valid'].sum(axis=1) >= npi + 1
            f, kin, truth, d6 = features(P, T, ev, npi)
            for k, v in list(f.items()) + list(kin.items()) + list(truth.items()):
                add(f'n{npi}_{k}', np.where(ok, v, np.nan))
            add(f'n{npi}_ok', ok); add(f'n{npi}_d6', ok & d6)
        add('run', np.full(len(b), run)); add('kind', np.full(len(b), KIND[kind]))
        add('fid', np.full(len(b), file_id)); add('entry', rep.tree_entry_start + np.nonzero(pre)[0])
        add('test_dirt_only', np.full(len(b), tag in TEST_DIRT_ONLY))
        add('w', ww[pre])
        add('sig2', sig2[pre]); add('sig3', sig3[pre])
        add('sig2_1p', (ak.to_numpy(b[E + 'MC_Signal_1p']).astype(bool) & sig2[pre]) if kind == 'mc' else np.zeros(len(b), bool))
        add('topology', ak.to_numpy(b[E + 'mc_topology']) if kind == 'mc' else np.full(len(b), 10))
        for name, (sb, _) in COMPARE.items():
            add('sel_' + name, ak.to_numpy(b[sb]).astype(bool))
        add('sel_1pi', ak.to_numpy(b['CC1mu1piXp_Selected']).astype(bool))
        add('sb_multipi_1pi', ak.to_numpy(b['CC1mu1piXp_sb_multipi']).astype(bool))
        mpm = b[E + 'mc_pi_mom']
        add('true_pi0_mom', ak.to_numpy(ak.fill_none(ak.firsts(mpm), -1.)))
        add('true_pi1_mom', ak.to_numpy(ak.fill_none(ak.firsts(mpm[:, 1:]), -1.)))


def build():
    on, ext, mc, dirt = B.parse_conf()
    sc, checks = B.scales(on, ext, mc, dirt)
    bad = [c for c in checks if abs(c[3] - c[4]) > 1e-5 * abs(c[4])]
    if bad: sys.exit(f'summed_pot differs from the release copies: {bad}')
    cols = {}
    sums = {'generated': {}, 'compare': {}, 'pot': {str(r): on[r][1] for r in on}, 'files': []}
    def add(name, x):
        cols.setdefault(name, []).append(np.asarray(x))
    for (run, kind), files in sorted(sc.items()):
        for path, s in files:
            sums['files'].append([run, kind, path, s])
            scan_file(path, kind, run, s, str(run), len(sums['files']) - 1, add, sums)
    out = {k: np.concatenate(v) for k, v in cols.items()}
    np.savez_compressed(NEW + 'events.npz', **out)
    json.dump(sums, open(NEW + 'sums.json', 'w'), indent=1)
    print(f'{len(out["w"])} rows written to {NEW}events.npz')


DETVAR = ['CV', 'LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']


def build_detvar():
    """The Run-4 detector-variation samples (FHC and RHC: CV and eight variations), each normalised to
    1e20 POT of its own exposure, into detvar_events.npz and detvar_sums.json."""
    cols = {}
    sums = {'generated': {}, 'compare': {}, 'files': []}
    def add(name, x):
        cols.setdefault(name, []).append(np.asarray(x))
    for horn, run in (('fhc', 4), ('rhc', 14)):
        for v in DETVAR:
            path = f'{NEW}detvar/xsec-ana-detvar_run4{horn}_{v}.root'
            scale = 1e20 / B.summed_pot(path)
            sums['files'].append([run, v, path, scale])
            scan_file(path, 'mc', run, scale, f'{horn}_{v}', len(sums['files']) - 1, add, sums)
    out = {k: np.concatenate(v) for k, v in cols.items()}
    np.savez_compressed(NEW + 'detvar_events.npz', **out)
    json.dump(sums, open(NEW + 'detvar_sums.json', 'w'), indent=1)
    print(f'{len(out["w"])} detector-variation rows written')


def detvar_evt():
    """Event number of every detector-variation row (detvar_events.npz), so that a split by event number puts the same
    event of the CV and of each variation on the same side -> detvar_evt.npy. Taken from the track dumps of the same raw
    inputs (/data/uboone/processed/mp_pid_detvar, macros/mp_pid/dump_pid_tracks.C: run, subrun, event and raw entry of
    every event with a pool track, which every row has); the processed FHC files were made from part 1 followed by part
    2, so part-2 entries are offset by the entries of part 1 (read from the raw tree header only)."""
    V = np.load(NEW + 'detvar_events.npz')
    fid_col, ent_col = V['fid'], V['entry']
    vs = json.load(open(NEW + 'detvar_sums.json'))
    man = {}
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'slurm', 'mp_evt_detvar_manifest.list')):
        raw, ft, sel, out = line.strip().split('|')
        man[os.path.basename(out)] = raw.split('+')
    evt = np.full(len(fid_col), -1, dtype=np.int64)
    ids = np.full((len(fid_col), 2), -1, dtype=np.int64)      # run and subrun, for the event identity of the bootstrap
    for fid, (run, var, path, scale) in enumerate(vs['files']):
        raws = man[os.path.basename(path)]
        stem = os.path.basename(path)[len('xsec-ana-'):-len('.root')]
        dumps = [f'/data/uboone/processed/mp_pid_detvar/pidtrk-{stem}' + (f'_part{k + 1}' if len(raws) > 1 else '') + '.root'
                 for k in range(len(raws))]
        ents, evts, rs, off = [], [], [], 0
        for k, (r, dp) in enumerate(zip(raws, dumps)):
            a = uproot.open(dp)['trk'].arrays(['entry', 'evt', 'run', 'sub'], library='np')
            ents.append(a['entry'].astype(np.int64) + off); evts.append(a['evt'].astype(np.int64))
            rs.append(np.stack([a['run'].astype(np.int64), a['sub'].astype(np.int64)], axis=1))
            if k + 1 < len(raws):                 # entries of this part, from the tree header (ROOT reads it fast)
                import ROOT
                f_ = ROOT.TFile.Open(r); off += int(f_.Get('nuselection/NeutrinoSelectionFilter').GetEntries()); f_.Close()
        ent_all, evt_all, rs_all = np.concatenate(ents), np.concatenate(evts), np.concatenate(rs)
        ue, ui = np.unique(ent_all, return_index=True)
        m = np.flatnonzero(fid_col == fid)
        pos = np.searchsorted(ue, ent_col[m]); pos = np.clip(pos, 0, len(ue) - 1)
        found = ue[pos] == ent_col[m]
        evt[m] = np.where(found, evt_all[ui][pos], -1)
        ids[m] = np.where(found[:, None], rs_all[ui][pos], -1)
        print(f'  {stem}: {len(m)} rows, {int((evt[m] < 0).sum())} without an event number', flush=True)
    np.save(NEW + 'detvar_evt.npy', evt)
    np.save(NEW + 'detvar_runsub.npy', ids)


# ---------------------------------------------------------------- train and evaluate
FEAT_PID = ['dL_up', 'dL_down', 'dL_runnerup', 'L_best', 'mu_Pmu', 'mu_Ppi', 'mu_Pp', 'mu_Pother', 'mu_is_cand',
            'rest_n', 'rest_maxPp', 'rest_maxPpi', 'rest_nPp05', 'sum_Pmu', 'sum_Ppi', 'sum_Pp', 'sum_Pother',
            'n_Pmu05', 'n_Ppi05', 'n_Pp05', 'n_Pother05', 'n_dump', 'n_short', 'evt_ntrack', 'evt_nshower',
            'topo', 'cosmicip']


def feat_names(npi, full):
    f = FEAT_PID + [f'pi{k}_P{c}' for k in range(npi) for c in ('mu', 'pi', 'p')]
    if full:
        f += [f'{t}_{x}' for t in ['mu'] + [f'pi{k}' for k in range(npi)] for x in ('len', 'dist', 'endcont', 'ts')]
        f += [f'cos_mu_pi{k}' for k in range(npi)] + [f'cos_pi{a}_pi{b}' for a, b in itertools.combinations(range(npi), 2)]
    return [x for x in f if x not in DROP]


def wauc(y, s, w):
    o = np.argsort(-s)
    y, w = y[o], w[o]
    tp = np.cumsum(w * y); fp = np.cumsum(w * (1 - y))
    tp = np.concatenate([[0.], tp / tp[-1]]); fp = np.concatenate([[0.], fp / fp[-1]])
    return float(np.sum((fp[1:] - fp[:-1]) * (tp[1:] + tp[:-1]) / 2.))


def proxy(S, Bk):
    stat = np.sqrt(S + Bk) / S if S > 0 else np.inf
    flux = FLUX * (1. + Bk / S) if S > 0 else np.inf
    return stat, flux, float(np.hypot(stat, flux))


def train():
    import xgboost as xgb
    npi = 3 if '--npi' in sys.argv and sys.argv[sys.argv.index('--npi') + 1] == '3' else 2
    D = dict(np.load(NEW + 'events.npz'))
    sums = json.load(open(NEW + 'sums.json'))
    pot = {int(r): v for r, v in sums['pot'].items()}
    lab = D[f'sig{npi}'].astype(int)
    dom = D[f'n{npi}_d6'].astype(bool)
    run, kind, w = D['run'], D['kind'], D['w']
    tr_all = dom & np.isin(run, TRAIN_RUNS) & ~D['test_dirt_only']
    te = dom & np.isin(run, TEST_RUNS)
    rng = np.random.default_rng(20261006)
    val = tr_all & (rng.random(len(w)) < 0.2)
    tr = tr_all & ~val
    res = {'npi': npi, 'train_runs': TRAIN_RUNS, 'test_runs': TEST_RUNS, 'flux_coefficient': FLUX,
           'rows': {'domain': int(dom.sum()), 'train': int(tr.sum()), 'val': int(val.sum()), 'test': int(te.sum())},
           'raw_signal': {'train': int((tr & (lab == 1)).sum()), 'test': int((te & (lab == 1)).sum())}}
    # projection of the test periods to the full exposure of each horn mode
    def horn_runs(r): return FHC_RUNS if r in FHC_RUNS else RHC_RUNS
    proj = {r: sum(pot[x] for x in horn_runs(r)) / sum(pot[x] for x in TEST_RUNS if x in horn_runs(r)) for r in TEST_RUNS}
    scores = {}
    for fs in ('pid', 'full'):
        names = feat_names(npi, fs == 'full')
        X = np.stack([D[f'n{npi}_{k}'] for k in names], axis=1)
        def bal(m):
            # classes balanced, total weight = number of rows (so min_child_weight counts rows)
            ww = w[m].copy(); y = lab[m]; nr = 0.5 * len(ww)
            ww[y == 1] *= nr / ww[y == 1].sum(); ww[y == 0] *= nr / ww[y == 0].sum()
            return ww
        Xtr, ytr, wtr = X[tr], lab[tr], bal(tr)
        if EAUG > 0.:
            # the Run-4 detector-variation events (CV and the eight variations) with an even event number, classes
            # balanced within them, weighted to a fraction EAUG of the nominal training weight
            Va = np.load(NEW + 'detvar_events.npz'); ev_a = np.load(NEW + 'detvar_evt.npy')
            ma = Va[f'n{npi}_d6'].astype(bool) & (ev_a % 2 == 0)
            Xa = np.stack([Va[f'n{npi}_{k}'][ma] for k in names], axis=1); ya = Va[f'sig{npi}'][ma].astype(int)
            wa = Va['w'][ma].astype(float).copy()
            wa[ya == 1] *= 0.5 / wa[ya == 1].sum(); wa[ya == 0] *= 0.5 / wa[ya == 0].sum()
            wa *= EAUG * wtr.sum()
            Xtr = np.concatenate([Xtr, Xa]); ytr = np.concatenate([ytr, ya]); wtr = np.concatenate([wtr, wa])
            print(f'  {fs}: {int(ma.sum())} detector-variation events added (weight fraction {EAUG})', flush=True)
        dtr = xgb.DMatrix(Xtr, label=ytr, weight=wtr, feature_names=names)
        dva = xgb.DMatrix(X[val], label=lab[val], weight=bal(val), feature_names=names)
        par = dict(objective='binary:logistic', eval_metric='auc', max_depth=5, eta=0.05, subsample=0.8,
                   colsample_bytree=0.8, min_child_weight=10., tree_method='hist', seed=1)
        bst = xgb.train(par, dtr, 3000, evals=[(dva, 'val')], early_stopping_rounds=100, verbose_eval=False)
        bst.save_model(NEW + f'evt_clf_n{npi}{VTAG}_{fs}.json')
        s = bst.predict(xgb.DMatrix(X, feature_names=names), iteration_range=(0, bst.best_iteration + 1))
        scores[fs] = s
        np.save(NEW + f'scores_n{npi}{VTAG}_{fs}.npy', s.astype(np.float32))
        gain = bst.get_score(importance_type='gain')
        res[fs] = {'features': names, 'best_iteration': int(bst.best_iteration),
                   'auc_test': wauc(lab[te], s[te], w[te]), 'auc_val': wauc(lab[val], s[val], w[val]),
                   'importance_gain': dict(sorted(gain.items(), key=lambda x: -x[1])[:15])}
    def combine(per_run, runs, f=1.):
        """Per-run yields {run: (signal, background)} of a sample of the given run periods, drawn with
        fraction f, projected to the full exposure of each horn mode (by data POT)."""
        tot = {}
        for cfg, hr in (('FHC', FHC_RUNS), ('RHC', RHC_RUNS)):
            rr = [x for x in runs if x in hr]
            if not rr:
                tot[cfg] = [0., 0., 0.]; continue
            k = sum(pot[x] for x in hr) / sum(pot[x] for x in rr)
            tot[cfg] = [sum(per_run[x][0] for x in rr) * k / f, sum(per_run[x][1] for x in rr) * k / f,
                        sum(sums['generated'][str(x)][f'{npi}pi'] for x in rr) * k]
        tot['COMB'] = [tot['FHC'][q] + tot['RHC'][q] for q in range(3)]
        r = {}
        for cfg, (S, Bk, G) in tot.items():
            stat, flux, t = proxy(S, Bk)
            r[cfg] = dict(signal=S, background=Bk, efficiency=S / G if G else np.nan,
                          purity=S / (S + Bk) if S + Bk else np.nan, b_over_s=Bk / S if S else np.inf,
                          stat=stat, flux=flux, total=t)
        return r
    def rows_yields(m, runs):
        return {x: (float(w[m & (run == x) & (lab == 1)].sum()), float(w[m & (run == x) & (lab == 0)].sum())) for x in runs}
    f_val = val.sum() / tr_all.sum()
    # comparison selections from the full-sample sums (not only the classifier domain), test periods
    comp = {}
    for name in (['CC1mu2pi', 'CC1mu2piNewPID'] if npi == 2 else ['CC1mu3pi', 'CC1mu3piNewPID']):
        c = sums['compare'][name]
        comp[name] = combine({x: (c[str(x)]['sig'], c[str(x)]['mc_bkg'] + c[str(x)]['ext'] + c[str(x)]['dirt'])
                              for x in TEST_RUNS}, TEST_RUNS)
    res['compare'] = comp
    # threshold scan on the validation sample of the training periods (which chooses the threshold) and on
    # the test periods (which measure the performance)
    thresholds = np.round(np.concatenate([np.arange(0.0, 0.9, 0.01), np.arange(0.9, 0.999, 0.002)]), 3)
    scan = {}
    for fs in ('pid', 'full'):
        rows_v, rows_t = [], []
        for thr in thresholds:
            m = dom & (scores[fs] > thr)
            rows_v.append(dict(threshold=float(thr), **combine(rows_yields(m & val, TRAIN_RUNS), TRAIN_RUNS, f_val)))
            rows_t.append(dict(threshold=float(thr), **combine(rows_yields(m & te, TEST_RUNS), TEST_RUNS)))
        iv = int(np.argmin([x['COMB']['total'] for x in rows_v]))
        it = int(np.argmin([x['COMB']['total'] for x in rows_t]))
        scan[fs] = dict(validation=rows_v, test=rows_t, chosen=rows_t[iv], best_on_test=rows_t[it])
    res['scan'] = scan
    res['validation_fraction'] = float(f_val)
    # the working point: the feature set with the smaller total at its validation-chosen threshold, or the one
    # given with --set and --threshold (the diagnostics below are then quoted there)
    fs = 'full' if scan['full']['chosen']['COMB']['total'] <= scan['pid']['chosen']['COMB']['total'] else 'pid'
    thr = scan[fs]['chosen']['threshold']
    if '--threshold' in sys.argv:
        fs = sys.argv[sys.argv.index('--set') + 1] if '--set' in sys.argv else 'pid'
        thr = float(sys.argv[sys.argv.index('--threshold') + 1])
    sel = dom & (scores[fs] > thr)
    res['working_point'] = dict(feature_set=fs, threshold=thr, given='--threshold' in sys.argv,
                                **combine(rows_yields(sel & te, TEST_RUNS), TEST_RUNS))
    tem = sel & te
    comp_bkg = {}
    for i, nm in enumerate(TOPO):           # class 0 (two-pion signal) is background of the 3pi classifier
        comp_bkg[nm] = float(sum(w[tem & (lab == 0) & (kind == 0) & (D['topology'] == i) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
    comp_bkg['dirt'] = float(sum(w[tem & (kind == 2) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
    comp_bkg['beam-off'] = float(sum(w[tem & (kind == 1) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
    res['working_point']['background_by_topology_COMB'] = comp_bkg
    # identity of the assigned tracks in selected signal and background (test periods, at data exposure)
    def frac(m, col, pdg):
        tot = w[m].sum()
        return float(w[m & (np.abs(D[col]) == pdg)].sum() / tot) if tot else np.nan
    ident = {}
    for part, m in (('signal', tem & (lab == 1)), ('background', tem & (lab == 0) & (kind == 0))):
        ident[part] = {'muon is a muon': frac(m, f'n{npi}_mu_true_pdg', 13)}
        for k in range(npi):
            ident[part][f'pion {k} is a pion'] = frac(m, f'n{npi}_pi{k}_true_pdg', 211)
            ident[part][f'pion {k} is a proton'] = frac(m, f'n{npi}_pi{k}_true_pdg', 2212)
    res['working_point']['assigned_identity_test'] = ident
    # overlaps (test periods): current selection, the blind 1pi signal region, the opened multi-pi region
    def ov(mask):
        return dict(signal=float(w[tem & (lab == 1) & mask].sum() / max(w[tem & (lab == 1)].sum(), 1e-12)),
                    all=float(w[tem & mask].sum() / max(w[tem].sum(), 1e-12)))
    cur = 'sel_CC1mu2pi' if npi == 2 else 'sel_CC1mu3pi'
    res['working_point']['overlap_test'] = {'current selection': ov(D[cur].astype(bool)),
                                             'single-pion signal region (blind)': ov(D['sel_1pi'].astype(bool)),
                                             'single-pion multi-pi control region (opened)': ov(D['sb_multipi_1pi'].astype(bool))}
    # efficiency against the true pion momenta, working point and current selection (test periods)
    if npi == 2:
        effs = {}
        allsig = np.isin(run, TEST_RUNS) & (lab == 1)
        for key, col in (('leading', 'true_pi0_mom'), ('subleading', 'true_pi1_mom')):
            e_new, e_cur = [], []
            for i, (lo, hi) in enumerate(zip(PEDGES[:-1], PEDGES[1:])):
                mb = allsig & (D[col] >= lo) & (D[col] < hi)
                den = sum(sums['generated'][str(r)][key][i] for r in TEST_RUNS)
                e_new.append(float(w[mb & sel].sum() / den) if den else np.nan)
                e_cur.append(float(w[mb & D[cur].astype(bool)].sum() / den) if den else np.nan)
            effs[key] = dict(edges=PEDGES, working_point=e_new, current=e_cur)
        res['working_point']['efficiency_vs_true_pion_momentum'] = effs
    # no classifier: the events whose most likely pion count is npi (L_npi above L_npi-1 and L_npi+1)
    ml = dom & (D[f'n{npi}_dL_up'] > 0) & (D[f'n{npi}_dL_down'] > 0)
    res['assignment_only'] = combine(rows_yields(ml & te, TEST_RUNS), TEST_RUNS)
    # first look at the control sample of decision D8 (two-pion analysis only): events of the domain in
    # which three pions are more likely than two, and their overlap with the working point
    if npi == 2:
        cr = dom & (D['n2_dL_down'] < 0) & (D['n2_n_dump'] >= 4)
        crt = cr & te
        cls = {}
        for i, nm in enumerate(TOPO):
            cls[nm] = float(sum(w[crt & (kind == 0) & (D['topology'] == i) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
        cls['dirt'] = float(sum(w[crt & (kind == 2) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
        cls['beam-off'] = float(sum(w[crt & (kind == 1) & (run == r)].sum() * proj[r] for r in TEST_RUNS))
        res['three_pion_like_sample'] = dict(events_COMB=sum(cls.values()), by_class_COMB=cls,
                                             overlap_with_working_point=float(w[crt & sel].sum() / max(w[crt].sum(), 1e-12)))
    tag = '' if npi == 2 else '_3pi'
    json.dump(res, open(os.path.join(OUTDIR, f'phase2_event{tag}{VTAG}.json'), 'w'), indent=1, default=float)
    write_md(res, tag)
    print(json.dumps(res['working_point'], indent=1, default=float)[:3000])


# ---------------------------------------------------------------- systematic terms against the threshold
# weight branch -> (term, multisim?, convention), as configs/ccpi_systcalc_numi.conf and
# UniverseMaker's apply_cv_correction_weights: UBGenie weights replace the tune, the others multiply the CV
UNIV = {'weight_All_UBGenie': ('xsec_multi', True, 'genie'), 'weight_RPA_CCQE_UBGenie': ('xsec_RPA_CCQE', True, 'genie'),
        'weight_xsr_scc_Fa3_SCC': ('xsec_Fa3_SCC', True, 'cv'), 'weight_xsr_scc_Fv3_SCC': ('xsec_Fv3_SCC', True, 'cv'),
        'weight_reint_all': ('reint', True, 'cv')}
for _k in ('AxFFCCQEshape', 'DecayAngMEC', 'NormCCCOH', 'NormNCCOH', 'ThetaDelta2NRad', 'Theta_Delta2Npi',
           'VecFFCCQEshape', 'XSecShape_CCMEC'):
    UNIV[f'weight_{_k}_UBGenie'] = (f'xsec_{_k}', False, 'genie')
CMP = {2: ['CC1mu2pi', 'CC1mu2piNewPID'], 3: ['CC1mu3pi', 'CC1mu3piNewPID']}
SYST_THR = np.round(np.concatenate([np.arange(0.50, 0.95, 0.01), np.arange(0.95, 0.99, 0.0025)]), 4)


def syst():
    """Cross-section and reinteraction terms of the one-bin total against the event-score threshold, from the
    universe weights of the test periods. The extracted total in universe u relative to the CV, with the
    pseudo-data N fixed at the CV prediction, is (N - B_u) G_u / (S_u G_CV): the background B_u (beam-off
    unvaried) and the efficiency S_u / G_u both move."""
    npi = 3 if '--npi' in sys.argv and sys.argv[sys.argv.index('--npi') + 1] == '3' else 2
    D = dict(np.load(NEW + 'events.npz'))
    sums = json.load(open(NEW + 'sums.json'))
    pot = {int(r): v for r, v in sums['pot'].items()}
    lab = D[f'sig{npi}'].astype(bool); dom = D[f'n{npi}_d6'].astype(bool)
    sigbr = E + 'MC_Signal' if npi == 2 else 'CC1mu3pi_MC_Signal'
    # the classifier (rows outside its domain never pass) and the current selections as score 1 / 0
    FS = ('pid', 'full') + tuple('cmp_' + c for c in CMP[npi])
    scores = {fs: np.where(dom, np.load(NEW + f'scores_n{npi}{VTAG}_{fs}.npy'), -1.) for fs in ('pid', 'full')}
    scores.update({'cmp_' + c: D['sel_' + c].astype(float) for c in CMP[npi]})
    T = len(SYST_THR)
    def tbin(x):                 # number of thresholds the score exceeds
        return np.searchsorted(SYST_THR, x, side='left')
    test_files = [(fid, f) for fid, f in enumerate(sums['files']) if f[0] in TEST_RUNS]
    first_mc = next(f[2] for _, f in test_files if f[1] == 'mc')
    tt = uproot.open(first_mc)['stv_tree']
    nuniv = {wb: int(ak.max(ak.num(tt[wb].array(entry_stop=5000)))) for wb in UNIV if wb in tt.keys()}
    nuniv = {k: v for k, v in nuniv.items() if v > 0}
    print('universes:', nuniv, flush=True)
    z = lambda *sh: np.zeros(sh)
    cv = {r: dict(S={fs: z(T + 1) for fs in FS}, Bmc={fs: z(T + 1) for fs in FS}, Bext={fs: z(T + 1) for fs in FS},
                  W2={fs: z(T + 1) for fs in FS}, G=0.) for r in TEST_RUNS}
    acc = {(r, wb): dict(S={fs: z(T + 1, n) for fs in FS}, B={fs: z(T + 1, n) for fs in FS}, G=z(n))
           for r in TEST_RUNS for wb, n in nuniv.items()}
    for fid, (run, kind, path, scale) in test_files:
        rows = np.nonzero(D['fid'] == fid)[0]
        if kind == 'ext':
            for fs in FS:
                b = tbin(scores[fs][rows])
                np.add.at(cv[run]['Bext'][fs], b, D['w'][rows]); np.add.at(cv[run]['W2'][fs], b, D['w'][rows] ** 2)
            continue
        ent = D['entry'][rows]
        t = uproot.open(path)['stv_tree']
        keys = set(t.keys())
        br = W + [sigbr] + [wb for wb in nuniv if wb in keys]
        print(f'  run {run} {kind} {os.path.basename(path)}: {len(rows)} domain rows', flush=True)
        for a, rep in t.iterate(br, step_size=100000, library='ak', report=True):
            lo, hi = rep.tree_entry_start, rep.tree_entry_stop
            ne = hi - lo
            comp = {x: ak.to_numpy(a[x]).astype(float) for x in W}
            wcv = safe(comp['spline_weight'] * comp['tuned_cv_weight'] * comp['ppfx_cv_weight'] * comp['normalisation_weight']) * scale
            sig = ak.to_numpy(a[sigbr]).astype(bool) if kind == 'mc' else np.zeros(ne, bool)
            inr = (ent >= lo) & (ent < hi)
            rr = rows[inr]; loc = ent[inr] - lo
            cv[run]['G'] += wcv[sig].sum()
            bins = {fs: tbin(scores[fs][rr]) for fs in FS}
            s_ = lab[rr]
            for fs in FS:
                b = bins[fs]
                np.add.at(cv[run]['S'][fs], b[s_], wcv[loc[s_]]); np.add.at(cv[run]['Bmc'][fs], b[~s_], wcv[loc[~s_]])
                np.add.at(cv[run]['W2'][fs], b, wcv[loc] ** 2)
            need = np.zeros(ne, bool); need[sig] = True; need[loc] = True
            idx = np.nonzero(need)[0]
            pos = np.full(ne, -1); pos[idx] = np.arange(len(idx))
            for wb, n in nuniv.items():
                term, multi, conv = UNIV[wb]
                U = np.ones((len(idx), n)); full = np.zeros(len(idx), bool)
                if wb in a.fields:
                    x = a[wb][idx]
                    full = ak.to_numpy(ak.num(x)) == n
                    if full.any(): U[full] = ak.to_numpy(x[full]).astype(float)
                base = (comp['ppfx_cv_weight'] * comp['normalisation_weight'] if conv == 'genie'
                        else comp['tuned_cv_weight'] * comp['ppfx_cv_weight'] * comp['normalisation_weight'])[idx]
                WU = safe(U * base[:, None]) * scale
                WU[~full] = wcv[idx][~full, None]          # no weights for this event: as the CV
                A = acc[(run, wb)]
                A['G'] += WU[pos[sig]].sum(axis=0)
                for fs in FS:
                    b = bins[fs]
                    np.add.at(A['S'][fs], b[s_], WU[pos[loc[s_]]]); np.add.at(A['B'][fs], b[~s_], WU[pos[loc[~s_]]])
    def above(h):                # yields above each threshold from the per-bin sums
        return np.cumsum(h[::-1], axis=0)[::-1][1:]
    proj = {r: sum(pot[x] for x in (FHC_RUNS if r in FHC_RUNS else RHC_RUNS)) / pot[r] for r in TEST_RUNS}
    out = {}
    with np.errstate(divide='ignore', invalid='ignore'):
        for fs in FS:
            out[fs] = {}
            for cfg, runs in (('FHC', (5,)), ('RHC', (13,)), ('COMB', (5, 13))):
                P = lambda key: sum(above(cv[r][key][fs]) * proj[r] for r in runs)
                S, Bmc, Bext = P('S'), P('Bmc'), P('Bext')
                G = sum(cv[r]['G'] * proj[r] for r in runs)
                Bk = Bmc + Bext; N = S + Bk
                terms = {}
                for wb, n in nuniv.items():
                    term, multi, conv = UNIV[wb]
                    Su = sum(above(acc[(r, wb)]['S'][fs]) * proj[r] for r in runs)
                    Bu = sum(above(acc[(r, wb)]['B'][fs]) * proj[r] for r in runs) + Bext[:, None]
                    Gu = sum(acc[(r, wb)]['G'] * proj[r] for r in runs)
                    ru = (N[:, None] - Bu) * Gu[None, :] / (Su * G)
                    d = ru - 1.
                    terms[term] = np.sqrt(np.nanmean(d ** 2, axis=1)) if multi else np.abs(d[:, 0])
                xsec = np.sqrt(sum(v ** 2 for k, v in terms.items() if k.startswith('xsec')))
                raw = sum(above(cv[r]['S'][fs]) for r in runs)
                out[fs][cfg] = dict(signal=S.tolist(), background=Bk.tolist(), beam_off=Bext.tolist(), efficiency=(S / G).tolist(),
                                    stat=(np.sqrt(N) / S).tolist(), flux=(FLUX * (1 + Bk / S)).tolist(),
                                    pot_targets=(np.hypot(0.02, 0.01) * (1 + Bk / S)).tolist(),
                                    mcstat=(np.sqrt(sum(above(cv[r]['W2'][fs]) for r in runs)) / raw).tolist(),
                                    xsec=xsec.tolist(), reint=terms.get('reint', np.zeros(T)).tolist(),
                                    terms={k: v.tolist() for k, v in terms.items()})
    json.dump({'thresholds': SYST_THR.tolist(), 'nuniv': nuniv, 'curves': out}, open(NEW + f'syst_n{npi}{VTAG}.json', 'w'), indent=1)
    print('universe terms written to', NEW + f'syst_n{npi}{VTAG}.json', flush=True)


def detsyst():
    """Detector term of the one-bin total against the event-score threshold, from the Run-4 detector-variation
    samples (build_detvar). Each variation v changes, per horn mode, the selected signal by e_v = S_v / S_CV
    and the selected neutrino background by b_v = B_v / B_CV (each sample at its own POT, the same generated
    events); applied to the nominal test-period yields with the beam-off unchanged, the extracted total moves
    by (N - B_mc b_v - B_ext) / (S e_v) relative to the CV. The variations are added in quadrature, FHC and RHC
    fully correlated in the combination."""
    import xgboost as xgb
    npi = 3 if '--npi' in sys.argv and sys.argv[sys.argv.index('--npi') + 1] == '3' else 2
    V = dict(np.load(NEW + 'detvar_events.npz'))
    vs = json.load(open(NEW + 'detvar_sums.json'))
    SY = json.load(open(NEW + f'syst_n{npi}{VTAG}.json'))
    thr = np.array(SY['thresholds'])
    lab = V[f'sig{npi}'].astype(bool); dom = V[f'n{npi}_d6'].astype(bool)
    if ODD:
        dom &= np.load(NEW + 'detvar_evt.npy') % 2 == 1
    out = {}
    for fs in ('pid', 'full') + tuple('cmp_' + c for c in CMP[npi]):
        if fs.startswith('cmp_'):
            sc = V['sel_' + fs[4:]].astype(float)
            if ODD: sc = np.where(np.load(NEW + 'detvar_evt.npy') % 2 == 1, sc, -1.)
        else:
            names = feat_names(npi, fs == 'full')
            bst = xgb.Booster(); bst.load_model(NEW + f'evt_clf_n{npi}{VTAG}_{fs}.json')
            X = np.stack([V[f'n{npi}_{k}'] for k in names], axis=1)
            sc = np.where(dom, bst.predict(xgb.DMatrix(X, feature_names=names)), -1.)
        # yields above each threshold per sample (file id -> horn, variation)
        Y = {}
        for fid, (run, var, path, scale) in enumerate(vs['files']):
            m = V['fid'] == fid
            horn = 'FHC' if run == 4 else 'RHC'
            ss, ww, ll = sc[m], V['w'][m], lab[m]
            Y[(horn, var)] = dict(S=np.array([ww[ll & (ss > t)].sum() for t in thr]),
                                  B=np.array([ww[~ll & (ss > t)].sum() for t in thr]),
                                  nS=np.array([int((ll & (ss > t)).sum()) for t in thr]))
        nom = SY['curves'][fs]
        res = {}
        for cfg, horns in (('FHC', ('FHC',)), ('RHC', ('RHC',)), ('COMB', ('FHC', 'RHC'))):
            per_var = {}
            for var in DETVAR[1:]:
                num = 0.; Ssum = 0.; Sv = 0.; N = 0.
                for h in horns:
                    S = np.array(nom[h]['signal']); Bk = np.array(nom[h]['background'])
                    Bmc = Bk * 0.                     # filled below from the composition of the nominal background
                    cv_, v_ = Y[(h, 'CV')], Y[(h, var)]
                    with np.errstate(divide='ignore', invalid='ignore'):
                        e = v_['S'] / cv_['S']; bb = v_['B'] / cv_['B']
                    ext = np.array(SY['curves'][fs][h].get('beam_off', np.zeros(len(thr))))
                    Bmc = Bk - ext
                    N += S + Bk
                    num = num + (S + Bk - Bmc * bb - ext)
                    Ssum = Ssum + S * e
                    Sv = Sv + S
                with np.errstate(divide='ignore', invalid='ignore'):
                    per_var[var] = (num / Ssum) - 1.
            tot = np.sqrt(sum(np.nan_to_num(v) ** 2 for v in per_var.values()))
            res[cfg] = dict(total=tot.tolist(), per_variation={k: v.tolist() for k, v in per_var.items()},
                            cv_selected_signal_raw={h: Y[(h, 'CV')]['nS'].tolist() for h in horns})
        out[fs] = res
    otag = VTAG + ('_odd' if ODD else '')
    json.dump({'thresholds': thr.tolist(), 'curves': out}, open(NEW + f'detsyst_n{npi}{otag}.json', 'w'), indent=1)
    print('detector term written to', NEW + f'detsyst_n{npi}{otag}.json', flush=True)


def detboot():
    """Statistical part of the detector term. Bootstrap in which every physical event (run, subrun, event, per horn
    mode) gets one Poisson(1) weight shared by the CV and the eight variation samples, so the correlation of the
    samples through their common generated events is kept; for each replica the change of the extracted one-bin
    total per variation is recomputed as in detsyst. The spread over replicas is the statistical uncertainty of each
    variation's shift; sqrt(sum of the squared spreads) is what statistics alone add to the quadrature sum."""
    import xgboost as xgb
    npi = 3 if '--npi' in sys.argv and sys.argv[sys.argv.index('--npi') + 1] == '3' else 2
    R = 200
    V = dict(np.load(NEW + 'detvar_events.npz'))
    ev = np.load(NEW + 'detvar_evt.npy'); rs = np.load(NEW + 'detvar_runsub.npy')
    vs = json.load(open(NEW + 'detvar_sums.json'))
    SY = json.load(open(NEW + f'syst_n{npi}{VTAG}.json'))
    thr = np.array(SY['thresholds'])
    lab = V[f'sig{npi}'].astype(bool); dom = V[f'n{npi}_d6'].astype(bool); w = V['w']
    horn_of = {fid: ('FHC' if run == 4 else 'RHC') for fid, (run, var, path, scale) in enumerate(vs['files'])}
    var_of = {fid: var for fid, (run, var, path, scale) in enumerate(vs['files'])}
    hrow = np.array([horn_of[f] for f in V['fid']])
    rng = np.random.default_rng(20261007)
    inv = np.zeros(len(w), dtype=np.int64); PW = {}
    for h in ('FHC', 'RHC'):
        m = np.flatnonzero(hrow == h)
        uk, iv = np.unique(np.stack([rs[m, 0], rs[m, 1], ev[m]], axis=1), axis=0, return_inverse=True)
        inv[m] = iv.ravel(); PW[h] = rng.poisson(1., size=(len(uk), R)).astype(np.float64)
    cases = [('pid', t) for t in (0.5, 0.8, 0.9, 0.92, 0.95)] + [('full', t) for t in (0.9, 0.93, 0.94)]
    cases += [('cmp_' + c, 0.5) for c in CMP[npi]]
    models = {}
    out = {}
    for fs, t in cases:
        if fs.startswith('cmp_'):
            sc = V['sel_' + fs[4:]].astype(float)
        else:
            if fs not in models:
                names = feat_names(npi, fs == 'full')
                bst = xgb.Booster(); bst.load_model(NEW + f'evt_clf_n{npi}{VTAG}_{fs}.json')
                X = np.stack([V[f'n{npi}_{k}'] for k in names], axis=1)
                models[fs] = np.where(dom, bst.predict(xgb.DMatrix(X, feature_names=names)), -1.)
            sc = models[fs]
        i = int(np.argmin(abs(thr - t)))
        Y = {}
        for fid in horn_of:
            h = horn_of[fid]; m = (V['fid'] == fid) & (sc > t)
            for part, mm in (('S', m & lab), ('B', m & ~lab)):
                k = np.flatnonzero(mm)
                Y[(h, var_of[fid], part)] = (w[k].sum(), (w[k][:, None] * PW[h][inv[k]]).sum(axis=0))
        nom = SY['curves'][fs]
        res = {}
        for cfg, horns in (('FHC', ('FHC',)), ('RHC', ('RHC',)), ('COMB', ('FHC', 'RHC'))):
            per = {}
            for var in DETVAR[1:]:
                num_c = den_c = 0.; num_r = den_r = 0.
                for h in horns:
                    S = nom[h]['signal'][i]; Bk = nom[h]['background'][i]; ext = nom[h]['beam_off'][i]; Bmc = Bk - ext
                    for central in (True, False):
                        j = 0 if central else 1
                        e = Y[(h, var, 'S')][j] / Y[(h, 'CV', 'S')][j]; b = Y[(h, var, 'B')][j] / Y[(h, 'CV', 'B')][j]
                        if central: num_c += S + Bk - Bmc * b - ext; den_c += S * e
                        else: num_r = num_r + (S + Bk - Bmc * b - ext); den_r = den_r + S * e
                shift = num_c / den_c - 1.; reps = num_r / den_r - 1.
                per[var] = dict(shift=float(shift), stat=float(np.std(reps)))
            tot = float(np.sqrt(sum(v['shift'] ** 2 for v in per.values())))
            noise = float(np.sqrt(sum(v['stat'] ** 2 for v in per.values())))
            res[cfg] = dict(per_variation=per, total=tot, stat_part=noise,
                            stat_subtracted=float(np.sqrt(max(tot ** 2 - noise ** 2, 0.))))
        out[f'{fs} > {t}' if not fs.startswith('cmp_') else fs[4:]] = res
    json.dump(out, open(NEW + f'detboot_n{npi}{VTAG}.json', 'w'), indent=1)
    L = [f'# Statistical part of the detector term ({npi}π, one-bin total)', '',
         f'Generated by `scripts/mp_evt.py detboot`. Paired bootstrap ({R} replicas): one Poisson weight per physical event (run, subrun, '
         'event) shared by the CV and the eight Run-4 variation samples. Shift = change of the extracted total per variation; '
         'stat. = its bootstrap spread. Statistics alone add sqrt(sum of the squared spreads) to the quadrature sum.', '',
         '| Selection (COMB) | Detector term | From statistics | Statistics subtracted | Recomb2 | WMAngleYZ | WMX | WMYZ | WMAngleXZ | SCE | LYdown | LYrayl |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for name, res in out.items():
        c = res['COMB']; pv = c['per_variation']
        L.append(f'| {name} | {100 * c["total"]:.1f}% | {100 * c["stat_part"]:.1f}% | {100 * c["stat_subtracted"]:.1f}% | '
                 + ' | '.join(f'{100 * pv[v]["shift"]:+.1f} ± {100 * pv[v]["stat"]:.1f}' for v in
                              ('Recomb2', 'WMAngleYZ', 'WMX', 'WMYZ', 'WMAngleXZ', 'SCE', 'LYdown', 'LYrayl')) + ' |')
    L += ['', '| Selection | FHC: term / from statistics | RHC: term / from statistics |', '|---|---|---|']
    for name, res in out.items():
        L.append(f'| {name} | {100 * res["FHC"]["total"]:.1f}% / {100 * res["FHC"]["stat_part"]:.1f}% | '
                 f'{100 * res["RHC"]["total"]:.1f}% / {100 * res["RHC"]["stat_part"]:.1f}% |')
    open(os.path.join(OUTDIR, f'phase2_detboot{"" if npi == 2 else "_3pi"}{VTAG}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


# ---------------------------------------------------------------- Phase 4: sideband constraint (feasibility)
REGIONS = ['SR', 'SB 3pi-like', 'SB p-like', 'SB score 0.5-T', 'SB score < 0.5']
FLUXW = 'weight_ppfx_all'


def constraint():
    """Feasibility of the sideband constraint of Phase 4, in the formalism of the framework (ConstrainedCalculator):
    every systematic is a change of the predicted signal + background in each region; for the cross-section and
    reinteraction universes the signal changes through the efficiency (smearceptance x CV true rate), for the flux
    universes through the selected signal rate (smearceptance numerator varied, denominator CV); the detector variations
    as in detsyst. The covariance of the regions (all sources, data statistics of every region, MC statistics) is
    conditioned on the sidebands, and the one-bin uncertainty is sqrt(C_SR,SR | sidebands) / S.
    Regions (exclusive, inside the domain): SR = score > T; below it the three-pion-like events (four or more tracks,
    three pions more likely than two), then the events with a proton-like assigned pion (P(p) > 0.5), then the score band
    0.5-T, then the rest; the sidebands exclude the events of the blind single-pion signal region (--with-1pi-sr keeps
    them). Test periods projected to the full exposure, as everywhere else.
        python3 scripts/mp_evt.py constraint [--tag X] [--threshold T]"""
    import xgboost as xgb
    T_ = float(sys.argv[sys.argv.index('--threshold') + 1]) if '--threshold' in sys.argv else 0.92
    D = dict(np.load(NEW + 'events.npz'))
    sums = json.load(open(NEW + 'sums.json'))
    pot = {int(r): v for r, v in sums['pot'].items()}
    lab = D['sig2'].astype(bool); dom = D['n2_d6'].astype(bool)
    sc = np.load(NEW + f'scores_n2{VTAG}_pid.npy')
    reg = np.full(len(lab), -1)
    sr = dom & (sc > T_)
    # the sidebands leave out the events of the blind single-pion signal region, so that they can be opened now
    ok1 = ~D['sel_1pi'].astype(bool) if '--with-1pi-sr' not in sys.argv else np.ones(len(lab), bool)
    three = dom & ~sr & ok1 & (D['n2_n_dump'] >= 4) & (D['n2_dL_down'] < 0)
    plike = dom & ~sr & ok1 & ~three & (np.fmax(D['n2_pi0_Pp'], D['n2_pi1_Pp']) > 0.5)
    band = dom & ~sr & ok1 & ~three & ~plike & (sc > 0.5)
    low = dom & ~sr & ok1 & ~three & ~plike & ~band
    for g, m in enumerate((sr, three, plike, band, low)): reg[m] = g
    G_ = len(REGIONS)
    univ = dict(UNIV); univ[FLUXW] = ('flux', True, 'ppfx')
    test_files = [(fid, f) for fid, f in enumerate(sums['files']) if f[0] in TEST_RUNS]
    tt = uproot.open(next(f[2] for _, f in test_files if f[1] == 'mc'))['stv_tree']
    nuniv = {wb: int(ak.max(ak.num(tt[wb].array(entry_stop=5000)))) for wb in univ if wb in tt.keys()}
    nuniv = {k: v for k, v in nuniv.items() if v > 0}
    z = lambda *sh: np.zeros(sh)
    cv = {r: dict(S=z(G_), Bmc=z(G_), Bext=z(G_), W2=z(G_), G=0.) for r in TEST_RUNS}
    acc = {(r, wb): dict(S=z(G_, n), B=z(G_, n), G=z(n)) for r in TEST_RUNS for wb, n in nuniv.items()}
    for fid, (run, kind, path, scale) in test_files:
        rows = np.flatnonzero((D['fid'] == fid) & (reg >= 0))
        if kind == 'ext':
            np.add.at(cv[run]['Bext'], reg[rows], D['w'][rows]); np.add.at(cv[run]['W2'], reg[rows], D['w'][rows] ** 2)
            continue
        ent = D['entry'][rows]
        t = uproot.open(path)['stv_tree']; keys = set(t.keys())
        br = W + [E + 'MC_Signal'] + [wb for wb in nuniv if wb in keys]
        print(f'  run {run} {kind} {os.path.basename(path)}', flush=True)
        for a, rep in t.iterate(br, step_size=100000, library='ak', report=True):
            lo, hi = rep.tree_entry_start, rep.tree_entry_stop; ne = hi - lo
            comp = {x: ak.to_numpy(a[x]).astype(float) for x in W}
            wcv = safe(comp['spline_weight'] * comp['tuned_cv_weight'] * comp['ppfx_cv_weight'] * comp['normalisation_weight']) * scale
            sig = ak.to_numpy(a[E + 'MC_Signal']).astype(bool) if kind == 'mc' else np.zeros(ne, bool)
            inr = (ent >= lo) & (ent < hi); rr = rows[inr]; loc = ent[inr] - lo
            g_ = reg[rr]; s_ = lab[rr]
            cv[run]['G'] += wcv[sig].sum()
            np.add.at(cv[run]['S'], g_[s_], wcv[loc[s_]]); np.add.at(cv[run]['Bmc'], g_[~s_], wcv[loc[~s_]])
            np.add.at(cv[run]['W2'], g_, wcv[loc] ** 2)
            need = np.zeros(ne, bool); need[sig] = True; need[loc] = True
            idx = np.flatnonzero(need); pos = np.full(ne, -1); pos[idx] = np.arange(len(idx))
            for wb, n in nuniv.items():
                term, multi, conv = univ[wb]
                U = np.ones((len(idx), n)); full = np.zeros(len(idx), bool)
                if wb in a.fields:
                    x = a[wb][idx]; full = ak.to_numpy(ak.num(x)) == n
                    if full.any(): U[full] = ak.to_numpy(x[full]).astype(float)
                if conv == 'genie': base = comp['ppfx_cv_weight'] * comp['normalisation_weight']
                elif conv == 'ppfx': base = comp['tuned_cv_weight'] * comp['normalisation_weight']
                else: base = comp['tuned_cv_weight'] * comp['ppfx_cv_weight'] * comp['normalisation_weight']
                WU = safe(U * base[idx][:, None]) * scale
                WU[~full] = wcv[idx][~full, None]
                A = acc[(run, wb)]
                A['G'] += WU[pos[sig]].sum(axis=0)
                np.add.at(A['S'], g_[s_], WU[pos[loc[s_]]]); np.add.at(A['B'], g_[~s_], WU[pos[loc[~s_]]])
    proj = {r: sum(pot[x] for x in (FHC_RUNS if r in FHC_RUNS else RHC_RUNS)) / pot[r] for r in TEST_RUNS}
    P = lambda key: sum(cv[r][key] * proj[r] for r in TEST_RUNS)
    S, Bmc, Bext = P('S'), P('Bmc'), P('Bext'); Gcv = sum(cv[r]['G'] * proj[r] for r in TEST_RUNS)
    N = S + Bmc + Bext
    # change of the predicted signal + background per universe and region (framework conventions)
    covs, deltas, UNI = {}, {}, {}
    for wb, n in nuniv.items():
        term, multi, conv = univ[wb]
        Su = sum(acc[(r, wb)]['S'] * proj[r] for r in TEST_RUNS); Bu = sum(acc[(r, wb)]['B'] * proj[r] for r in TEST_RUNS)
        Gu = sum(acc[(r, wb)]['G'] * proj[r] for r in TEST_RUNS)
        if conv == 'ppfx': dS = Su - S[:, None]
        else: dS = S[:, None] * ((Su / Gu[None, :]) / (S / Gcv)[:, None] - 1.)
        d = np.nan_to_num(dS) + (Bu - Bmc[:, None])
        d = np.vstack([d, (Bu[0] - Bmc[0])[None, :]])          # row G_: background only, signal region
        C = (d @ d.T) / n if multi else np.outer(d[:, 0], d[:, 0])
        UNI[wb] = (Su, Bu, Gu)
        covs[term] = covs.get(term, 0.) + C
    # detector variations: ratios of the Run-4 variation samples, per horn, applied to the nominal yields of that horn
    Vd = dict(np.load(NEW + 'detvar_events.npz')); vs = json.load(open(NEW + 'detvar_sums.json'))
    import xgboost as xgb
    names = feat_names(2, False); bst = xgb.Booster(); bst.load_model(NEW + f'evt_clf_n2{VTAG}_pid.json')
    scd = bst.predict(xgb.DMatrix(np.stack([Vd[f'n2_{k}'] for k in names], 1), feature_names=names))
    domd = Vd['n2_d6'].astype(bool); labd = Vd['sig2'].astype(bool)
    srd = domd & (scd > T_)
    ok1d = ~Vd['sel_1pi'].astype(bool) if '--with-1pi-sr' not in sys.argv else np.ones(len(labd), bool)
    threed = domd & ~srd & ok1d & (Vd['n2_n_dump'] >= 4) & (Vd['n2_dL_down'] < 0)
    plid = domd & ~srd & ok1d & ~threed & (np.fmax(Vd['n2_pi0_Pp'], Vd['n2_pi1_Pp']) > 0.5)
    bandd = domd & ~srd & ok1d & ~threed & ~plid & (scd > 0.5)
    lowd = domd & ~srd & ok1d & ~threed & ~plid & ~bandd
    regd = np.full(len(labd), -1)
    for g, m in enumerate((srd, threed, plid, bandd, lowd)): regd[m] = g
    Y = {}
    for fid, (run, var, path, scale) in enumerate(vs['files']):
        h = 'FHC' if run == 4 else 'RHC'; m = (Vd['fid'] == fid) & (regd >= 0)
        Sv = np.zeros(G_); Bv = np.zeros(G_)
        np.add.at(Sv, regd[m & labd], Vd['w'][m & labd]); np.add.at(Bv, regd[m & ~labd], Vd['w'][m & ~labd])
        Y[(h, var)] = (Sv, Bv)
    per_h = {h: dict(S=cv[r]['S'] * proj[r], Bmc=cv[r]['Bmc'] * proj[r]) for h, r in (('FHC', 5), ('RHC', 13))}
    Cdet = 0.
    for var in DETVAR[1:]:
        d = np.zeros(G_ + 1)
        for h in ('FHC', 'RHC'):
            with np.errstate(divide='ignore', invalid='ignore'):
                e = np.nan_to_num(Y[(h, var)][0] / Y[(h, 'CV')][0], nan=1.); b = np.nan_to_num(Y[(h, var)][1] / Y[(h, 'CV')][1], nan=1.)
            d[:G_] += per_h[h]['S'] * (e - 1.) + per_h[h]['Bmc'] * (b - 1.)
            d[G_] += per_h[h]['Bmc'][0] * (b[0] - 1.)
        Cdet = Cdet + np.outer(d, d)
    covs['detector'] = Cdet
    norm = np.hypot(0.02, 0.01); NB = np.append(N, Bmc[0]); covs['pot_targets'] = np.outer(norm * NB, norm * NB)
    covs['data_stat'] = np.diag(np.append(N, 0.))
    rel_mc = np.sqrt(sum(cv[r]['W2'] for r in TEST_RUNS)) / np.maximum(sum(cv[r]['S'] + cv[r]['Bmc'] + cv[r]['Bext'] for r in TEST_RUNS), 1e-12)
    mcs = np.diag(np.append((rel_mc * N) ** 2, (rel_mc[0] * Bmc[0]) ** 2)); mcs[0, G_] = mcs[G_, 0] = (rel_mc[0] * Bmc[0]) ** 2
    covs['mc_stat'] = mcs
    Ctot = sum(covs.values())
    groups = {'flux': ['flux'], 'cross section': [k for k in covs if k.startswith('xsec')], 'reinteraction': ['reint'],
              'detector': ['detector'], 'POT, targets': ['pot_targets'], 'data statistics': ['data_stat'], 'MC statistics': ['mc_stat']}
    def residual(Q):
        """Per-group residual variance of SR after the constraint with the sideband set Q (K from the total covariance)."""
        if not Q: K = np.zeros((1, 0))
        else: K = Ctot[np.ix_([0], Q)] @ np.linalg.inv(Ctot[np.ix_(Q, Q)])
        out = {}
        for gname, ks in groups.items():
            C = sum(covs[k] for k in ks)
            v = C[0, 0]
            if Q: v = v - 2 * (K @ C[np.ix_(Q, [0])])[0, 0] + (K @ C[np.ix_(Q, Q)] @ K.T)[0, 0]
            out[gname] = float(np.sqrt(max(v, 0.)) / S[0])
        tot = Ctot[0, 0] - (0. if not Q else (Ctot[np.ix_([0], Q)] @ np.linalg.inv(Ctot[np.ix_(Q, Q)]) @ Ctot[np.ix_(Q, [0])])[0, 0])
        out['total'] = float(np.sqrt(tot) / S[0])
        return out
    sets = {'none': [], 'three-pion-like': [1], 'proton-like': [2], 'score 0.5-T': [3], 'score < 0.5': [4],
            'all sidebands': [1, 2, 3, 4]}
    res = {'threshold': T_, 'regions': REGIONS, 'signal': S.tolist(), 'background_mc': Bmc.tolist(), 'beam_off': Bext.tolist(),
           'purity': (S / N).tolist(), 'constraint': {k: residual(v) for k, v in sets.items()}}
    # correlation of each sideband with the SR background-driven prediction, per group
    res['sr_sb_correlation_total'] = (Ctot[0] / np.sqrt(np.diag(Ctot) * Ctot[0, 0])).tolist()
    # composition of each region by true class (test periods projected)
    comp_ = {}
    for g, name in enumerate(REGIONS):
        m = reg == g
        c = {}
        for i, nm in enumerate(TOPO):
            c[nm] = float(sum(D['w'][m & np.isin(D['run'], [r]) & (D['kind'] == 0) & (D['topology'] == i)].sum() * proj[r] for r in TEST_RUNS))
        c['dirt'] = float(sum(D['w'][m & np.isin(D['run'], [r]) & (D['kind'] == 2)].sum() * proj[r] for r in TEST_RUNS))
        c['beam-off'] = float(sum(D['w'][m & np.isin(D['run'], [r]) & (D['kind'] == 1)].sum() * proj[r] for r in TEST_RUNS))
        comp_[name] = c
    res['composition'] = comp_
    # alternative-model closure: pseudo-data from a reweighted model in every region, extracted with the CV model,
    # unconstrained and with each sideband set (background estimate B_CV + K_B (N_Q - N_Q,CV), K_B from the total
    # covariance), against that model's true signal rate
    def extract(Nm, Q):
        if not Q: Bh = Bmc[0]
        else: Bh = Bmc[0] + (Ctot[np.ix_([G_], Q)] @ np.linalg.inv(Ctot[np.ix_(Q, Q)]) @ (Nm[Q] - N[Q]))[0]
        return (Nm[0] - Bh - Bext[0]) / S[0]
    alt = {}
    models = [('GENIE universe 545', 'weight_All_UBGenie', 545), ('Delta->N pi angular', 'weight_Theta_Delta2Npi_UBGenie', 0)]
    for mname, wb, u in models:
        if wb not in UNI: continue
        Su, Bu, Gu = UNI[wb]
        Nm = Su[:, u] + Bu[:, u] + Bext
        truth = Gu[u] / Gcv
        alt[mname] = {k: dict(ratio=float(extract(Nm, v) / truth), unc=res['constraint'][k]['total']) for k, v in sets.items()}
    # coverage over the GENIE multisim universes: fraction with |bias| below the total uncertainty
    if 'weight_All_UBGenie' in UNI:
        Su, Bu, Gu = UNI['weight_All_UBGenie']
        cov = {}
        for k, v in sets.items():
            b_ = np.array([extract(Su[:, u] + Bu[:, u] + Bext, v) / (Gu[u] / Gcv) - 1. for u in range(Su.shape[1])])
            cov[k] = dict(rms_bias=float(np.sqrt(np.mean(b_ ** 2))), within_1sigma=float(np.mean(np.abs(b_) < res['constraint'][k]['total'])))
        res['genie_multisim_coverage'] = cov
    res['alternative_models'] = alt
    json.dump(res, open(os.path.join(OUTDIR, f'phase4_constraint{VTAG}.json'), 'w'), indent=1)
    L = [f'# Phase 4 feasibility: sideband constraint of the two-pion one-bin total ({("deployed" if not VTAG else VTAG[1:])} particle classifier, PID set, score > {T_})', '',
         'Generated by `scripts/mp_evt.py constraint`. Formalism of the framework\'s ConstrainedCalculator: covariance of the predicted '
         'signal + background in every region from all sources (flux from the PPFX universes, cross section and reinteraction from the '
         'universe weights, detector from the eight Run-4 variations, POT and targets, data and MC statistics), conditioned on the '
         'sidebands. Test periods projected to the full exposure (COMB).', '',
         '| Region | Events | Signal | Purity | Largest classes |', '|---|---|---|---|---|']
    for g, name in enumerate(REGIONS):
        top = sorted(comp_[name].items(), key=lambda x: -x[1])[:3]
        L.append(f'| {name} | {N[g]:.0f} | {S[g]:.0f} | {100 * S[g] / N[g]:.0f}% | ' + ', '.join(f'{k} {100 * v / N[g]:.0f}%' for k, v in top) + ' |')
    L += ['', '| Constraint | Total | Flux | Cross section | Reinteraction | Detector | POT, targets | Data statistics | MC statistics |',
          '|---|---|---|---|---|---|---|---|---|']
    for k, r in res['constraint'].items():
        L.append(f'| {k} | {100 * r["total"]:.1f}% | ' + ' | '.join(f'{100 * r[gn]:.1f}%' for gn in
                 ('flux', 'cross section', 'reinteraction', 'detector', 'POT, targets', 'data statistics', 'MC statistics')) + ' |')
    L += ['', '### Alternative-model closure (extracted over true signal rate; uncertainty = total of the row above)', '',
          '| Constraint | ' + ' | '.join(alt) + ' | GENIE multisim: rms bias | within 1σ |', '|---|' + '---|' * (len(alt) + 2)]
    for k in sets:
        cg = res.get('genie_multisim_coverage', {}).get(k, {})
        L.append(f'| {k} | ' + ' | '.join(f'{alt[m][k]["ratio"]:.3f} ({(alt[m][k]["ratio"] - 1) / alt[m][k]["unc"]:+.2f}σ)' for m in alt)
                 + f' | {100 * cg.get("rms_bias", float("nan")):.1f}% | {100 * cg.get("within_1sigma", float("nan")):.0f}% |')
    open(os.path.join(OUTDIR, f'phase4_constraint{VTAG}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


TERMS = ['stat', 'flux', 'pot_targets', 'xsec', 'reint', 'detector', 'mcstat']


def wp():
    """Expected total uncertainty of the one-bin total against the threshold, all terms combined, and the
    working point at its minimum (COMB); the current selections with the same terms for comparison."""
    npi = 3 if '--npi' in sys.argv and sys.argv[sys.argv.index('--npi') + 1] == '3' else 2
    tag = '' if npi == 2 else '_3pi'
    SY = json.load(open(NEW + f'syst_n{npi}{VTAG}.json'))
    dpath = NEW + f'detsyst_n{npi}{VTAG}{"_odd" if ODD else ""}.json'
    DS = json.load(open(dpath)) if os.path.exists(dpath) else None
    thr = np.array(SY['thresholds'])
    out = {'thresholds': thr.tolist(), 'detector_included': DS is not None, 'sets': {}}
    for fs, cur in SY['curves'].items():
        o = {}
        for cfg in ('FHC', 'RHC', 'COMB'):
            c = {k: np.array(cur[cfg][k], dtype=float) for k in ('signal', 'background', 'efficiency', 'stat', 'flux',
                                                                  'pot_targets', 'xsec', 'reint', 'mcstat')}
            c['detector'] = np.array(DS['curves'][fs][cfg]['total']) if DS else np.zeros(len(thr))
            with np.errstate(invalid='ignore'):
                c['total'] = np.sqrt(sum(np.nan_to_num(c[k], nan=np.inf) ** 2 for k in TERMS))
                c['purity'] = c['signal'] / (c['signal'] + c['background'])
            o[cfg] = {k: v.tolist() for k, v in c.items()}
            if DS: o[cfg]['detvar_cv_selected_signal_raw'] = DS['curves'][fs][cfg]['cv_selected_signal_raw']
        i = int(np.nanargmin(o['COMB']['total'])) if not fs.startswith('cmp_') else 0
        o['best_index'] = i
        out['sets'][fs] = o
    json.dump(out, open(os.path.join(OUTDIR, f'phase2_wp{tag}{VTAG}{"_odd" if ODD else ""}.json'), 'w'), indent=1)
    L = [f'# Phase 2, item 4: working point of the {npi}π event classifier on the expected total uncertainty', '',
         f'Generated by `scripts/mp_evt.py wp{" --npi 3" if npi == 3 else ""}` from the universe step (`syst`) and the detector step (`detsyst`). '
         'One-bin total, test periods (FHC Run 5, RHC Run 3) projected to the full exposure. Terms: data statistics; flux 17% × (1 + B/S); '
         'POT and targets (2% and 1%) × (1 + B/S); cross section (GENIE multisim, RPA and SCC multisims, eight GENIE unisims, as '
         '`configs/ccpi_systcalc_numi.conf`); reinteraction; detector (eight Run-4 variations, FHC and RHC correlated)'
         + ('' if DS else ' NOT YET INCLUDED') + '; MC statistics of the test samples. The current selections carry the same terms.', '',
         '| Selection (COMB) | Score > | Efficiency | Purity | Stat. | Flux | POT, targets | Cross section | Reint. | Detector | MC stat. | Total |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    def row(name, o, i, cfg='COMB'):
        c = o[cfg]
        th = '–' if name.startswith('CC1mu') else f'{thr[i]:.4g}'
        return (f'| {name} | {th} | {100 * c["efficiency"][i]:.2f}% | {100 * c["purity"][i]:.1f}% | '
                + ' | '.join(f'{100 * c[k][i]:.1f}%' for k in TERMS) + f' | {100 * c["total"][i]:.1f}% |')
    for fs, o in out['sets'].items():
        if fs.startswith('cmp_'):
            L.append(row(fs[4:], o, 0))
    for fs in ('pid', 'full'):
        o = out['sets'][fs]
        L.append(row(f'classifier, {fs} set, minimum', o, o['best_index']))
    for fs in ('pid',):
        o = out['sets'][fs]
        L += ['', f'Scan, {fs} set (COMB):', '', '| Score > | Efficiency | Purity | Stat. | Flux | POT, targets | Cross section | Reint. | Detector | MC stat. | Total | Det. CV signal (raw, FHC/RHC) |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for i, t in enumerate(thr):
            if not (abs(t * 20 - round(t * 20)) < 1e-6 or t >= 0.95): continue
            c = o['COMB']
            raw = c.get('detvar_cv_selected_signal_raw', {})
            rw = '/'.join(str(raw[h][i]) for h in ('FHC', 'RHC') if h in raw) if raw else '–'
            L.append(f'| {t:.4g} | {100 * c["efficiency"][i]:.2f}% | {100 * c["purity"][i]:.1f}% | '
                     + ' | '.join(f'{100 * c[k][i]:.1f}%' for k in TERMS) + f' | {100 * c["total"][i]:.1f}% | {rw} |')
    open(os.path.join(OUTDIR, f'phase2_wp{tag}{VTAG}{"_odd" if ODD else ""}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:12]))


def write_md(res, tag):
    npi = res['npi']
    L = [f'# Phase 2: track assignment and {npi}π event classifier', '',
         f'Generated by `scripts/mp_evt.py train{" --npi 3" if npi == 3 else ""}`. Trained on run periods {list(res["train_runs"])}, '
         f'tested on {list(res["test_runs"])} (held out of the particle classifier as well). Test-period yields projected to the '
         'full exposure of each horn mode. Total = data statistics in quadrature with 17% × (1 + B/S). The threshold is chosen '
         f'on a random {100 * res["validation_fraction"]:.0f}% of the training periods, projected the same way; the tables quote the test periods.', '',
         f'Rows in the domain {res["rows"]["domain"]}; training {res["rows"]["train"]} ({res["raw_signal"]["train"]} signal), '
         f'test {res["rows"]["test"]} ({res["raw_signal"]["test"]} signal).', '',
         '| Feature set | Trees | AUC (validation) | AUC (test runs) |', '|---|---|---|---|']
    for fs in ('pid', 'full'):
        r = res[fs]
        L.append(f'| {fs} ({len(r["features"])} features) | {r["best_iteration"] + 1} | {r["auc_val"]:.4f} | {r["auc_test"]:.4f} |')
    L += ['', '## Working point against the current selections (COMB, full exposure)', '',
          '| Selection | Signal | Background | Efficiency | Purity | B/S | Stat. | Flux term | Total |', '|---|---|---|---|---|---|---|---|---|']
    def row(name, c):
        return (f'| {name} | {c["signal"]:.1f} | {c["background"]:.1f} | {100 * c["efficiency"]:.2f}% | {100 * c["purity"]:.1f}% | '
                f'{c["b_over_s"]:.2f} | {100 * c["stat"]:.1f}% | {100 * c["flux"]:.1f}% | {100 * c["total"]:.1f}% |')
    for name, c in res['compare'].items():
        L.append(row(name, c['COMB']))
    L.append(row(f'assignment only (most likely pion count {npi})', res['assignment_only']['COMB']))
    for fs in ('pid', 'full'):
        b = res['scan'][fs]['chosen']
        L.append(row(f'event classifier, {fs} set, score > {b["threshold"]:.3f}', b['COMB']))
        b = res['scan'][fs]['best_on_test']
        L.append(row(f'(threshold best on the test periods: score > {b["threshold"]:.3f})', b['COMB']))
    wp = res['working_point']
    L += ['', f'Diagnostics below at the working point: {wp["feature_set"]} set, score > {wp["threshold"]:.4g}' + (' (given)' if wp.get('given') else ' (chosen on the validation sample, proxy total)') + '.', '',
          '| Configuration | Signal | Background | Efficiency | Purity | B/S | Total |', '|---|---|---|---|---|---|---|']
    for cfg in ('FHC', 'RHC', 'COMB'):
        c = wp[cfg]
        L.append(f'| {cfg} | {c["signal"]:.1f} | {c["background"]:.1f} | {100 * c["efficiency"]:.2f}% | {100 * c["purity"]:.1f}% | {c["b_over_s"]:.2f} | {100 * c["total"]:.1f}% |')
    L += ['', '### Background at the working point (COMB, full exposure)', '', '| Class | Events |', '|---|---|']
    for k, v in sorted(wp['background_by_topology_COMB'].items(), key=lambda x: -x[1]):
        L.append(f'| {k} | {v:.1f} |')
    L += ['', '### Identity of the assigned tracks (test periods)', '', '| | ' + ' | '.join(wp['assigned_identity_test']['signal']) + ' |',
          '|---|' + '---|' * len(wp['assigned_identity_test']['signal'])]
    for part, d in wp['assigned_identity_test'].items():
        L.append(f'| {part} | ' + ' | '.join(f'{100 * v:.1f}%' for v in d.values()) + ' |')
    L += ['', '### Overlaps (test periods)', '', '| With | Fraction of selected signal | Fraction of all selected |', '|---|---|---|']
    for k, v in wp['overlap_test'].items():
        L.append(f'| {k} | {100 * v["signal"]:.1f}% | {100 * v["all"]:.1f}% |')
    if 'efficiency_vs_true_pion_momentum' in wp:
        L += ['', '### Efficiency against the true pion momenta (test periods)', '']
        for key, d in wp['efficiency_vs_true_pion_momentum'].items():
            e = d['edges']
            L += [f'{key.capitalize()} pion:', '', '| p (GeV/c) | ' + ' | '.join(f'{lo:g}–{hi:g}' if hi < 5 else f'> {lo:g}' for lo, hi in zip(e[:-1], e[1:])) + ' |',
                  '|---|' + '---|' * (len(e) - 1),
                  '| working point | ' + ' | '.join(f'{100 * x:.1f}%' for x in d['working_point']) + ' |',
                  '| current | ' + ' | '.join(f'{100 * x:.1f}%' for x in d['current']) + ' |', '']
    if 'three_pion_like_sample' in res:
        t3 = res['three_pion_like_sample']
        L += ['', '### Three-pion-like events (decision D8; COMB, full exposure)', '',
              f'Events of the domain with at least four tracks in which three pions are more likely than two: {t3["events_COMB"]:.1f}; '
              f'{100 * t3["overlap_with_working_point"]:.1f}% of them pass the working point.', '', '| Class | Events |', '|---|---|']
        for k, v in sorted(t3['by_class_COMB'].items(), key=lambda x: -x[1]):
            L.append(f'| {k} | {v:.1f} |')
    L += ['', '### Threshold scan (COMB, full exposure)', '', '| Set | Score > | Efficiency | Purity | B/S | Stat. | Flux term | Total |', '|---|---|---|---|---|---|---|---|']
    for fs in ('pid', 'full'):
        for r in res['scan'][fs]['test'][::5]:
            c = r['COMB']
            L.append(f'| {fs} | {r["threshold"]:.2f} | {100 * c["efficiency"]:.2f}% | {100 * c["purity"]:.1f}% | {c["b_over_s"]:.2f} | '
                     f'{100 * c["stat"]:.1f}% | {100 * c["flux"]:.1f}% | {100 * c["total"]:.1f}% |')
    L += ['', '### Largest feature importances (gain)', '']
    for fs in ('pid', 'full'):
        L.append(f'- {fs}: ' + ', '.join(f'{k} ({v:.1f})' for k, v in list(res[fs]['importance_gain'].items())[:10]))
    open(os.path.join(OUTDIR, f'phase2_event{tag}{VTAG}.md'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    {'build': build, 'build_detvar': build_detvar, 'detvar_evt': detvar_evt, 'train': train, 'syst': syst, 'detsyst': detsyst, 'wp': wp, 'detboot': detboot, 'constraint': constraint}[sys.argv[1]]()
