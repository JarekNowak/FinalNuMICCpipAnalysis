#!/usr/bin/env python3
"""mp_highn_feasibility.py -- can the multi-pion analysis go to four and five charged pions?

Reads the Phase 0 trees (/data/uboone/processed/mp_phase0: 14 overlays, 7 beam-off, 7 dirt, release file
list, normalised as the framework by scripts/mp_phase0_baseline.py) and evaluates, for N = 3, 4, 5:

  truth  the N-pion signal of CC1mu1piXp::define_signal with required_charged_pions() = N: true vertex in
         the FV, CC, |nu_pdg| = 14, one muon above threshold with p_mu > 0.15 GeV/c, exactly N charged pions
         with the softest above 0.10 GeV/c (and, as a variant, 0.175), no pi0, kaon or heavier meson,
         theta(mu, leading pi) < 2.6 rad;
  reco   the selection a CC1mu<N>pi subclass of CC1mu3pi would apply, emulated from the CC1mu3pi branches
         (its pion identification is the one such a subclass inherits): software trigger, vertex, topology,
         muon (cut-flow bits 0-3), the reco opening angle (bit 8), N counted candidates with at least one
         contained, and the multiplicity cuts relaxed by one per pion (non-protons < N + 3, primary tracks
         < N + 4).

The emulation is validated by reproducing CC1mu3pi_MC_Signal and CC1mu3pi_Selected exactly for N = 3.

    python3 scripts/mp_highn_feasibility.py  ->  report/multipion/highn_feasibility.{json,md}
"""
import json, os, sys
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_phase0_baseline as B

S = 'CC1mu3pi_'
TRUTH = ['sig_truevertex_in_fv', 'sig_ccnc', 'sig_is_numu', 'sig_one_muon_above_thresh', 'mc_n_threshold_pionpm', 'mc_n_threshold_pion0',
         'mc_n_kaons', 'mc_n_heaviermeson', 'candidate_muon_mom_true', 'mc_pionpm_min_mom', 'true_mu_pi_opening_angle', 'MC_Signal']
RECO = ['Selected', 'cutflow_bits', 'pion_number_reco', 'n_contained_pion', 'sb_nnonproton', 'sb_nprimtrk']
WEIGHTS = ['spline_weight', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
NS = (3, 4, 5)
PRE = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 8)     # trigger, vertex, topology, muon, opening angle


def signal(a, n, thr=0.10):
    g = lambda k: a[S + k]
    return (g('sig_truevertex_in_fv').astype(bool) & g('sig_ccnc').astype(bool) & g('sig_is_numu').astype(bool)
            & g('sig_one_muon_above_thresh').astype(bool) & (g('mc_n_threshold_pionpm') == n) & (g('mc_n_threshold_pion0') == 0)
            & (g('mc_n_kaons') == 0) & (g('mc_n_heaviermeson') == 0) & (g('candidate_muon_mom_true') > 0.15)
            & (g('mc_pionpm_min_mom') > thr) & (g('true_mu_pi_opening_angle') < 2.6))


def selected(a, n):
    g = lambda k: a[S + k]
    return (((g('cutflow_bits') & PRE) == PRE) & (g('n_contained_pion') >= 1) & (g('pion_number_reco') == n)
            & (g('sb_nnonproton') < n + 3) & (g('sb_nprimtrk') < n + 4))


def main():
    B.NEW = '/data/uboone/processed/mp_phase0/'
    on, ext, mc, dirt = B.parse_conf()
    sc, _ = B.scales(on, ext, mc, dirt)
    keys = [('gen', n) for n in NS] + [('gen175', n) for n in NS] + [(k, n) for k in ('sig', 'bkg', 'ext', 'dirt') for n in NS]
    acc = {c: {k: 0. for k in keys} for c in B.CONFIGS}
    bkgcls = {c: {n: np.zeros(4) for n in NS} for c in B.CONFIGS}        # true N-1, N+1 or more, other numu CC, NC/other
    check = {'signal_mismatch': 0, 'selected_mismatch': 0, 'events': 0}
    # COMB migration in the pion multiplicity: true N = 1..5+ (all pions > 0.10 GeV/c) against reco N = 0..6+,
    # for events passing the common preselection with the multiplicity cuts relaxed by the reco N
    MIG = np.zeros((6, 7)); MIGB = np.zeros(7); GEN = np.zeros(6)
    for (run, kind), files in sorted(sc.items()):
        cfgs = [c for c, runs in B.CONFIGS.items() if run in runs]
        for path, s in files:
            t = uproot.open(path)['stv_tree']
            br = [S + k for k in RECO] + ([S + k for k in TRUTH] + WEIGHTS if kind != 'ext' else [])
            a = t.arrays(br, library='np')
            w = s * (np.ones(t.num_entries) if kind == 'ext' else
                     np.where(np.isfinite(x := np.prod([a[k].astype(float) for k in WEIGHTS], axis=0)) & (x >= 0) & (x <= 30), x, 1.))
            sel3 = selected(a, 3)
            check['selected_mismatch'] += int((sel3 != a[S + 'Selected'].astype(bool)).sum()); check['events'] += t.num_entries
            if kind != 'ext':
                check['signal_mismatch'] += int((signal(a, 3) != a[S + 'MC_Signal'].astype(bool)).sum()) if kind == 'mc' else 0
            rn = np.minimum(a[S + 'pion_number_reco'], 6)
            pre = (((a[S + 'cutflow_bits'] & PRE) == PRE) & (a[S + 'n_contained_pion'] >= 1) & (a[S + 'pion_number_reco'] >= 1)
                   & (a[S + 'sb_nnonproton'] < a[S + 'pion_number_reco'] + 3) & (a[S + 'sb_nprimtrk'] < a[S + 'pion_number_reco'] + 4))
            if kind == 'mc':
                tn = np.zeros(t.num_entries, int)
                for n in range(1, 6):
                    sig = signal(a, n) if n < 5 else signal(a, 5) | signal(a, 6) | signal(a, 7)
                    tn[sig] = n; GEN[n] += w[sig].sum()
                for n in range(0, 6):
                    m = pre & (tn == n); np.add.at(MIG[n] if n else MIGB, rn[m], w[m])
            else:
                m = pre; np.add.at(MIGB, rn[m], w[m])
            for n in NS:
                sel = selected(a, n)
                if kind == 'mc':
                    sig = signal(a, n); npi = a[S + 'mc_n_threshold_pionpm']
                    cc = a[S + 'sig_ccnc'].astype(bool) & a[S + 'sig_is_numu'].astype(bool)
                    for c in cfgs:
                        acc[c][('gen', n)] += w[sig].sum(); acc[c][('gen175', n)] += w[signal(a, n, 0.175)].sum()
                        acc[c][('sig', n)] += w[sel & sig].sum(); acc[c][('bkg', n)] += w[sel & ~sig].sum()
                        m = sel & ~sig
                        bkgcls[c][n] += [w[m & cc & (npi == n - 1)].sum(), w[m & cc & (npi >= n + 1)].sum(),
                                         w[m & cc & (npi != n - 1) & (npi < n + 1)].sum(), w[m & ~cc].sum()]
                else:
                    for c in cfgs: acc[c][(kind, n)] += w[sel].sum()
            print(f'  run {run} {kind}: {os.path.basename(path)}', flush=True)
    out = {'check': check, 'configs': {}}
    for c in B.CONFIGS:
        for n in NS:
            a = acc[c]; tot = a[('sig', n)] + a[('bkg', n)] + a[('ext', n)] + a[('dirt', n)]
            out['configs'][f'{c} {n}pi'] = dict(
                generated=a[('gen', n)], generated_thr0175=a[('gen175', n)], selected_signal=a[('sig', n)], selected_total=tot,
                nu_bkg=a[('bkg', n)], ext=a[('ext', n)], dirt=a[('dirt', n)],
                efficiency=a[('sig', n)] / a[('gen', n)] if a[('gen', n)] else float('nan'),
                purity=a[('sig', n)] / tot if tot else float('nan'),
                stat_one_bin=np.sqrt(tot) / a[('sig', n)] if a[('sig', n)] else float('nan'),
                bkg_classes=dict(zip(['true N-1 pi', 'true N+1 or more pi', 'other numuCC', 'NC/nue/other'], [float(x) for x in bkgcls[c][n]])))
    out['migration_comb'] = dict(true_rows='N=1..5+', reco_cols='N=0..6+', generated=GEN[1:].tolist(), matrix=MIG[1:].tolist(),
                                 background_and_beamoff=MIGB.tolist())
    diag = {}
    for edges in ((1, 2, 3), (1, 2, 3, 4)):
        R = np.zeros((len(edges), len(edges)))
        for n in range(1, 6):
            for r in range(1, 7):
                R[edges.index(min(n, edges[-1])), edges.index(min(r, edges[-1]))] += MIG[n, r]
        diag[','.join(map(str, edges[:-1])) + f',>={edges[-1]}'] = dict(reco_given_true=(np.diag(R) / R.sum(axis=1)).tolist(),
                                                                       true_given_reco=(np.diag(R) / R.sum(axis=0)).tolist())
    out['migration_diagonals'] = diag
    R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'report', 'multipion'))
    json.dump(out, open(os.path.join(R, 'highn_feasibility.json'), 'w'), indent=1, default=float)
    L = ['# Four and five charged pions: feasibility at the full exposure', '',
         f"Phase 0 trees (current pion identification of the three-pion selection). Emulation check on {check['events']} events: "
         f"signal flag mismatches {check['signal_mismatch']}, selection mismatches {check['selected_mismatch']} (N = 3 against the stored flags).", '',
         '| Configuration | N | Generated signal | (all pions > 0.175 GeV/c) | Selected signal | Selected total | Efficiency | Purity | Stat. unc. (one bin) | Background: N-1 pi / N+1 or more / other CC / NC |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for k, v in out['configs'].items():
        c, n = k.split(); b = v['bkg_classes']
        L.append(f"| {c} | {n} | {v['generated']:.1f} | {v['generated_thr0175']:.1f} | {v['selected_signal']:.1f} | {v['selected_total']:.1f} | "
                 f"{100 * v['efficiency']:.1f}% | {100 * v['purity']:.1f}% | {100 * v['stat_one_bin']:.0f}% | "
                 f"{b['true N-1 pi']:.1f} / {b['true N+1 or more pi']:.1f} / {b['other numuCC']:.1f} / {b['NC/nue/other']:.1f} |")
    L += ['', '## Pion-multiplicity migration, combined (selected events, central-value weights, data exposure)', '',
          '| True N | Generated | reco 1 | reco 2 | reco 3 | reco 4 | reco 5 | reco 6+ |', '|---|---|---|---|---|---|---|---|']
    for i, n in enumerate(['1', '2', '3', '4', '5+']):
        L.append(f'| {n} | {GEN[i + 1]:.1f} | ' + ' | '.join(f'{x:.1f}' for x in MIG[i + 1, 1:]) + ' |')
    L.append('| no signal (incl. beam-off, dirt) | | ' + ' | '.join(f'{x:.1f}' for x in MIGB[1:]) + ' |')
    L += [''] + [f"Binning {k}: P(reco bin | true bin) {', '.join(f'{x:.2f}' for x in v['reco_given_true'])}; "
                 f"P(true bin | reco bin) {', '.join(f'{x:.2f}' for x in v['true_given_reco'])}" for k, v in diag.items()]
    open(os.path.join(R, 'highn_feasibility.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
