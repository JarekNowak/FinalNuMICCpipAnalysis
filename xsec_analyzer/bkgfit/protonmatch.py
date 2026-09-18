"""protonmatch.py -- what the proton candidate actually is, and whether it matters.

The proton-tagged observables are defined at truth level with the LEADING true proton, while the
reconstruction picks a candidate by PID, momentum and containment. Where the two disagree the
response absorbs the difference, so the question is not "how often is the pick right" on its own but
"how much of the smearing comes from the pick rather than from momentum resolution".

Both are measured here, on selected signal events of processed/w, using the friend trees of
macros/proton_angle_friend.C (prmatch: -1 no truth match, 0 not a proton, 1 a proton but not the
leading one, 2 the leading one):

  the breakdown of the candidate's truth by category;
  the reco-minus-true residual of theta_p and theta_pi-p IN EACH CATEGORY, so the smearing of a
  correct pick and of a wrong pick can be compared directly. If the residual of category 2 is as
  wide as the rest, the pick is not what limits these observables.

  python3 -m bkgfit.protonmatch
"""
import glob, os
import numpy as np, uproot

S = 'CC1mu1pi1p'
W = '/data/uboone/processed/w/'
FR = W + 'friends_pa/'
LAB = {2: 'the leading proton', 1: 'a proton, not leading', 0: 'not a proton', -1: 'no truth match'}


def load():
    """prmatch and the two angles, reco and true, for selected signal events."""
    out = {k: [] for k in ('prm', 'thp_r', 'thp_t', 'tpp_r', 'tpp_t')}
    for fr in sorted(glob.glob(FR + '*.pa.root')):
        main = W + os.path.basename(fr).replace('.pa.root', '.root')
        if not os.path.exists(main):
            print('   no main tree for', os.path.basename(fr)); continue
        a = uproot.open(fr)['pa'].arrays(['prmatch', 'thpipr_reco', 'thpipr_true',
                                          'costhp_reco', 'costhp_true'], library='np')
        b = uproot.open(main)['stv_tree'].arrays([f'{S}_Selected', f'{S}_MC_Signal'], library='np')
        if len(a['prmatch']) != len(b[f'{S}_Selected']):
            print('   LENGTH MISMATCH', os.path.basename(fr)); continue
        m = b[f'{S}_Selected'].astype(bool) & b[f'{S}_MC_Signal'].astype(bool)
        # the friend fills -9999 for events it does not treat; keep only filled angles
        ok = m & (a['costhp_reco'] > -2) & (a['costhp_true'] > -2)
        out['prm'].append(a['prmatch'][ok])
        out['thp_r'].append(np.arccos(np.clip(a['costhp_reco'][ok], -1, 1)))
        out['thp_t'].append(np.arccos(np.clip(a['costhp_true'][ok], -1, 1)))
        out['tpp_r'].append(a['thpipr_reco'][ok]); out['tpp_t'].append(a['thpipr_true'][ok])
        print('   read', os.path.basename(fr), ok.sum(), flush=True)
    return {k: np.concatenate(v) for k, v in out.items()}


def main():
    d = load()
    p = d['prm']; n = len(p)
    print(f'\nselected SIGNAL events with a proton candidate: {n}\n')
    print('%-24s %8s %7s | %-26s | %-26s' % ('candidate truth', 'events', 'share',
                                             'theta_p residual [rad]', 'theta_pi-p residual [rad]'))
    for k in (2, 1, 0, -1):
        m = p == k
        if not m.any():
            print('%-24s %8d %6.1f%%' % (LAB[k], 0, 0)); continue
        r1, r2 = d['thp_r'][m] - d['thp_t'][m], d['tpp_r'][m] - d['tpp_t'][m]
        f = lambda r: 'med %+.3f  IQR %.3f  |r|>0.3: %4.1f%%' % (
            np.median(r), np.subtract(*np.percentile(r, [75, 25])), 100 * np.mean(np.abs(r) > 0.3))
        print('%-24s %8d %6.1f%% | %s | %s' % (LAB[k], m.sum(), 100 * m.sum() / n, f(r1), f(r2)))
    good = p == 2
    print('\nany true proton: %.1f%%   leading: %.1f%%' % (100 * np.mean(p >= 1), 100 * np.mean(good)))
    for nm, r, t in (('theta_p', d['thp_r'], d['thp_t']), ('theta_pi-p', d['tpp_r'], d['tpp_t'])):
        res = r - t
        w = np.subtract(*np.percentile(res, [75, 25]))
        wg = np.subtract(*np.percentile(res[good], [75, 25]))
        # The core width is the wrong question for a binned measurement: what moves an event between
        # bins is a large residual, so attribute those instead. A wrong pick is not a mild smearing --
        # it is essentially a random angle -- and that is what feeds the migration off the diagonal.
        big = np.abs(res) > 0.3
        share = 100 * np.sum(big & ~good) / np.sum(big)
        print(f'{nm:11s} IQR all {w:.3f}, correct-pick only {wg:.3f} ({100 * (1 - wg / w):.0f}% of the core '
              f'width); |res|>0.3 rad in {100 * big.mean():.1f}% of events, of which '
              f'{share:.0f}% are wrong picks')


if __name__ == '__main__':
    main()
