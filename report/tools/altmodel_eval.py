"""altmodel_eval.py -- the alternative-model closure: how much bias survives a wrong model?

The closure reported with every released extraction throws its pseudo-data from the same simulation
that builds the response, so it tests the implementation and is silent about model robustness. These
extractions throw the pseudo-data from a reweighted interaction model and leave the response and the
background prediction nominal. Two variations are used:

  altgenie   GENIE multisim universe 545, a coherent all-parameter variation: rate x 1.04 on signal
             (so this is a SHAPE test, not a normalisation one) with true-level shape distortions of
             up to 29% in cos(theta_pi) and 24% in p_pi;
  altdelta   the Delta -> N pi decay angular unisim, a targeted change of the pion angular
             distribution in resonance events.

The framework's FakeData universe is filled from the thrown files in truth and reco, so the truth it
reports IS the alternative model's truth, and the closure it computes is the model-mismatch bias.
What matters is the size of that bias against the uncertainty the measurement quotes: a bias well
inside the systematic uncertainty is a measurement that survives the model being wrong; a bias
comparable to it is a limit on the result.

    python3 report/tools/altmodel_eval.py [FHC5|RHCFULL] [incl|1p]
    (1p = the proton-tagged family: W_had, the two-bin TKI variables and cos theta_pi, FHC only)
"""
import os, sys
import numpy as np, uproot

RB = '/data/uboone/processed/rebuild_alt'
LIVE = '/data/uboone/processed'
CFG = sys.argv[1] if len(sys.argv) > 1 else 'FHC5'     # FHC5 | RHCFULL
FAM = sys.argv[2] if len(sys.argv) > 2 else 'incl'     # incl | 1p
OBS = ['pmu', 'costhmu', 'costhpi', 'ppi2bin'] if FAM == 'incl' else ['Whad', 'dpt2bin', 'dalphat2bin', 'dphit2bin', 'pn2bin', 'costhpi']
PFX = 'ccpi' if FAM == 'incl' else 'ccpi1p'
LIVEPFX = '' if FAM == 'incl' else 'ccpi1p_'
TAGS = [('altgenie', 'GENIE multisim u545'), ('altdelta', 'Delta->N pi angular')]
LAB = {'pmu': 'p_mu', 'costhmu': 'cos th_mu', 'costhpi': 'cos th_pi', 'ppi2bin': 'p_pi (2 region)',
       'Whad': 'W_had', 'dpt2bin': 'delta p_T (2 bin)', 'dalphat2bin': 'delta alpha_T (2 bin)',
       'dphit2bin': 'delta phi_T (2 bin)', 'pn2bin': 'p_n (2 bin)'}


def load(path):
    if not os.path.exists(path):
        return None
    f = uproot.open(path)
    keys = {k.split(';')[0] for k in f.keys()}
    if 'h_unfolded_nuwro' not in keys or 'h_fakedata_truth' not in keys:
        return None
    h = f['h_unfolded_nuwro']
    return dict(v=h.values(), e=h.errors(), t=f['h_fakedata_truth'].values())


def main():
    print('%-16s %-22s %5s %8s %8s %9s %9s' % ('observable', 'variation', 'bins', 'ratio',
                                               'max dev', 'bias/sigma', 'worst bin'))
    rows = []
    for obs in OBS:
        nom = load(f'{LIVE}/closure_hists_xsec_{LIVEPFX}{CFG}_{obs}.root')
        for tag, desc in TAGS:
            d = load(f'{RB}/closure_hists_xsec_{PFX}_{CFG}_{obs}_{tag}.root')
            if d is None:
                print('%-16s %-22s %s' % (LAB[obs], desc, 'pending')); continue
            ok = d['t'] > 0
            # the ratio and its max deviation are meaningless in nearly empty bins (proton-tagged
            # W_had has three with < 10% of the peak truth and unfolded values consistent with zero),
            # so they use only bins holding >= 10% of the largest truth; the pull uses every bin
            full = ok & (d['t'] >= 0.1 * d['t'].max())
            r = d['v'][full] / d['t'][full]                  # unfolded / alternative-model truth
            pull = (d['v'][ok] - d['t'][ok]) / d['e'][ok]    # bias in units of the quoted uncertainty
            i = int(np.argmax(np.abs(pull)))
            print('%-16s %-22s %5d %8.3f %8.1f%% %9.2f %9d' %
                  (LAB[obs], desc, full.sum(), np.mean(r), 100 * np.max(np.abs(r - 1)),
                   pull[i], np.flatnonzero(ok)[i] + 1))
            rows.append((obs, tag, np.mean(r), np.max(np.abs(r - 1)), np.abs(pull).max()))
    if rows:
        print('\nlargest |bias| over all bins and variations: %.2f sigma of the quoted uncertainty'
              % max(x[4] for x in rows))
        print('mean unfolded/truth ratio spans %.3f--%.3f'
              % (min(x[2] for x in rows), max(x[2] for x in rows)))


if __name__ == '__main__':
    main()
