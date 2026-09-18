"""farsb.py -- far sidebands: control regions made purer in background by one extra cut.

The four control regions are defined in CC1mu1piXp.cxx (sb_cc0pi, sb_pi0, sb_multipi, sb_cosmic).
Each still contains some signal, and signal in a control region is what makes a background
constraint bite its own tail: the fit lowers a background component, the signal region loses the
same events twice. A "far" sideband adds one cut that pushes the region away from the signal
topology, trading statistics for a smaller signal fraction. Whether that trade is worth taking is
decided by three numbers per variant, which this script measures:

  signal fraction   how much of the region is signal (the reason for the exercise);
  retained events   what the extra cut costs in statistics;
  composition       the background mix by class -- a far sideband that constrains a DIFFERENT
                    mixture than the near one is not a purer version of it, it is another region.

Everything is at data exposure, with the per-period MC scales, beam-off gate ratios and dirt
scaling of macros/bkgfit_templates.C (kept in sync by hand; see sb_samples).

  python3 -m bkgfit.farsb [--quick]      (--quick: first MC file per period only)
"""
import sys, os
import numpy as np
import uproot, awkward as ak

S = 'CC1mu1piXp'
P = '/data/uboone/processed/sb_pi0/'
BO = '/data/uboone/processed/beamon_skim/'
PNAME = ['FHC_R1', 'FHC_R2', 'FHC_R4', 'FHC_R5', 'RHC_R1', 'RHC_R2', 'RHC_R3', 'RHC_R4']
POT = [3.283, 1.268, 2.075, 2.231, 0.6053, 2.591, 5.003, 2.883]
G1, G3, G4, G5, OCCX = 4582248.27, 32649128.65, 34831148.625, 19256341.475, 0.98
DIRT = P + 'xsec-ana-prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root'


def samples():
    """(mc, ext, data) per period: lists of (path, scale). Scales are bkgfit_templates.C's."""
    mc = [[] for _ in PNAME]; ext = [[] for _ in PNAME]; data = [[] for _ in PNAME]
    f = lambda n: P + 'xsec-ana-' + n + '.root'
    mc[0] = [(f('Run1_fhc_new_numi_flux_fhc_pandora_ntuple'), 0.14101)]
    mc[1] = [(f('Run2_fhc_new_numi_flux_fhc_pandora_ntuple'), 0.05085)]
    mc[2] = [(f('Run4_fhc_new_numi_flux_fhc_pandora_ntuple'), 0.07323)]
    mc[3] = [(f('reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc'), 0.11560)]
    mc[4] = [(f('Run1_rhc_new_numi_flux_rhc_pandora_ntuple'), 0.06728)]
    mc[5] = [(f('Run2_rhc_new_numi_flux_rhc_pandora_ntuple'), 0.04478)]
    mc[6] = [(f(f'Run3_rhc_new_numi_flux_rhc_pandora_ntuple_{s}'), 0.09066) for s in 'aa ab ac ad ae'.split()]
    mc[7] = [(f(f'{s}_rhc_new_numi_flux_rhc_pandora_ntuple'), 0.08847) for s in ('Run4a', 'Run4b', 'Run4c')]
    R4 = ['numi_pelee_ntuple_beam_off_run4a_rhc_ana.root', 'numi_pelee_ntuple_beam_off_run4b_rhc_ana.root',
          'numi_pelee_ntuple_beam_off_run4c_fhc_ana.root', 'numi_pelee_ntuple_beam_off_run4d_fhc_ana.root']
    E = P + 'xsec-ana-'
    ext[0] = [(E + 'neutrinoselection_filt_run1_beamoff.root', OCCX * 9846635. / G1)]
    ext[2] = [(E + r, OCCX * 4131149. / G4) for r in R4]
    ext[3] = [(E + 'numi_pelee_ntuple_beam_off_run5_fhc_ana.root', OCCX * 5154196. / G5)]
    ext[4] = [(E + 'neutrinoselection_filt_run1_beamoff.root', OCCX * 1458253. / G1)]
    ext[6] = [(E + 'neutrinoselection_filt_run3b_beamoff.root', OCCX * 10349610. / G3)]
    ext[7] = [(E + r, OCCX * 6304167. / G4) for r in R4]
    data[0] = [BO + 'xsec-ana-beamon_fhc_run1.root']
    data[2] = [BO + 'xsec-ana-beamon_fhc_run4c.root', BO + 'xsec-ana-beamon_fhc_run4d.root']
    data[3] = [BO + 'xsec-ana-beamon_fhc_run5.root']
    data[4] = [BO + 'xsec-ana-beamon_rhc_run1.root']
    data[6] = [BO + 'xsec-ana-beamon_rhc_run3b.root']
    data[7] = [BO + 'xsec-ana-beamon_rhc_run4a.root', BO + 'xsec-ana-beamon_rhc_run4b.root']
    return mc, ext, data


# Each variant is (region, label, predicate). The near definitions are the selection's own flags;
# the far ones add exactly one cut, so the comparison isolates that cut.
def variants():
    v = []
    add = lambda r, lab, fn: v.append((r, lab, fn))
    # sb_nnonproton counts primary tracks with LLR > 0.1 INCLUDING the muon, so == 1 means the muon is
    # the only non-proton-like track; sb_nprimtrk >= 2 always in these regions (muon + at least one
    # proton candidate). pion_number_reco saturates at 2, so ">= 3 pions" selects nothing and the
    # multi-pi region has to be pushed on the track counts instead.
    add('cc0pi', 'near (sb_cc0pi)', lambda d: d['cc0pi'])
    add('cc0pi', 'far: muon is only non-proton', lambda d: d['cc0pi'] & (d['nnonproton'] == 1))
    add('cc0pi', 'far: + exactly one proton', lambda d: d['cc0pi'] & (d['nnonproton'] == 1) & (d['nprimtrk'] == 2))
    add('pi0', 'near (sb_pi0)', lambda d: d['pi0'])
    add('pi0', 'far: >=2 showers', lambda d: d['pi0'] & (d['nshr'] >= 2))
    add('pi0', 'far: >=3 showers', lambda d: d['pi0'] & (d['nshr'] >= 3))
    add('multipi', 'near (sb_multipi)', lambda d: d['multipi'])
    # ">=3 non-proton tracks" is implied by sb_multipi (muon + two pion candidates) and returned the
    # region unchanged, so the multi-pi region can only be pushed on the total track count.
    add('multipi', 'far: >=4 primary tracks', lambda d: d['multipi'] & (d['nprimtrk'] >= 4))
    add('cosmic', 'near (sb_cosmic)', lambda d: d['cosmic'])
    add('cosmic', 'far: theta_mupi > 2.9', lambda d: d['cosmic'] & (d['thmupi'] > 2.9))
    add('cosmic', 'far: theta_mupi > 3.0', lambda d: d['cosmic'] & (d['thmupi'] > 3.0))
    return v


CLASSES = ['signal', 'CC0pi', 'pi0', 'multipi', 'other']


def klass(cat, npi, npi0):
    """The five classes of bkgfit_templates.C, by true final state."""
    c = np.full(len(cat), 4)
    c[(cat == 3) & (npi == 0) & (npi0 == 0)] = 1
    c[((cat == 3) | (cat == 5)) & (npi0 >= 1)] = 2
    c[(cat == 3) & (npi >= 2) & (npi0 == 0)] = 3
    c[cat == 0] = 0
    return c


def safe(w):
    w = np.asarray(w, float)
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.0)


BASE = [f'{S}_sb_cc0pi', f'{S}_sb_pi0', f'{S}_sb_multipi', f'{S}_sb_cosmic',
        f'{S}_sb_nprimtrk', f'{S}_sb_nnonproton', f'{S}_pion_number_reco',
        f'{S}_mu_pi_opening_angle', 'pfp_generation_v', 'trk_score_v']
MCB = [f'{S}_EventCategory', f'{S}_mc_n_threshold_pionpm', f'{S}_mc_n_threshold_pion0', 'mc_nu_pdg',
       'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']


def read(path, is_mc):
    """The variables the variants need, plus the class and the CV weight."""
    want = BASE + (MCB if is_mc else [])
    with uproot.open(path) as f:
        t = f['stv_tree']
        have = {k.split(';')[0] for k in t.keys()}
        a = t.arrays([b for b in want if b in have], library='ak')
    g = lambda b, dflt: np.asarray(a[b]) if b in a.fields else np.full(len(a[a.fields[0]]), dflt)
    gen, score = a['pfp_generation_v'], a['trk_score_v']
    nshr = np.asarray(ak.sum((gen == 2) & (score < 0.5), axis=1))
    d = dict(cc0pi=g(f'{S}_sb_cc0pi', False).astype(bool), pi0=g(f'{S}_sb_pi0', False).astype(bool),
             multipi=g(f'{S}_sb_multipi', False).astype(bool), cosmic=g(f'{S}_sb_cosmic', False).astype(bool),
             nprimtrk=g(f'{S}_sb_nprimtrk', -1), nnonproton=g(f'{S}_sb_nnonproton', -1),
             npi=g(f'{S}_pion_number_reco', -1), thmupi=g(f'{S}_mu_pi_opening_angle', -9.), nshr=nshr)
    if is_mc:
        d['cls'] = klass(g(f'{S}_EventCategory', -1), g(f'{S}_mc_n_threshold_pionpm', 0),
                         g(f'{S}_mc_n_threshold_pion0', 0))
        d['w'] = safe(g('tuned_cv_weight', 1.) * g('ppfx_cv_weight', 1.) * g('normalisation_weight', 1.))
    else:
        d['cls'] = np.full(len(nshr), 4); d['w'] = np.ones(len(nshr))
    return d


def main():
    quick = '--quick' in sys.argv
    mc, ext, data = samples()
    V = variants()
    nu = np.zeros((len(V), len(CLASSES)))   # neutrino MC by class
    eo = np.zeros(len(V))                   # beam-off
    dt = np.zeros(len(V))                   # dirt
    ob = np.zeros(len(V))                   # beam-on data
    for p, name in enumerate(PNAME):
        files = [(f, s, 'mc') for f, s in (mc[p][:1] if quick else mc[p])]
        files += [(f, s, 'ext') for f, s in ext[p]]
        sdirt = 0.092402 * 0.65 * POT[p] / 8.857 if p < 4 else 0.071666 * 0.65 * POT[p] / 11.082
        files += [(DIRT, sdirt, 'dirt')] + [(f, 1.0, 'data') for f in data[p]]
        for path, scale, kind in files:
            if not os.path.exists(path):
                print(f'   MISSING {os.path.basename(path)}'); continue
            d = read(path, kind in ('mc', 'dirt'))
            for i, (_, _, fn) in enumerate(V):
                m = fn(d)
                if kind == 'mc':
                    w = d['w'][m] * scale
                    for c in range(len(CLASSES)):
                        nu[i, c] += w[d['cls'][m] == c].sum()
                elif kind == 'ext': eo[i] += m.sum() * scale
                elif kind == 'dirt': dt[i] += (d['w'][m] * scale).sum()
                else: ob[i] += m.sum()
        print(f'== {name} done', flush=True)

    tot = nu.sum(axis=1) + eo + dt
    print('\n%-28s %9s %8s %8s %9s %8s   %s' %
          ('variant', 'pred', 'data', 'd/p', 'signal %', 'kept %', 'background composition [%]'))
    for r in ('cc0pi', 'pi0', 'multipi', 'cosmic'):
        idx = [i for i, (rr, _, _) in enumerate(V) if rr == r]
        near = idx[0]
        for i in idx:
            bkg = tot[i] - nu[i, 0]
            comp = [100 * nu[i, c] / bkg for c in range(1, len(CLASSES))] + [100 * (eo[i] + dt[i]) / bkg]
            print('%-28s %9.1f %8.0f %8.2f %8.1f%% %7.1f%%   %s' % (
                V[i][1], tot[i], ob[i], ob[i] / tot[i] if tot[i] else 0,
                100 * nu[i, 0] / tot[i] if tot[i] else 0,
                100 * tot[i] / tot[near] if tot[near] else 0,
                ' '.join(f'{n}={v:.0f}' for n, v in zip(CLASSES[1:] + ['EXT+dirt'], comp))))
        print()


if __name__ == '__main__':
    main()
