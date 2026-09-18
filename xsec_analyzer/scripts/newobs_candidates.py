# Observables added on the user's decision of 2026-09-16 (theta_p chosen over cos theta_p, theta_pipr
# reprocessed in, and ALL screened 2D pairs built). Edges MUST match
# generator_predictions/gen2d/obs_ext.h (checked by scripts/check_ext_bins.py).
#  * proton-tagged 1D, LIVE family (configs/, processed/): binnings that meet the 0.68 criterion
#    (0.50 is study-only) -- theta_p 4 bins (FHC/RHC min diagonal 0.73/0.68, at the limit like
#    cos theta_mu), theta_pipr 3 bins (0.77/0.75).
#  * 2D pairs, STUDY family (configs/d2/, processed/rebuild_2d/): the best-screened grid of each
#    pair, including those that fail both criteria, for a complete comparison.
PI = 3.1416
# theta_p edges: the RHC-derived maximin scheme adopted 2026-09-18. The FHC-derived edges it replaces
# ([0, 0.597, 0.942, 1.287, PI]) condition three to four times worse in both usable configurations
# (FHC 0.082 vs 0.383, COMB 0.125 vs 0.408); see report/newobs_note.tex, section 3.5.
P1 = 'proton-tagged'
CANDS = [
    dict(name='thetap', pfx='ccpi1p', xvar='thetap', xedges=[0, 0.659, 1.068, 1.539, PI], xopen=False,
         subdir='', suffix='', univdir='/data/uboone/processed', cfgs=('fhc5', 'rhcfull', 'comb'), gen_obs='thetap',
         note='CC1mu1pi1p proton polar angle theta_p about the neutrino direction, 4 bins (maximin scan, min diag 0.73/0.68).'),
    dict(name='thpipr', pfx='ccpi1p', xvar='thpipr', xedges=[0, 1.193, 2.010, PI], xopen=False,
         subdir='', suffix='', univdir='/data/uboone/processed', cfgs=('fhc5', 'rhcfull', 'comb'), gen_obs='thpipr',
         note='CC1mu1pi1p pion-proton opening angle theta_pipr, 3 bins (maximin scan, min diag 0.77/0.75).'),
]
D2 = dict(subdir='d2', suffix='_2d', univdir='/data/uboone/processed/rebuild_2d', cfgs=('fhc5', 'rhcfull', 'comb'))
for nm, pfx, xv, xe, xo, yv, ye, yo, ylo, diag in [
    ('costhpi_costhmu', 'ccpi', 'costhpi', [-1, 0.35, 1], False, 'costhmu', [-1, 0.65, 0.85, 1], False, False, '0.69/0.68'),
    ('costhmu_pmu', 'ccpi', 'costhmu', [-1, 0.8, 1], False, 'pmu', [0.15, 0.55, 3.0], True, False, '0.70/0.68'),
    ('costhpi_thmupi', 'ccpi', 'costhpi', [-1, 0.35, 1], False, 'thmupi', [0, 0.85, 1.5, 2.6], False, False, '0.66/0.60'),
    ('costhpi_ppi', 'ccpi', 'costhpi', [-1, 0.35, 1], False, 'ppi', [0.175, 0.205, 1.0], True, False, '0.53/0.51'),
    ('thetap_dpt', 'ccpi1p', 'thetap', [0, 0.8, PI], False, 'dpt', [0, 0.3, 2.5], True, True, '0.53/0.52'),
    ('thpipr_dpt', 'ccpi1p', 'thpipr', [0, 1.2, PI], False, 'dpt', [0, 0.3, 2.5], True, True, '0.56/0.55'),
    ('thetap_Wpipr', 'ccpi1p', 'thetap', [0, 0.8, PI], False, 'Wpipr', [1.08, 1.19, 2.9], True, True, '0.40/0.39'),
    ('thpipr_Wpipr', 'ccpi1p', 'thpipr', [0, 1.2, PI], False, 'Wpipr', [1.08, 1.19, 2.9], True, True, '0.21/0.20'),
]:
    CANDS.append(dict(name=nm, pfx=pfx, xvar=xv, xedges=xe, xopen=xo, yvar=yv, yedges=ye, yopen=yo, ylow_open=ylo,
                      gen_obs=nm, note=f'2D study {nm} (X | Y), screened min diagonal FHC/RHC {diag}.', **D2))
