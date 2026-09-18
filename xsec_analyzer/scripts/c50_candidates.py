# Candidates of the 0.50-criterion study (2026-09-16) sent to full extraction, chosen from
# logs/c50/screen_{fhc,rhc}_{incl,1p}.log: 1D binnings that pass 0.50 (not 0.68) in BOTH FHC
# and RHC and were not already built in the 2026-09-05 re-check, plus the proton angles
# (branch already in processed/w). theta_pipr and the 2D pairs wait for the user's decision.
PI = 3.1416
CANDS = [
    # inclusive
    dict(name='costhmu6', pfx='ccpi', xvar='costhmu', xedges=[-1, 0.3, 0.55, 0.7, 0.8, 0.9, 1], xopen=False, note='min diag 0.64/0.63'),
    dict(name='costhmu7', pfx='ccpi', xvar='costhmu', xedges=[-1, 0.2, 0.45, 0.6, 0.7, 0.8, 0.9, 1], xopen=False, note='min diag 0.61/0.61'),
    dict(name='thetamu7', pfx='ccpi', xvar='thetamu', xedges=[0, 0.35, 0.5, 0.65, 0.8, 1.0, 1.3, 3.15], xopen=False, note='min diag 0.63/0.62'),
    dict(name='thetamu8', pfx='ccpi', xvar='thetamu', xedges=[0, 0.3, 0.42, 0.55, 0.68, 0.84, 1.05, 1.4, 3.15], xopen=False, note='min diag 0.61/0.61'),
    dict(name='ppi3b', pfx='ccpi', xvar='ppi', xedges=[0.175, 0.205, 0.26, 1.0], xopen=True, note='min diag 0.51/0.53'),
    # proton-tagged
    dict(name='dpt3', pfx='ccpi1p', xvar='dpt', xedges=[0, 0.225, 0.45, 2.5], xopen=True, note='min diag 0.59/0.55'),
    dict(name='dalphat3', pfx='ccpi1p', xvar='dalphat', xedges=[0, 85, 150, 180], xopen=False, note='min diag 0.52/0.50'),
    dict(name='Wpipr2c50', pfx='ccpi1p', xvar='Wpipr', xedges=[1.08, 1.19, 2.90], xopen=True, note='min diag 0.56/0.56'),
    dict(name='Whad2', pfx='ccpi1p', xvar='Whad', xedges=[0, 1.18, 2.74], xopen=True, note='min diag 0.71/0.72, passes 0.68'),
    dict(name='costhp4', pfx='ccpi1p', xvar='costhp', xedges=[-1, 0.275, 0.575, 0.825, 1], xopen=False, note='min diag 0.73/0.66'),
    dict(name='costhp6', pfx='ccpi1p', xvar='costhp', xedges=[-1, 0.275, 0.575, 0.7, 0.8, 0.9, 1], xopen=False, note='min diag 0.58/0.55'),
    dict(name='thetap5', pfx='ccpi1p', xvar='thetap', xedges=[0, 0.377, 0.691, 0.942, 1.225, PI], xopen=False, note='min diag 0.67/0.64'),
    dict(name='thetap7', pfx='ccpi1p', xvar='thetap', xedges=[0, 0.283, 0.471, 0.659, 0.816, 1.005, 1.287, PI], xopen=False, note='min diag 0.56/0.55'),
    # 2026-09-18: the released four-bin theta_p is rank-deficient in RHC (s_min/s1 = 0.002) although its
    # diagonal is fine, so the RHC-optimal schemes are tested as well. Edges from the maximin scan run
    # on the RHC chain (macros/tki_binning_scan.C with MODE="rhc"), which the earlier screen never did.
    # Built for all three configurations: a per-mode binning would be a last resort.
    dict(name='thetap3r', pfx='ccpi1p', xvar='thetap', xedges=[0, 0.785, 1.539, PI], xopen=False,
         cfgs=('fhc5', 'rhcfull', 'comb'), note='RHC maximin, min diag 0.78 RHC'),
    dict(name='thetap4r', pfx='ccpi1p', xvar='thetap', xedges=[0, 0.659, 1.068, 1.539, PI], xopen=False,
         cfgs=('fhc5', 'rhcfull', 'comb'), note='RHC maximin, min diag 0.73 RHC'),
]
