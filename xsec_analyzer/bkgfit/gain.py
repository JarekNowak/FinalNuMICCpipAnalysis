"""bkgfit.gain -- expected effect of a control-region constraint on the one-bin total cross section.

Per beam mode b the one-bin total is x = (D - B) * H / S, with D the selected data, B the selected
background (neutrino MC + beam-off + dirt), S the selected signal and H = G / Phi the generated signal
per unit flux. In a flux universe the integrated flux moves with the generated signal, so H stays at
its central value there (this is what makes the extraction's flux term 17% x (1 + B/S)); in every other
universe H = G_u. With Asimov data D = S + B.

The joint Gaussian of z = [control-region fit bins, B_b, S_b, H_b] is built from the same universes,
detector knobs, POT/target and statistical terms as bkgfit.model. The constraint is the Gaussian
posterior of (B, S, H) given the control-region bins, with an optional set of FREE class
normalisations (flat prior), in closed form (generalised least squares / universal kriging):
  Var = V_yy - V_yx Vxx^-1 V_xy + R (Jx^T Vxx^-1 Jx)^-1 R^T,  R = Jy - V_yx Vxx^-1 Jx.
No free parameters is exactly the conditional-covariance constraint (ConstrainedCalculator).
x is propagated linearly.
"""
import numpy as np
from .model import PERIODS, NC, SR, SIGGEN, SYST_AVG, KNOBS, CC0PI, PI0, COSMIC

SIG = [0, 1]


def build(T, periods, fit_bins, detector=True, pot=True, mcstat=True, dvstat=True, decorrelate=(), exclude=None, det='matched'):
    """decorrelate: families whose control-region x signal-region blocks are zeroed (the constraint may not
    act through them). exclude: {family: [universe indices]} or {'detector': [knob indices]} left out
    (leave-one-out closure)."""
    beams = [m for m in ('FHC', 'RHC') if any(p.startswith(m) for p in periods)]
    groups = [[PERIODS.index(p) for p in periods if p.startswith(m)] for m in beams]
    Bf = np.array(fit_bins); nx = len(groups) * len(Bf); ny = 3 * len(groups)
    bkgc = [c for c in range(NC) if c not in SIG]

    def vec(cvP, extP, dirtP, flux=False, cv_for_H=None):
        """cvP [P,C,B] (universe or CV), ext/dirt [P,B] -> z vector."""
        x = np.concatenate([cvP[g][:, :, Bf].sum((0, 1)) + extP[g][:, Bf].sum(0) + dirtP[g][:, Bf].sum(0) for g in groups])
        y = []
        for gi, g in enumerate(groups):
            Bv = cvP[g][:, bkgc, SR].sum() + extP[g][:, SR].sum() + dirtP[g][:, SR].sum()
            Sv = cvP[g][:, SIG, SR].sum()
            Hv = cv_for_H[gi] if flux else cvP[g][:, SIG, SIGGEN].sum()
            y += [Bv, Sv, Hv]
        return x, np.array(y)

    zx0, zy0 = vec(T.cv, T.ext, T.dirt)
    H0 = [T.cv[g][:, SIG, SIGGEN].sum() for g in groups]
    n = nx + ny; V = np.zeros((n, n)); fam = {}
    exclude = exclude or {}
    for s, U in T.univ.items():
        keep = [u for u in range(U.shape[1]) if u not in exclude.get(s, ())]
        nu = len(keep); D = np.zeros((nu, n))
        for iu, u in enumerate(keep):
            cvu = T.cv + U[:, u]
            ux, uy = vec(cvu, T.ext, T.dirt, flux=(s == 'flux'), cv_for_H=H0)
            D[iu] = np.concatenate([ux - zx0, uy - zy0])
        Vs = D.T @ D / (nu if SYST_AVG[s] else 1.0); V += Vs; fam[s] = Vs
    # detector: fractional knob shifts on the neutrino MC of each beam (H unaffected). det='matched' (default)
    # uses the POT-normalised event-matched shifts with the bootstrap noise subtracted (bkgfit.detvar);
    # 'pot' leaves the noise in; 'raw' is the unscaled convention of the sideband macros (+ independent noise).
    if detector:
        Vd = np.zeros((n, n))
        # element descriptors: (template bin, neutrino-MC amount, mode) for every z entry; H entries carry none
        desc = []
        for gi, g in enumerate(groups):
            m = 'run4fhc' if PERIODS[g[0]].startswith('FHC') else 'run4rhc'
            mcx = T.cv[g][:, :, Bf].sum((0, 1))
            desc += [(int(b), float(a), m) for b, a in zip(Bf, mcx)]
        ydesc = []
        for gi, g in enumerate(groups):
            m = 'run4fhc' if PERIODS[g[0]].startswith('FHC') else 'run4rhc'
            ydesc += [(SR, float(T.cv[g][:, bkgc, SR].sum()), m), (SR, float(T.cv[g][:, SIG, SR].sum()), m), (None, 0.0, m)]
        desc = desc + ydesc
        binv = np.array([d[0] if d[0] is not None else 0 for d in desc]); amt = np.array([d[1] for d in desc])
        modev = np.array([d[2] for d in desc]); isH = np.array([d[0] is None for d in desc])
        for k in range(len(KNOBS)):
            if k in exclude.get('detector', ()): continue
            if det == 'raw':
                fk = np.array([T.dv[md][k][b] for b, md in zip(binv, modev)])
            else:
                fk = np.array([T.dvm[md]['f'][k][b] for b, md in zip(binv, modev)])
            d = np.where(isH, 0.0, fk * amt); Vd += np.outer(d, d)
            if det == 'matched':
                for md in ('run4fhc', 'run4rhc'):
                    sel = np.nonzero((modev == md) & ~isH)[0]
                    if len(sel) == 0: continue
                    Cn = T.dvm[md]['noise'][k][np.ix_(binv[sel], binv[sel])]
                    Vd[np.ix_(sel, sel)] -= Cn * np.outer(amt[sel], amt[sel])
            if det == 'raw' and dvstat:
                cvv, cv2 = T.dvraw[modev[0]]['CV']
        if det == 'raw' and dvstat:
            for i in np.nonzero(~isH)[0]:
                md = modev[i]; cvv, cv2 = T.dvraw[md]['CV']; b = binv[i]
                fv = sum((1 + T.dv[md][k][b]) ** 2 * ((T.dvraw[md][KNOBS[k]][1][b] / T.dvraw[md][KNOBS[k]][0][b] ** 2 if T.dvraw[md][KNOBS[k]][0][b] > 0 else 0)
                         + (cv2[b] / cvv[b] ** 2 if cvv[b] > 0 else 0)) for k in range(len(KNOBS)) if k not in exclude.get('detector', ()))
                Vd[i, i] += fv * amt[i] ** 2
        if det == 'matched':
            w, Q = np.linalg.eigh(Vd); Vd = (Q * np.clip(w, 0, None)) @ Q.T
        V += Vd; fam['detector'] = Vd
    if pot:
        m = np.concatenate([np.concatenate([T.cv[g][:, :, Bf].sum((0, 1)) for g in groups]),
                            np.concatenate([[T.cv[g][:, bkgc, SR].sum(), T.cv[g][:, SIG, SR].sum(), 0.0] for g in groups])])
        Vp = (0.02 ** 2 + 0.01 ** 2) * np.outer(m, m); V += Vp; fam['POT+targets'] = Vp
    if mcstat:
        dvar = np.concatenate([np.concatenate([T.w2[g][:, :, Bf].sum((0, 1)) + T.ext2[g][:, Bf].sum(0) + T.dirt2[g][:, Bf].sum(0) for g in groups]),
                               np.concatenate([[T.w2[g][:, bkgc, SR].sum() + T.ext2[g][:, SR].sum() + T.dirt2[g][:, SR].sum(),
                                                T.w2[g][:, SIG, SR].sum(), T.w2[g][:, SIG, SIGGEN].sum()] for g in groups])])
        Vm = np.diag(dvar); V += Vm; fam['MC stat'] = Vm
    for s in decorrelate:
        if s not in fam: continue
        blk = fam[s].copy(); blk[:nx, nx:] = 0; blk[nx:, :nx] = 0
        V += blk - fam[s]; fam[s] = blk
    # class Jacobians for free normalisations: d z / d theta_c
    Jx = {c: np.concatenate([T.cv[g][:, c, Bf].sum(0) for g in groups]) for c in range(NC)}
    Jy = {c: np.concatenate([[T.cv[g][:, c, SR].sum() if c not in SIG else 0.0, 0.0, 0.0] for g in groups]) for c in range(NC)}
    return dict(beams=beams, nx=nx, ny=ny, V=V, fam=fam, zx=zx0, zy=zy0, Jx=Jx, Jy=Jy)


def x_and_grad(zy, gi):
    B, S, H = zy[3 * gi:3 * gi + 3]; D = S + B
    x = (D - B) * H / S
    # D is data (fixed): dx/dB = -H/S, dx/dS = -(D-B)H/S^2, dx/dH = (D-B)/S
    g = np.zeros(len(zy)); g[3 * gi] = -H / S; g[3 * gi + 1] = -(D - B) * H / S ** 2; g[3 * gi + 2] = (D - B) / S
    return x, g, D * (H / S) ** 2


def posterior_yy(G, free_classes=(), use_cr=True):
    nx = G['nx']; V = G['V']
    Vyy = V[nx:, nx:]
    if not use_cr: return Vyy
    Vxx = V[:nx, :nx] + np.diag(np.maximum(G['zx'], 1.0))      # Asimov data Poisson term on the CR bins
    Vyx = V[nx:, :nx]; Vxx_i = np.linalg.inv(Vxx)
    post = Vyy - Vyx @ Vxx_i @ Vyx.T
    groups = [list(c) for c in free_classes]
    if groups:
        Jx = np.array([sum(G['Jx'][c] for c in grp) for grp in groups]).T
        Jy = np.array([sum(G['Jy'][c] for c in grp) for grp in groups]).T
        R = Jy - Vyx @ Vxx_i @ Jx
        post = post + R @ np.linalg.inv(Jx.T @ Vxx_i @ Jx) @ R.T
    return post


def summary(G, free_sets):
    """Relative uncertainty on x per beam: prior, conditional constraint, and each free-parameter set;
    also the prior breakdown by family (the flux term is the check against the extraction)."""
    out = {}
    nx = G['nx']
    for gi, b in enumerate(G['beams']):
        x, g, dstat = x_and_grad(G['zy'], gi)
        row = {}
        for name, Vfam in G['fam'].items():
            row[f'prior {name}'] = float(np.sqrt(g @ Vfam[nx:, nx:] @ g) / x)
        row['prior data stat'] = float(np.sqrt(dstat) / x)
        row['prior total'] = float(np.sqrt(g @ posterior_yy(G, use_cr=False) @ g + dstat) / x)
        row['conditional (A)'] = float(np.sqrt(g @ posterior_yy(G) @ g + dstat) / x)
        for nm, fs in free_sets.items():
            row[f'fit {nm} (B)'] = float(np.sqrt(g @ posterior_yy(G, fs) @ g + dstat) / x)
        B, S, H = G['zy'][3 * gi:3 * gi + 3]
        row['B/S'] = float(B / S)
        out[b] = row
    return out


def protocol_style(G, free_classes):
    """The frozen procedure's propagation (item 4(ii)): the class normalisation is updated in every universe,
    so the flux, cross-section and detector terms are UNCHANGED, and the fitted normalisation uncertainty is
    added as an independent term. The fit uncertainty is the Fisher uncertainty of the class norms from the
    control regions alone."""
    nx = G['nx']; V = G['V']; out = {}
    Vxx = V[:nx, :nx] + np.diag(np.maximum(G['zx'], 1.0)); Vi = np.linalg.inv(Vxx)
    Jx = np.array([sum(G['Jx'][c] for c in grp) for grp in free_classes]).T
    Jy = np.array([sum(G['Jy'][c] for c in grp) for grp in free_classes]).T
    Ct = np.linalg.inv(Jx.T @ Vi @ Jx)
    for gi, b in enumerate(G['beams']):
        x, g, dstat = x_and_grad(G['zy'], gi)
        prior = g @ V[nx:, nx:] @ g + dstat
        fitterm = g @ Jy @ Ct @ Jy.T @ g
        out[b] = dict(prior=float(np.sqrt(prior) / x), with_fit=float(np.sqrt(prior + fitterm) / x),
                      sigma_theta=np.sqrt(np.diag(Ct)).tolist())
    return out


def closure(T, periods, fit_bins, family, indices, free_classes=()):
    """Leave-one-out closure: for each universe (or detector knob) the truth is that variation, the covariance is
    built WITHOUT it, the control-region Asimov data and the signal-region data come from it, and the constrained
    one-bin total is compared with the truth. Returns the prior and posterior biases and pulls per beam."""
    rows = []
    for idx in indices:
        G = build(T, periods, fit_bins, exclude={family: [idx]})
        nx = G['nx']
        # truth vectors from the excluded variation
        if family == 'detector':
            zx = G['zx'].copy(); zy = G['zy'].copy()
            beams = G['beams']; nb = nx // len(beams)
            for gi, bm in enumerate(beams):
                m = 'run4fhc' if bm == 'FHC' else 'run4rhc'; fr = T.dvm[m]['f'][idx]
                g = [PERIODS.index(p) for p in periods if p.startswith(bm)]
                mcx = T.cv[g][:, :, fit_bins].sum((0, 1)); zx[gi * nb:(gi + 1) * nb] += fr[fit_bins] * mcx
                zy[3 * gi] += fr[SR] * T.cv[g][:, [c for c in range(NC) if c not in SIG], SR].sum(); zy[3 * gi + 1] += fr[SR] * T.cv[g][:, SIG, SR].sum()
        else:
            U = T.univ[family][:, idx]; H0 = [T.cv[[PERIODS.index(p) for p in periods if p.startswith(bm)]][:, SIG, SIGGEN].sum() for bm in G['beams']]
            groups = [[PERIODS.index(p) for p in periods if p.startswith(bm)] for bm in G['beams']]
            Bf = np.array(fit_bins); cvu = T.cv + U; bkgc = [c for c in range(NC) if c not in SIG]
            zx = np.concatenate([cvu[g][:, :, Bf].sum((0, 1)) + T.ext[g][:, Bf].sum(0) + T.dirt[g][:, Bf].sum(0) for g in groups])
            zy = np.concatenate([[cvu[g][:, bkgc, SR].sum() + T.ext[g][:, SR].sum() + T.dirt[g][:, SR].sum(), cvu[g][:, SIG, SR].sum(),
                                  H0[gi] if family == 'flux' else cvu[g][:, SIG, SIGGEN].sum()] for gi, g in enumerate(groups)])
        Vxx = G['V'][:nx, :nx] + np.diag(np.maximum(G['zx'], 1.0)); Vyx = G['V'][nx:, :nx]
        K = Vyx @ np.linalg.inv(Vxx)
        post_y = G['zy'] + K @ (zx - G['zx'])
        Vpost = posterior_yy(G, free_classes)
        row = dict(index=idx)
        for gi, bm in enumerate(G['beams']):
            B, S, H = zy[3 * gi:3 * gi + 3]; D = S + B; x_true = H * (D - B) / S if family != 'flux' else H
            x_true = H  # (D - B) = S by construction: the truth is the generated signal per unit flux
            Bp, Sp, Hp = post_y[3 * gi:3 * gi + 3]; x_post = (D - Bp) * Hp / Sp
            B0, S0, H0v = G['zy'][3 * gi:3 * gi + 3]; x_prior = (D - B0) * H0v / S0
            _, gr, dstat = x_and_grad(G['zy'], gi)
            sp = np.sqrt(gr @ Vpost @ gr); s0 = np.sqrt(gr @ G['V'][nx:, nx:] @ gr)
            row[bm] = dict(prior_rel=(x_prior - x_true) / x_true, post_rel=(x_post - x_true) / x_true,
                           prior_pull=(x_prior - x_true) / s0, post_pull=(x_post - x_true) / sp)
        rows.append(row)
    return rows
