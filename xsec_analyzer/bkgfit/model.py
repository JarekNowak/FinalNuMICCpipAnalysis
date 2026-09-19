"""bkgfit.model -- prediction, covariance and test statistic for the joint FHC+RHC control-region fit.

Inputs are the templates of macros/bkgfit_templates.C. Conventions follow the extraction:
  * multisim families (flux, GENIE multisim, re-interaction, RPA/SCC pairs) average over universes,
    unisim knobs sum (configs/ccpi_systcalc_numi.conf);
  * the eight detector knobs are fractional shifts from the Run-4 samples of each beam mode, taken on
    the class-summed prediction (per class they are at the Monte-Carlo-statistical floor), fully
    correlated between run periods of a mode and between the two modes (same physical variation);
  * POT 2% and target count 1% fully correlated over the neutrino MC;
  * MC, beam-off and dirt statistics; data statistics in the combined Neyman-Pearson form.
Class parameters multiply the class templates IN EVERY UNIVERSE, so a fitted normalisation keeps its
correlation with the flux and cross-section terms (the (1+B/S) coherence of the extraction).
"""
import numpy as np
import uproot

PERIODS = ['FHC_R1', 'FHC_R2', 'FHC_R4', 'FHC_R5', 'RHC_R1', 'RHC_R2', 'RHC_R3', 'RHC_R4']
POT = dict(zip(PERIODS, [2.192, 1.268, 2.075, 2.231, 0.6053, 2.591, 5.003, 2.883]))
DATA_PERIODS = ['FHC_R1', 'FHC_R4', 'FHC_R5', 'RHC_R1', 'RHC_R3', 'RHC_R4']
NB, NC = 46, 10
SR, SIGGEN, MULTI, COSMIC = 0, 1, 44, 45
CC0PI = list(range(2, 23)); PI0 = list(range(23, 44))
CTH = [-1, 0.0, 0.45, 0.65, 0.8, 0.9, 0.95, 1.0]; PMU = [0.15, 0.35, 0.75, np.inf]
TOT_CC0PI, TOT_PI0 = 46, 47                                                                     # virtual total bins
CLASS_NAMES = ['signal', 'CC0pi', 'pi0', 'multipi', 'other']
SYST_AVG = {'flux': True, 'xsec_multi': True, 'reint': True, 'xsec_AxFFCCQEshape': False, 'xsec_DecayAngMEC': False,
            'xsec_NormCCCOH': False, 'xsec_NormNCCOH': False, 'xsec_RPA_CCQE': True, 'xsec_ThetaDelta2NRad': False,
            'xsec_Theta_Delta2Npi': False, 'xsec_VecFFCCQEshape': False, 'xsec_XSecShape_CCMEC': False,
            'xsec_xsr_scc_Fa3_SCC': True, 'xsec_xsr_scc_Fv3_SCC': True}
KNOBS = ['LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']

def cls(c, nubar):
    return 2 * c + nubar

def bin_label(i):
    if i == SR: return 'SR'
    if i == SIGGEN: return 'SIGGEN'
    if i == MULTI: return 'MULTI'
    if i == COSMIC: return 'COSMIC'
    if i == TOT_CC0PI: return 'CC0pi total'
    if i == TOT_PI0: return 'PI0 total'
    reg = 'CC0pi' if i in CC0PI else 'PI0'
    k = i - (2 if reg == 'CC0pi' else 23); ic, ip = divmod(k, 3)
    return f'{reg} cth[{CTH[ic]},{CTH[ic + 1]}) pmu[{PMU[ip]},{PMU[ip + 1]})'


class Templates:
    """All template arrays in memory."""
    def __init__(self, path='/data/uboone/processed/bkgfit/templates.root'):
        f = uproot.open(path)
        keys = set(k.split(';')[0] for k in f.keys())
        self.cv = np.array([f[f'cv_{p}'].values().reshape(NC, NB) for p in PERIODS])          # [P,C,B]
        self.w2 = np.array([f[f'w2_{p}'].values().reshape(NC, NB) for p in PERIODS])
        self.ext = np.array([f[f'ext_{p}'].values() for p in PERIODS])                          # [P,B]
        self.ext2 = np.array([f[f'ext2_{p}'].values() for p in PERIODS])
        self.dirt = np.array([f[f'dirt_{p}'].values() for p in PERIODS])
        self.dirt2 = np.array([f[f'dirt2_{p}'].values() for p in PERIODS])
        self.data = np.array([f[f'data_{p}'].values() for p in PERIODS])
        self.univ = {}                                                                          # s -> [P,U,C,B] deltas
        for s in SYST_AVG:
            if f'u_{s}_{PERIODS[0]}' not in keys: continue
            arr = []
            for pi, p in enumerate(PERIODS):
                u = f[f'u_{s}_{p}'].values()                                                    # [U, C*B]
                arr.append(u.reshape(u.shape[0], NC, NB) - self.cv[pi][None])
            nu = min(a.shape[0] for a in arr)
            self.univ[s] = np.array([a[:nu] for a in arr])
        self.dv = {}; self.dvraw = {}
        for m in ('run4fhc', 'run4rhc'):
            self.dvraw[m] = {k: (f[f'dv_{m}_{k}'].values().reshape(NC, NB).sum(0), f[f'dv2_{m}_{k}'].values().reshape(NC, NB).sum(0))
                             for k in ['CV'] + KNOBS}
            cv = f[f'dv_{m}_CV'].values().reshape(NC, NB).sum(0)
            fr = []
            for k in KNOBS:
                alt = f[f'dv_{m}_{k}'].values().reshape(NC, NB).sum(0)
                fr.append(np.divide(alt - cv, cv, out=np.zeros(NB), where=cv > 0))
            fr = np.array(fr); fr[:, SIGGEN] = 0.0                                              # truth-level count
            self.dv[m] = fr                                                                     # [K,B]
        self._add_virtual([CC0PI, PI0])
        # event-matched, POT-normalised shifts and their bootstrap noise (bkgfit.detvar); fall back to nothing
        import os
        dvp = os.path.join(os.path.dirname(path), 'detvar_matched.npz')
        self.dvm = {}
        if os.path.exists(dvp):
            z = np.load(dvp)
            for m in ('run4fhc', 'run4rhc'):
                self.dvm[m] = dict(f=z[f'{m}_f'], noise=z[f'{m}_noise'])

    def _add_virtual(self, groups):
        """Append virtual bins that are sums of template bins (region totals), consistently in every array."""
        cat = lambda a, ax=-1: np.concatenate([a] + [a.take(g, axis=ax).sum(axis=ax, keepdims=True) for g in groups], axis=ax)
        self.cv, self.w2 = cat(self.cv), cat(self.w2)
        self.ext, self.ext2, self.dirt, self.dirt2, self.data = cat(self.ext), cat(self.ext2), cat(self.dirt), cat(self.dirt2), cat(self.data)
        self.univ = {s: cat(u) for s, u in self.univ.items()}
        for m in self.dvraw:
            self.dvraw[m] = {k: (cat(v[0]), cat(v[1])) for k, v in self.dvraw[m].items()}
            cv = self.dvraw[m]['CV'][0]
            fr = np.array([np.divide(self.dvraw[m][k][0] - cv, cv, out=np.zeros_like(cv), where=cv > 0) for k in KNOBS])
            fr[:, SIGGEN] = 0.0; self.dv[m] = fr


class Model:
    """Prediction and covariance for a choice of periods, bins and class parameters.

    params: dict name -> list of class indices it scales (e.g. {'CC0pi': [2, 3]}).
    extra:  dict name -> list of periods whose WHOLE prediction (MC+EXT+dirt) it scales (exposure).
    merge:  sum the run periods of a beam mode into one set of bins.
    """
    def __init__(self, T, periods, bins, params=None, extra=None, merge=False, logdet=False,
                 stat='pearson', detector=True, pot=True, mcstat=True, dvstat=True, cov='nominal', det='matched', systs=None):
        """cov: 'nominal' holds the systematic covariance at theta = 1 during the fit (default: with a
        theta-dependent covariance a poor shape fit is 'improved' by inflating the covariance and the parameters run
        to their limits, which is what the beam-on data do); 'dependent' re-evaluates it at theta (then logdet=True
        is the consistent likelihood); an array fixes it at that theta (for iteration).
        stat: 'pearson' (default) data variance = prediction at the same theta as the systematic covariance; 'cnp' the
        combined Neyman-Pearson form, which depends on the data and leaves a +0.02 to +0.03 sigma pull bias in the 2D
        fits (100 000-toy ensembles, 2026-09-17, bkgfit.studies ensembles); 'neyman' data variance.
        det: 'matched' (default) POT-normalised event-matched knob shifts with the bootstrap statistical noise
        subtracted (f f^T - C, projected to PSD); 'pot' POT-normalised shifts, noise left in; 'raw' the unscaled
        count ratios of the sideband macros (with dvstat: independent-sample noise added)."""
        self.T, self.periods, self.bins = T, list(periods), list(bins)
        self.params = params or {}; self.extra = extra or {}
        self.names = list(self.params) + list(self.extra)
        self.merge, self.logdet, self.stat = merge, logdet, stat
        self.detector, self.pot, self.mcstat, self.dvstat, self.cov, self.det = detector, pot, mcstat, dvstat, cov, det
        self.pidx = [PERIODS.index(p) for p in self.periods]
        if merge:
            self.groups = [[p for p in self.pidx if PERIODS[p].startswith(m)] for m in ('FHC', 'RHC')]
            self.groups = [g for g in self.groups if g]
        else:
            self.groups = [[p] for p in self.pidx]
        self.nbin = len(self.groups) * len(self.bins)
        self.labels = [f"{'+'.join(PERIODS[p] for p in g)} {bin_label(b)}" for g in self.groups for b in self.bins]
        B = np.array(self.bins)
        # per-group arrays [G, C, b]
        self.cv = np.array([T.cv[g][:, :, B].sum(0) for g in self.groups])
        self.w2 = np.array([T.w2[g][:, :, B].sum(0) for g in self.groups])
        self.fixed = np.array([(T.ext[g][:, B] + T.dirt[g][:, B]).sum(0) for g in self.groups])
        self.fixed2 = np.array([(T.ext2[g][:, B] + T.dirt2[g][:, B]).sum(0) for g in self.groups])
        self.data = np.array([T.data[g][:, B].sum(0) for g in self.groups]).ravel()
        self.univ = {s: np.array([u[g][:, :, :, B].sum(0) for g in self.groups]).transpose(1, 0, 2, 3)   # [U,G,C,b]
                     for s, u in T.univ.items() if systs is None or s in systs}
        self.mode_of_group = ['run4fhc' if PERIODS[g[0]].startswith('FHC') else 'run4rhc' for g in self.groups]
        self.dvfrac = np.array([T.dv[m][:, B] for m in self.mode_of_group])                    # [G,K,b]
        # variance of each knob's fractional shift from the finite detector-variation samples, summed over knobs [G,b]
        dvv = {}
        for m in set(self.mode_of_group):
            cvv, cv2 = T.dvraw[m]['CV']; rel2 = np.divide(cv2, cvv ** 2, out=np.zeros_like(cvv), where=cvv > 0); tot = np.zeros(len(cvv))
            for k, kn in enumerate(KNOBS):
                av, a2 = T.dvraw[m][kn]
                tot += (1 + T.dv[m][k]) ** 2 * (np.divide(a2, av ** 2, out=np.zeros_like(av), where=av > 0) + rel2)
            dvv[m] = tot[B]
        self.dvfracvar = np.array([dvv[m] for m in self.mode_of_group])
        if det in ('matched', 'pot'):
            self.dvfrac = np.array([T.dvm[m]['f'][:, B] for m in self.mode_of_group])        # [G,K,b]
            self.dvnoise = {m: T.dvm[m]['noise'][:, B][:, :, B] for m in set(self.mode_of_group)}  # [K,b,b]
        # exposure parameters act per group
        self.extra_mask = {n: np.array([all(PERIODS[p] in per for p in g) for g in self.groups]) for n, per in self.extra.items()}

    # ---- class and exposure scale factors
    def scales(self, theta):
        """Class scale per group [G, C] and exposure scale per group [G]. A parameter is a list of classes (acting in
        both beams) or a dict {'classes': [...], 'beam': 'FHC'|'RHC'} (acting in that beam only)."""
        th = dict(zip(self.names, theta))
        cs = np.ones((len(self.groups), NC))
        for n, spec in self.params.items():
            cl, beam = (spec['classes'], spec['beam']) if isinstance(spec, dict) else (spec, None)
            for gi, g in enumerate(self.groups):
                if beam is None or PERIODS[g[0]].startswith(beam):
                    cs[gi, cl] = th[n]
        gs = np.ones(len(self.groups))
        for n, mask in self.extra_mask.items():
            gs = np.where(mask, gs * th[n], gs)
        return cs, gs

    def mc(self, theta):
        cs, gs = self.scales(theta)
        return np.einsum('gcb,gc,g->gb', self.cv, cs, gs)

    def predict(self, theta, split=False):
        cs, gs = self.scales(theta)
        mc = np.einsum('gcb,gc,g->gb', self.cv, cs, gs); fx = self.fixed * gs[:, None]
        return (mc, fx) if split else (mc + fx).ravel()

    def gram(self):
        """Sum over systematic families of (1/N) sum_u U_u,c,i U_u,d,j, shape [C, n, C, n]. The universe part of the
        covariance is quadratic in the class scales, so this is computed once per model."""
        if not hasattr(self, '_gram'):
            n = self.nbin; G = np.zeros((NC, n, NC, n))
            for s, U in self.univ.items():
                A = U.transpose(0, 2, 1, 3).reshape(U.shape[0], NC, n)              # [U, C, n] (groups flattened)
                G += np.tensordot(A, A, axes=(0, 0)) / (U.shape[0] if SYST_AVG[s] else 1.0)
            self._gram = G
        return self._gram

    def cov_syst(self, theta):
        cs, gs = self.scales(theta)
        n = self.nbin
        Sb = np.repeat(cs * gs[:, None], len(self.bins), axis=0).T                  # [C, n]: class scale of each bin
        V = np.einsum('ci,cidj,dj->ij', Sb, self.gram(), Sb, optimize=True)
        mc = self.mc(theta)
        if self.detector: V += self.cov_detector(mc)
        if self.pot:
            m = mc.ravel(); V += (0.02 ** 2 + 0.01 ** 2) * np.outer(m, m)
        if self.mcstat:
            V += np.diag((np.einsum('gcb,gc,g->gb', self.w2, cs ** 2, gs ** 2) + self.fixed2 * gs[:, None] ** 2).ravel())
        return V

    def cov_stat(self, mu, data):
        if self.stat == 'pearson': return np.diag(np.maximum(mu, 1e-9))
        if self.stat == 'neyman': return np.diag(np.maximum(data, 1.0))
        v = np.where(data > 0, 3.0 / (1.0 / np.maximum(data, 1e-9) + 2.0 / np.maximum(mu, 1e-9)), mu / 2.0)
        return np.diag(np.maximum(v, 1e-9))

    def chi2(self, theta, data=None):
        data = self.data if data is None else data
        mu = self.predict(theta); r = data - mu
        # the CNP variance grows with mu(theta): evaluated at theta it lets the fit lower chi2 by raising the
        # prediction (a Pearson-type upward bias, ~1 sigma over 86 bins in toys). Unless the covariance is fully
        # theta-dependent it is therefore evaluated at the SAME theta as the systematic covariance.
        if isinstance(self.cov, str) and self.cov == 'dependent':
            mu_s = mu
        else:
            mu_s = self.predict(np.ones(len(self.names)) if isinstance(self.cov, str) else np.array(self.cov))
        V = self.cov_for_fit(theta) + self.cov_stat(mu_s, data)
        c, low = np.linalg.slogdet(V)
        x = np.linalg.solve(V, r)
        return float(r @ x + (low if self.logdet else 0.0))

    def cov_detector(self, mc):
        n = self.nbin; nb = len(self.bins)
        D = (self.dvfrac * mc[:, None, :]).transpose(1, 0, 2).reshape(len(KNOBS), n)
        V = D.T @ D
        if self.det == 'raw':
            if self.dvstat: V += np.diag((self.dvfracvar * mc ** 2).ravel())
            return V
        if self.det == 'matched':
            for g, mg in enumerate(self.mode_of_group):
                for h, mh in enumerate(self.mode_of_group):
                    if mg != mh: continue                                    # the two modes use independent samples
                    V[g * nb:(g + 1) * nb, h * nb:(h + 1) * nb] -= self.dvnoise[mg].sum(0) * np.outer(mc[g], mc[h])
            w, Q = np.linalg.eigh(V); V = (Q * np.clip(w, 0, None)) @ Q.T
        return V

    def cov_for_fit(self, theta):
        if isinstance(self.cov, str) and self.cov == 'dependent': return self.cov_syst(theta)
        key = 'nominal' if isinstance(self.cov, str) else tuple(np.round(self.cov, 12))
        if getattr(self, '_covkey', None) != key:
            self._covfix = self.cov_syst(np.ones(len(self.names)) if key == 'nominal' else np.array(self.cov)); self._covkey = key
        return self._covfix

    def fisher(self, theta, eps=1e-4):
        """Fisher information of the parameters at theta for Asimov data (Gaussian, V fixed at theta)."""
        mu0 = self.predict(theta); V = self.cov_syst(theta) + self.cov_stat(mu0, mu0); Vi = np.linalg.inv(V)
        J = []
        for k in range(len(self.names)):
            t = np.array(theta, float); t[k] += eps; J.append((self.predict(t) - mu0) / eps)
        J = np.array(J)
        return J @ Vi @ J.T
