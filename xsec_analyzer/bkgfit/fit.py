"""bkgfit.fit -- minimisation, intervals, tests, conditional constraint and MCMC for bkgfit.model.Model."""
import numpy as np
from iminuit import Minuit
from scipy import stats


def fit(model, data=None, start=None, minos=False, fix=None):
    """Minuit2 MIGRAD + HESSE (+ MINOS). Returns dict with values, errors, covariance, chi2."""
    n = len(model.names)
    if n == 0:
        return dict(values={}, errors={}, cov=np.zeros((0, 0)), chi2=model.chi2([], data), valid=True, minos={})
    f = lambda *t: model.chi2(np.array(t), data)
    m = Minuit(f, *(start if start is not None else np.ones(n)), name=model.names)
    m.errordef = Minuit.LEAST_SQUARES
    for nm in model.names: m.limits[nm] = (0.0, 5.0)
    for nm in (fix or []): m.fixed[nm] = True
    m.migrad(); m.hesse()
    out = dict(values=dict(zip(model.names, m.values)), errors=dict(zip(model.names, m.errors)),
               cov=np.array(m.covariance) if m.covariance is not None else np.full((n, n), np.nan),
               chi2=float(m.fval), valid=bool(m.valid), minos={})
    if minos and m.valid:
        try:
            m.minos()
            out['minos'] = {nm: (m.merrors[nm].lower, m.merrors[nm].upper) for nm in model.names if not m.fixed[nm]}
        except Exception as e:  # noqa: BLE001
            out['minos_error'] = str(e)
    return out


def fit_iterated(model, data=None, n=5, tol=1e-3):
    """Fit with the covariance held fixed, re-evaluate the covariance at the result, repeat."""
    old = model.cov; model.cov = 'nominal'; r = fit(model, data)
    for i in range(n):
        th = np.array([r['values'][k] for k in model.names])
        model.cov = th; r2 = fit(model, data, start=th)
        th2 = np.array([r2['values'][k] for k in model.names]); r = r2
        if np.max(np.abs(th2 - th)) < tol: break
    r['iterations'] = i + 1; model.cov = old
    return r


def lr_test(model_null, model_alt, data=None):
    """Likelihood-ratio test between nested models; Wilks with the difference in parameter count."""
    f0, f1 = fit(model_null, data), fit(model_alt, data)
    d = f0['chi2'] - f1['chi2']; k = len(model_alt.names) - len(model_null.names)
    return dict(delta_chi2=d, ndf=k, p=float(stats.chi2.sf(max(d, 0), k)), null=f0, alt=f1)


def gof(model, res, data=None):
    data = model.data if data is None else data
    ndf = len(data) - len([n for n in model.names])
    return dict(chi2=res['chi2'], ndf=ndf, p=float(stats.chi2.sf(res['chi2'], ndf)))


def parameter_gof(models, joint, data_list=None, joint_data=None):
    """Parameter goodness-of-fit (Maltoni-Schwetz): chi2_joint_min - sum chi2_i_min against
    ndf = sum n_i - n_joint, for the SAME parameters fitted to independent subsets (e.g. FHC, RHC)."""
    fi = [fit(m, None if data_list is None else data_list[i]) for i, m in enumerate(models)]
    fj = fit(joint, joint_data)
    chi_pg = fj['chi2'] - sum(f['chi2'] for f in fi)
    ndf = sum(len(m.names) for m in models) - len(joint.names)
    return dict(chi2=chi_pg, ndf=ndf, p=float(stats.chi2.sf(max(chi_pg, 0), max(ndf, 1))), separate=fi, joint=fj)


def profile(model, name, grid, data=None):
    k = model.names.index(name); out = []
    for v in grid:
        start = np.ones(len(model.names)); start[k] = v
        r = fit(model, data, start=start, fix=[name]); out.append(r['chi2'])
    return np.array(out)


def conditional(model_all, cr_index, y_index, data_cr):
    """Gaussian conditional constraint (DocDB 32672, as ConstrainedCalculator): posterior mean and
    covariance of the bins y_index given the observed control-region bins cr_index. model_all spans both."""
    th = np.ones(len(model_all.names))
    mu = model_all.predict(th); V = model_all.cov_syst(th)
    x, y = np.array(cr_index), np.array(y_index)
    Vxx = V[np.ix_(x, x)] + model_all.cov_stat(mu[x], data_cr)
    K = V[np.ix_(y, x)] @ np.linalg.inv(Vxx)
    mean = mu[y] + K @ (data_cr - mu[x])
    cov = V[np.ix_(y, y)] - K @ V[np.ix_(x, y)]
    return mean, cov, K


def mcmc(model, data=None, n=20000, burn=4000, seed=1, step=None, chains=4):
    """Adaptive Metropolis (Haario et al.) on exp(-chi2/2) with flat priors on [0,5]. Returns samples
    [chains, n, npar], Gelman-Rubin R-hat and an effective-sample-size estimate."""
    rng = np.random.default_rng(seed); k = len(model.names)
    f0 = fit(model, data)
    x0 = np.array([f0['values'][nm] for nm in model.names])
    C = f0['cov'] if np.all(np.isfinite(f0['cov'])) else np.eye(k) * 0.01
    out = np.zeros((chains, n, k))
    for ch in range(chains):
        x = x0 + rng.multivariate_normal(np.zeros(k), C); x = np.clip(x, 1e-3, 5)
        lp = -0.5 * model.chi2(x, data); samples = []; Cp = C * 2.38 ** 2 / k
        for i in range(n + burn):
            if i > 500 and i % 500 == 0 and len(samples) > 200:
                Cp = np.cov(np.array(samples[-2000:]).T).reshape(k, k) * 2.38 ** 2 / k + 1e-8 * np.eye(k)
            y = rng.multivariate_normal(x, Cp)
            if np.all((y > 0) & (y < 5)):
                ly = -0.5 * model.chi2(y, data)
                if np.log(rng.random()) < ly - lp: x, lp = y, ly
            if i >= burn: out[ch, i - burn] = x
            samples.append(x.copy())
    m = out.mean(1); W = out.var(1, ddof=1).mean(0); Bv = n * m.var(0, ddof=1)
    rhat = np.sqrt(((n - 1) / n * W + Bv / n) / W)
    ess = []
    for j in range(k):
        s = out[:, :, j].ravel() - out[:, :, j].mean(); ac = np.correlate(s[:5000], s[:5000], 'full')[4999:] / (s[:5000] @ s[:5000])
        tau = 1 + 2 * np.sum(ac[1:np.argmax(ac < 0.05) if np.any(ac < 0.05) else 200])
        ess.append(out.shape[0] * n / tau)
    return dict(samples=out, rhat=dict(zip(model.names, rhat)), ess=dict(zip(model.names, ess)),
                mean=dict(zip(model.names, out.reshape(-1, k).mean(0))), std=dict(zip(model.names, out.reshape(-1, k).std(0))),
                q16=dict(zip(model.names, np.percentile(out.reshape(-1, k), 16, axis=0))),
                q84=dict(zip(model.names, np.percentile(out.reshape(-1, k), 84, axis=0))))


def toy(model, theta_true, rng, inject=None):
    """Pseudo-data: prediction at theta_true (optionally times per-bin factors 'inject'), a systematic
    throw from the covariance at theta_true, then Poisson."""
    mu = model.predict(theta_true)
    if inject is not None: mu = mu * inject
    V = model.cov_syst(theta_true)
    w, Q = np.linalg.eigh(V); w = np.clip(w, 0, None)
    shift = Q @ (np.sqrt(w) * rng.standard_normal(len(w)))
    return rng.poisson(np.clip(mu + shift, 0, None)).astype(float)
