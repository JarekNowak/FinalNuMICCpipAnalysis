"""maximin_binning.py -- migration-diagonal binning scan on the current exposure, for the 0.68 and 0.50
criteria, with the edge conventions of the released bin configurations.

Step 1 (cache): read the processed overlays once per family (inclusive CC1mu1piXp, proton-tagged
CC1mu1pi1p) and horn mode, weight every event with the safe central-value weight times the per-run
exposure scale of Table tab:pot, and fill for each observable a fine true x reco matrix of selected
signal (reco under/overflow kept). Cached in logs/binreview/fine_<family>_<mode>.npz.
Step 2 (scan): for K = 2..KMAX, a maximin dynamic programme over the fine edges finds the binning
whose worst column-normalised diagonal, taken over FHC AND RHC, is largest, subject to a floor on the
selected signal per true bin (FHC, data exposure). Diagonal of a bin = selected signal reconstructed in
the same bin / selected signal reconstructed in any bin, with the first reconstructed bin open below
and the last open above for open-top observables (exactly as scripts/make_bin_configs.py and the
released configs), so it reproduces the unfolder's [DIAGDUMP] up to the fine-grid rounding.
Step 3 (evaluate): any explicit edge list (released, 0.50 candidates) is evaluated the same way.

    python3 report/tools/maximin_binning.py cache                 (reads the ntuples; ~10 min)
    python3 report/tools/maximin_binning.py scan [floor]          (default floor 50 selected signal, FHC)
    python3 report/tools/maximin_binning.py eval                  (released + candidate edges)
"""
import sys, os, json, numpy as np, uproot
OUT = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/binreview/'
P = '/data/uboone/processed/'
RUNS = {
 'fhc': [('Run1_fhc_new_numi_flux_fhc_pandora_ntuple', 0.09415), ('Run2_fhc_new_numi_flux_fhc_pandora_ntuple', 0.05085),
         ('Run4_fhc_new_numi_flux_fhc_pandora_ntuple', 0.07323), ('reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc', 0.11560)],
 'rhc': [('Run1_rhc_new_numi_flux_rhc_pandora_ntuple', 0.06728), ('Run2_rhc_new_numi_flux_rhc_pandora_ntuple', 0.04478)]
        + [(f'Run4{s}_rhc_new_numi_flux_rhc_pandora_ntuple', 0.08847) for s in 'abc']
        + [(f'Run3_rhc_new_numi_flux_rhc_pandora_ntuple_{s}', 0.09066) for s in ('aa', 'ab', 'ac', 'ad', 'ae')]}
PI = 3.1416
# name: (true expr, reco expr, lo, hi, step, open_top) ; X = selection prefix
OBS = {
 'pmu':     ('X_candidate_muon_mom_true', 'X_candidate_muon_mom_reco', 0.15, 3.0, 0.025, True),
 'ppi':     ('X_candidate_pion_mom_true', 'PPI', 0.175, 1.0, 0.005, True),
 'costhmu': ('X_candidate_muon_costh_true', 'X_candidate_muon_costh_reco', -1.0, 1.0, 0.025, False),
 'thetamu': ('ACOS:X_candidate_muon_costh_true', 'ACOS:X_candidate_muon_costh_reco', 0.0, 3.15, 0.02, False),
 'costhpi': ('X_candidate_pion_costh_true', 'X_candidate_pion_costh_reco', -1.0, 1.0, 0.025, False),
 'thmupi':  ('X_true_mu_pi_opening_angle', 'X_mu_pi_opening_angle', 0.0, 2.6, 0.025, False),
}
OBS1P = {
 'Whad':    ('X_W_had_true', 'X_W_had_reco', 0.0, 2.74, 0.02, True),
 'Wpipr':   ('X_W_pipr_true', 'X_W_pipr_reco', 1.08, 2.9, 0.01, True),
 'dpt':     ('X_deltaPt_true', 'X_deltaPt_reco', 0.0, 2.5, 0.025, True),
 'dalphat': ('X_deltaAlphaT_true', 'X_deltaAlphaT_reco', 0.0, 180.0, 2.5, False),
 'dphit':   ('X_deltaPhiT_true', 'X_deltaPhiT_reco', 0.0, 180.0, 2.5, False),
 'pn':      ('X_pn_true', 'X_pn_reco', 0.0, 2.0, 0.025, True),
 'thetap':  ('ACOS:X_proton_costh_true', 'ACOS:X_proton_costh_reco', 0.0, PI, 0.0314, False),
 'thpipr':  ('X_pi_pr_opening_angle_true', 'X_pi_pr_opening_angle_reco', 0.0, PI, 0.0314, False),
 # proton momentum (added 2026-09-26): the leading true proton, above the 0.3 GeV/c signal threshold; last bin open
 'pp':      ('X_proton_mom_true', 'X_proton_mom_reco', 0.3, 1.5, 0.02, True),
}
def fam_obs(fam): return dict(OBS, **OBS1P) if fam == '1p' else OBS
def grid(lo, hi, step): n = int(round((hi - lo) / step)); return np.linspace(lo, hi, n + 1)

def cache(fam, mode, only=None):
    """fine migration matrices of every observable of the family (or of ONLY, merged into the existing cache)"""
    X = 'CC1mu1pi1p' if fam == '1p' else 'CC1mu1piXp'
    obs = fam_obs(fam)
    if only: obs = {k: v for k, v in obs.items() if k in only}
    br = {X + '_MC_Signal', X + '_Selected', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight', X + '_candidate_pion_mom_reco'}
    for t, r, *_ in obs.values():
        for e in (t, r):
            if e != 'PPI': br.add(e.replace('ACOS:', '').replace('X_', X + '_'))
    mats = {k: None for k in obs}; gen = {k: None for k in obs}
    for name, sc in RUNS[mode]:
        f = P + ('w/' if fam == '1p' else '') + f'xsec-ana-{name}.root'
        for a in uproot.iterate(f + ':stv_tree', sorted(br), library='np', step_size='400 MB'):
            w = a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']
            w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.0) * sc
            sig = a[X + '_MC_Signal'].astype(bool); sel = a[X + '_Selected'].astype(bool)
            for k, (te, re_, lo, hi, st, ot) in obs.items():
                g = grid(lo, hi, st); n = len(g) - 1
                def val(e):
                    if e == 'PPI':
                        p = a[X + '_candidate_pion_mom_reco']; return np.sqrt((np.sqrt(p**2 + 0.011164) - 0.10566 + 0.13957)**2 - 0.019480)
                    v = a[e.replace('ACOS:', '').replace('X_', X + '_')].astype(float)
                    return np.arccos(np.clip(v, -1, 1)) if e.startswith('ACOS:') else v
                tv = val(te); rv = val(re_)
                ti = np.searchsorted(g, tv, side='right') - 1
                ti = np.where(ot & (tv >= hi), n - 1, ti)
                okt = sig & (ti >= 0) & (ti < n) & np.isfinite(tv)
                ri = np.searchsorted(g, rv, side='right') - 1           # -1 under, n over
                ri = np.clip(ri, -1, n); ri = np.where(np.isfinite(rv), ri, n) + 1   # 0 under ... n+1 over
                m = okt & sel
                M = np.zeros((n, n + 2)); np.add.at(M, (ti[m], ri[m]), w[m])
                G = np.zeros(n); np.add.at(G, ti[okt], w[okt])
                mats[k] = M if mats[k] is None else mats[k] + M
                gen[k] = G if gen[k] is None else gen[k] + G
        print('  read', name, flush=True)
    save = {f'M_{k}': v for k, v in mats.items()}; save.update({f'G_{k}': v for k, v in gen.items()})
    fn = OUT + f'fine_{fam}_{mode}.npz'
    if only and os.path.exists(fn):
        old = np.load(fn); save = {**{k: old[k] for k in old.files}, **save}
    np.savez_compressed(fn, **save)

def load(fam, mode):
    d = np.load(OUT + f'fine_{fam}_{mode}.npz'); obs = fam_obs(fam)
    return {k: (d[f'M_{k}'], d[f'G_{k}']) for k in obs}

def bin_stats(M, i, j, first, last, open_top):
    """diag, selected signal, reco-in-any-bin for merged fine true/reco range [i, j) (M has under at 0, over at n+1)"""
    n = M.shape[0]; rows = M[i:j]
    num = rows[:, 1 + i:1 + j].sum() + (rows[:, 0].sum() if first else 0) + (rows[:, n + 1].sum() if (last and open_top) else 0)
    den = rows[:, 1:1 + n].sum() + rows[:, 0].sum() + (rows[:, n + 1].sum() if open_top else 0)
    return (num / den if den > 0 else 0.0), den

def evaluate(fam, obsname, edges, mats):
    lo, hi, st, ot = fam_obs(fam)[obsname][2:]
    g = grid(lo, hi, st)
    idx = [int(np.argmin(np.abs(g - e))) for e in edges]
    out = {}
    for mode, M in mats.items():
        d, s = [], []
        for b in range(len(idx) - 1):
            dd, ss = bin_stats(M[obsname][0], idx[b], idx[b + 1], b == 0, b == len(idx) - 2, ot)
            d.append(dd); s.append(ss)
        out[mode] = (d, s)
    return out

def scan(fam, obsname, mats, floor, kmax=10):
    lo, hi, st, ot = fam_obs(fam)[obsname][2:]
    n = mats['fhc'][obsname][0].shape[0]
    D = {}
    def diag2(i, j, first, last):
        key = (i, j, first, last)
        if key not in D:
            dF, sF = bin_stats(mats['fhc'][obsname][0], i, j, first, last, ot)
            dR, _ = bin_stats(mats['rhc'][obsname][0], i, j, first, last, ot)
            D[key] = min(dF, dR) if sF >= floor else -1.0
        return D[key]
    res = {}
    for K in range(1, kmax + 1):
        # best[k][i]: best worst-diagonal for the first k bins ending at fine edge i
        best = np.full((K + 1, n + 1), -2.0); arg = np.zeros((K + 1, n + 1), dtype=int); best[0][0] = 9.0
        for k in range(1, K + 1):
            for i in range(1, n + 1):
                last = (k == K and i == n)
                if k == K and i != n: continue
                for j in range(k - 1, i):
                    if best[k - 1][j] < -1: continue
                    v = min(best[k - 1][j], diag2(j, i, j == 0, last))
                    if v > best[k][i]: best[k][i] = v; arg[k][i] = j
        if best[K][n] < 0: res[K] = None; continue
        e = [n]; i = n
        for k in range(K, 0, -1): i = arg[k][i]; e.append(i)
        g = grid(lo, hi, st); edges = [round(float(g[x]), 4) for x in e[::-1]]
        res[K] = (float(best[K][n]), edges)
    return res

if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'cache':
        fams = sys.argv[2:] or ['incl', '1p']
        for fam in fams:
            for mode in ('fhc', 'rhc'):
                print('caching', fam, mode, flush=True); cache(fam, mode)
    elif what == 'cache_only':          # cache_only <fam> <obs,obs>: add observables to an existing cache
        fam, only = sys.argv[2], sys.argv[3].split(',')
        for mode in ('fhc', 'rhc'):
            print('caching', fam, mode, only, flush=True); cache(fam, mode, only)
    elif what == 'scan':
        floor = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
        fams = sys.argv[3:] or ['incl', '1p']
        out = {}
        for fam in fams:
            if not (os.path.exists(OUT + f'fine_{fam}_fhc.npz') and os.path.exists(OUT + f'fine_{fam}_rhc.npz')):
                print('skip', fam, '(cache incomplete)'); continue
            mats = {m: load(fam, m) for m in ('fhc', 'rhc')}
            for o in fam_obs(fam):
                r = scan(fam, o, mats, floor, kmax=12); out[f'{fam}:{o}'] = r
                best68 = max([K for K, v in r.items() if v and v[0] >= 0.68] or [0])
                best50 = max([K for K, v in r.items() if v and v[0] >= 0.50] or [0])
                print(f'{fam:4s} {o:8s} finest K passing 0.68: {best68:2d} {r.get(best68, ("",""))[1] if best68 else ""}')
                print(f'{"":4s} {"":8s} finest K passing 0.50: {best50:2d} {r.get(best50, ("",""))[1] if best50 else ""}')
                print('      maximin worst diag by K: ' + ' '.join(f'{K}:{v[0]:.2f}' if v else f'{K}:--' for K, v in r.items()), flush=True)
        json.dump({k: {str(K): v for K, v in r.items()} for k, r in out.items()}, open(OUT + f'scan_floor{int(floor)}_{"_".join(fams)}.json', 'w'), indent=1)
