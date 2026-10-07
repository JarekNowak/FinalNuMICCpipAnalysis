#!/usr/bin/env python3
"""mp_pid_detsens.py -- detector sensitivity of the multi-pion particle classifier (follow-up of Phase 2 of
report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md: the detector term of the two-pion working point is led by Recomb2
and the wire-modification variations).

Input: /data/uboone/processed/mp_pid_detvar/pidtrk-detvar_run4{fhc,rhc}_<variation>[_partN].root, the Run-4
detector-variation samples dumped with macros/mp_pid/dump_pid_tracks.C (slurm/slurm_mp_pid_dump.sbatch,
slurm/mp_pid_detvar_manifest.list), the same track pool as the training sample.

Paired comparison: a track of the CV sample and a track of a variation are the same particle when they belong to the
same event (run, subrun, event) and carry the same backtracked PDG code and true momentum; particles that appear more
than once in an event in either sample are dropped. Tracks as in the training (backtracked purity > 0.5, track score
>= 0.5, length >= 5 cm). For every input: the median paired shift over the CV spread (half the 16-84% range), per true
class. For a classifier: the change of the fraction of true pions, protons and muons above the P(pi) threshold that
keeps 80% of the CV pions (paired, so the statistical noise of the samples cancels to first order), on the events with
an odd event number only (the even ones may enter an augmented training, scripts/mp_pid_train.py --aug).

    python3 scripts/mp_pid_detsens.py inputs            -> report/planning/multipion/detsens_inputs.{json,md}
    python3 scripts/mp_pid_detsens.py model <model.json> [--features common|<list file>]  (prints the metric)
"""
import glob, json, os, sys
import numpy as np
import uproot
import xgboost as xgb

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
OUTDIR = os.path.join(REPO, 'report', 'planning', 'multipion')
DV = '/data/uboone/processed/mp_pid_detvar/'
VARS = ['LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']
sys.argv_saved = list(sys.argv)
sys.argv = [sys.argv[0], '--common']           # the deployed model's input list
import mp_pid_train as T                        # noqa: E402  (FEATURES, MISSING, sentinel handling)
sys.argv = sys.argv_saved
COMMON = list(T.FEATURES)
META = ['run', 'sub', 'evt', 'pdg', 'true_p', 'label', 'purity', 'ts', 'len']
CLS = {0: 'muon', 1: 'pion', 2: 'proton'}
HORNS = sys.argv[sys.argv.index('--horns') + 1].split(',') if '--horns' in sys.argv else ['fhc', 'rhc']


def load(horn, var, cols):
    parts = sorted(glob.glob(f'{DV}pidtrk-detvar_run4{horn}_{var}*.root'))
    if not parts: raise SystemExit(f'no dump for {horn} {var}')
    arrs = [uproot.open(p)['trk'].arrays(cols, library='np') for p in parts]
    d = {c: np.concatenate([a[c] for a in arrs]) for c in cols}
    keep = (d['purity'] > 0.5) & (d['ts'] >= 0.5) & (d['len'] >= 5.)
    return {c: v[keep] for c, v in d.items()}


def pair(a, b):
    """Indices (ia, ib) of the tracks of a and b that are the same particle (unique in both)."""
    import pandas as pd
    K = ['run', 'sub', 'evt', 'pdg', 'p']
    def frame(d, col):
        return pd.DataFrame({'run': d['run'].astype(np.int64), 'sub': d['sub'].astype(np.int64),
                             'evt': d['evt'].astype(np.int64), 'pdg': d['pdg'].astype(np.int64),
                             'p': np.round(d['true_p'].astype(np.float64), 5), col: np.arange(len(d['run']))})
    da, db = frame(a, 'ia'), frame(b, 'ib')
    da = da[~da.duplicated(K, keep=False)]; db = db[~db.duplicated(K, keep=False)]
    m = da.merge(db, on=K)
    return m['ia'].to_numpy(), m['ib'].to_numpy()


def matrix(d, feats):
    X = np.stack([d[c].astype(np.float32) for c in feats], axis=1)
    X[~np.isfinite(X) | (X <= -9998) | (np.abs(X) > 1e30)] = T.MISSING
    return X


def inputs():
    cols = sorted(set(COMMON + META))
    res = {}
    for horn in HORNS:
        cv = load(horn, 'CV', cols)
        for var in VARS:
            v = load(horn, var, cols)
            ia, ib = pair(cv, v)
            r = {'paired_tracks': int(len(ia)), 'cv_tracks': int(len(cv['run'])), 'shift': {}}
            for f in COMMON:
                x0, x1 = cv[f][ia].astype(float), v[f][ib].astype(float)
                ok = np.isfinite(x0) & np.isfinite(x1) & (x0 > -9998) & (x1 > -9998) & (np.abs(x0) < 1e30) & (np.abs(x1) < 1e30)
                r['shift'][f] = {}
                for lab, name in CLS.items():
                    m = ok & (cv['label'][ia] == lab)
                    if m.sum() < 50: continue
                    q16, q84 = np.percentile(x0[m], [16, 84])
                    sp = max((q84 - q16) / 2., 1e-9)
                    r['shift'][f][name] = float(np.median(x1[m] - x0[m]) / sp)
            res[f'{horn}_{var}'] = r
            print(f'  {horn} {var}: {len(ia)} paired of {len(cv["run"])}', flush=True)
    json.dump(res, open(os.path.join(OUTDIR, 'detsens_inputs.json'), 'w'), indent=1)
    # table: inputs ranked by the largest |shift| over the variations, per class
    L = ['# Detector sensitivity of the particle-classifier inputs', '',
         'Generated by `scripts/mp_pid_detsens.py inputs`. Median paired shift (variation minus CV, same particle) in units of',
         'half the 16-84% range of the CV distribution, per true class; FHC and RHC Run-4 samples. Largest over the variations.', '',
         '| Input | Pions | Protons | Muons | Variation leading the pion shift |', '|---|---|---|---|---|']
    rows = []
    for f in COMMON:
        best = {}
        for name in CLS.values():
            vals = [(abs(r['shift'][f].get(name, 0.)), r['shift'][f].get(name, 0.), k) for k, r in res.items()]
            best[name] = max(vals)
        rows.append((best['pion'][0], f, best))
    for _, f, best in sorted(rows, key=lambda x: -x[0]):
        L.append(f'| {f} | {best["pion"][1]:+.3f} | {best["proton"][1]:+.3f} | {best["muon"][1]:+.3f} | {best["pion"][2]} |')
    open(os.path.join(OUTDIR, 'detsens_inputs.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:30]))


def model_metric(path, feats):
    """Change of the fraction above the 80%-pion threshold, per variation and class (paired)."""
    bst = xgb.Booster(); bst.load_model(path)
    cols = sorted(set(feats + META))
    out = {}
    for horn in HORNS:
        cv = load(horn, 'CV', cols)
        Xc = matrix(cv, feats)
        Pc = bst.predict(xgb.DMatrix(Xc, feature_names=feats, missing=np.nan))
        thr = np.percentile(Pc[(cv['label'] == 1) & (cv['evt'].astype(np.int64) % 2 == 1), 1], 20.)   # keeps 80% of the CV pions
        for var in VARS:
            v = load(horn, var, cols)
            ia, ib = pair(cv, v)
            odd = cv['evt'][ia].astype(np.int64) % 2 == 1          # events never used in an augmented training
            ia, ib = ia[odd], ib[odd]
            Pv = bst.predict(xgb.DMatrix(matrix(v, feats)[ib], feature_names=feats, missing=np.nan))
            r = {}
            for lab, name in CLS.items():
                m = cv['label'][ia] == lab
                f0 = (Pc[ia][m, 1] > thr).mean(); f1 = (Pv[m, 1] > thr).mean()
                r[name] = dict(cv=float(f0), var=float(f1), rel=float(f1 / f0 - 1.) if f0 > 0 else float('nan'))
            out[f'{horn}_{var}'] = r
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'inputs':
        inputs()
    elif sys.argv[1] == 'model':
        feats = COMMON
        if '--features' in sys.argv:
            fl = sys.argv[sys.argv.index('--features') + 1]
            if fl != 'common': feats = [x.strip() for x in open(fl) if x.strip()]
        res = model_metric(sys.argv[2], feats)
        od = os.path.join(OUTDIR, 'detrobust'); os.makedirs(od, exist_ok=True)
        json.dump(res, open(os.path.join(od, 'detsens_' + os.path.basename(sys.argv[2])), 'w'), indent=1)
        for k, r in res.items():
            print(f'{k:22s} ' + '  '.join(f'{n}: {100 * r[n]["rel"]:+5.1f}% (cv {100 * r[n]["cv"]:.1f}%)' for n in CLS.values()))
        tot = {n: float(np.sqrt(np.mean([[r[n]['rel'] ** 2 for k2, r in res.items() if k2.startswith(h)] for h in HORNS]))) for n in CLS.values()}
        print('rms over variations and horns:', {n: f'{100 * v:.1f}%' for n, v in tot.items()})
