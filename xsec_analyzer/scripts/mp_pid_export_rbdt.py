"""mp_pid_export_rbdt.py -- write the Phase 1 particle classifier (scripts/mp_pid_train.py) in the ROOT
format of TMVA::Experimental::RBDT, which the selection can evaluate in C++, and check that RBDT reproduces
XGBoost on held-out tracks.
    python3 scripts/mp_pid_export_rbdt.py [_unmatched]  ->  /data/uboone/processed/mp_pid/mp_pid_rbdt<tag>.root (key "mp_pid")
"""
import os, sys
import numpy as np
import uproot
import xgboost as xgb
import ROOT
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_pid_train as T

tag = sys.argv[1] if len(sys.argv) > 1 else ''
bst = xgb.Booster(); bst.load_model(T.DUMP + f'model_xgb{tag}.json')
clf = xgb.XGBClassifier(); clf._Booster = bst; clf.n_classes_ = 4; clf.objective = 'multi:softprob'
out = T.DUMP + f'mp_pid_rbdt{tag}.root'
ROOT.TMVA.Experimental.SaveXGBoost(clf, 'mp_pid', out, num_inputs=len(T.FEATURES))
a = uproot.open(T.DUMP + 'pidtrk-reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root')['trk'].arrays(T.FEATURES, library='np', entry_stop=20000)
X = np.stack([a[c].astype(np.float32) for c in T.FEATURES], axis=1)
X[~np.isfinite(X) | (X <= -9998) | (np.abs(X) > 1e30)] = T.MISSING
p_x = bst.predict(xgb.DMatrix(X, feature_names=T.FEATURES))
p_r = np.asarray(ROOT.TMVA.Experimental.RBDT('mp_pid', out).Compute(X)).reshape(len(X), -1)
print(f'wrote {out}; RBDT vs XGBoost on {len(X)} Run-5 tracks: max |difference| {np.abs(p_r - p_x).max():.2e}', flush=True)
os._exit(0)      # ROOT and XGBoost both register exit handlers; skip them
