"""mp_evt_export_rbdt.py -- write the two-pion event classifiers of scripts/mp_evt.py (XGBoost json) in the ROOT format of
TMVA::Experimental::RBDT (key "evt_clf", read by CC1mu1piXp::classify_event), and check RBDT against XGBoost on the
rows of the study dump. Two-fold tags (ending in 2f) export both fold models.

    MP_EVT_DIR=/data/uboone/processed/mp_evt_llronly python3 scripts/mp_evt_export_rbdt.py llronly2f pid
        ->  <MP_EVT_DIR>/evt_clf_n2_<tag>[_<fold>]_<set>_rbdt.root
"""
import os, sys
import numpy as np
import xgboost as xgb
import ROOT
sys.argv = [sys.argv[0]] + sys.argv[1:]
tag, fs = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_evt as M

D = np.load(M.NEW + 'events.npz')
names = M.feat_names(2, fs)
X = np.stack([D[f'n2_{k}'] for k in names], axis=1).astype(np.float32)
dom = D['n2_d6'].astype(bool)
Xc = X[dom][:20000]
for fold in (['even', 'odd'] if tag.endswith('2f') else ['']):
    stem = f'evt_clf_n2_{tag}' + (f'_{fold}' if fold else '') + f'_{fs}'
    bst = xgb.Booster(); bst.load_model(M.NEW + stem + '.json')
    clf = xgb.XGBClassifier(); clf._Booster = bst; clf.n_classes_ = 2; clf.objective = 'binary:logistic'
    out = M.NEW + stem + '_rbdt.root'
    ROOT.TMVA.Experimental.SaveXGBoost(clf, 'evt_clf', out, num_inputs=len(names))
    p_x = bst.predict(xgb.DMatrix(Xc, feature_names=names))
    p_r = np.asarray(ROOT.TMVA.Experimental.RBDT('evt_clf', out).Compute(Xc)).reshape(len(Xc), -1)[:, 0]
    print(f'wrote {out}: {bst.num_boosted_rounds()} trees; RBDT vs XGBoost on {len(Xc)} rows: '
          f'max |difference| {np.abs(p_r - p_x).max():.2e}', flush=True)
os._exit(0)      # ROOT and XGBoost both register exit handlers; skip them
