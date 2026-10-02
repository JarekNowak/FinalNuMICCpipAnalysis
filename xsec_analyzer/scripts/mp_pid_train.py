"""mp_pid_train.py -- Phase 1 of report/MULTIPION_BNB_ADAPTATION_PLAN.md: a multiclass particle
classifier (muon / pion / proton / other) for the primary tracks of the multi-pion selections, trained
with XGBoost on the dump of macros/mp_pid/dump_pid_tracks.C and compared, on the same tracks, with the
pion identification the selections use today.

Training sample: backtracked purity > 0.5, track score >= 0.5, length >= 5 cm (the plan's definition),
classes balanced by weight. Runs are held out, not events: the model is trained on run periods 1, 2
and 4 (FHC and RHC) and tested on periods 3 (RHC) and 5 (FHC), so the test measures how the classifier
carries over to detector conditions it has not seen.

Comparison population: the pion-candidate pool of the selection without the muon candidate (as the
selection counts pions), track score >= 0.5. Pion-vs-rest, pion-vs-proton and pion-vs-muon ROC areas
for the new P(pi) and for the current scores: the split mp_pion_bdt (soft below 20 cm, hard above) as
CC1mu2pi applies it, the single mp_pion_bdt, the inclusive pion and MIP BDTs and the LLR score.

    python3 scripts/mp_pid_train.py [--max-train N]
Outputs in report/multipion/: phase1_classifier.json (all numbers), phase1_classifier.md (tables);
the model in /data/uboone/processed/mp_pid/model_xgb.json.
"""
import glob, json, os, sys, time
import numpy as np
import uproot
import xgboost as xgb
from sklearn.metrics import roc_auc_score, confusion_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
OUTDIR = os.path.join(REPO, 'report', 'multipion')
DUMP = '/data/uboone/processed/mp_pid/'
CLASSES = ['muon', 'pion', 'proton', 'other']
TRAIN_PERIODS, TEST_PERIODS = (1, 2, 4), (3, 5)
MISSING = -9999.
FEATURES = [
    'len', 'ts', 'dist', 'contained', 'start_contained', 'llr', 'llr_u', 'llr_v', 'llr_y',
    'bragg_p', 'bragg_mu', 'bragg_mip', 'bragg_pion', 'bragg_p_u', 'bragg_p_v', 'bragg_mu_u', 'bragg_mu_v',
    'bragg_pion_u', 'bragg_pion_v', 'bragg_mip_u', 'bragg_mip_v', 'fwd_p', 'fwd_mu', 'fwd_pion',
    'pida', 'pida_u', 'pida_v', 'chipr', 'chimu', 'chipi', 'chika', 'chipr_u', 'chipr_v', 'chimu_u', 'chimu_v',
    'chipi_u', 'chipi_v', 'trunk_dedx_y', 'trunk_dedx_u', 'trunk_dedx_v', 'trunk_rr_dedx_y', 'trunk_rr_dedx_u',
    'trunk_rr_dedx_v', 'calo_e_y', 'calo_e_u', 'calo_e_v', 'defl_mean', 'defl_stdev', 'defl_sep_mean',
    'mcs_mom', 'range_mom_mu', 'mcs_over_range', 'e_proton', 'end_sp', 'nhits_u', 'nhits_v', 'nhits_y',
    'planehits_u', 'planehits_v', 'planehits_y', 'n_trk_daughters', 'n_shr_daughters', 'n_descendents']
OLD = ['s_mppi_soft', 's_mppi_hard', 's_mppi', 's_pion', 's_mip']
# inputs whose ntuple branches are missing from the older productions (beam-on, beam-off and dirt of
# Runs 1-3, and the Run-4b/c/d and Run-5 beam-off and dirt for some of them; all overlays carry them):
# the model must not use them, or it cannot be applied to data consistently (--common)
NOT_EVERYWHERE = ['bragg_pion', 'bragg_pion_u', 'bragg_pion_v', 'fwd_p', 'fwd_mu', 'fwd_pion', 'trunk_dedx_y',
                  'trunk_dedx_u', 'trunk_dedx_v', 'trunk_rr_dedx_y', 'trunk_rr_dedx_u', 'trunk_rr_dedx_v', 'defl_mean',
                  'defl_stdev', 'defl_sep_mean', 'end_sp', 'nhits_u', 'nhits_v', 'nhits_y', 'n_descendents']
if '--common' in sys.argv:
    FEATURES = [f for f in FEATURES if f not in NOT_EVERYWHERE]
META = ['period', 'horn', 'label', 'charge', 'purity', 'true_p', 'is_mu_cand', 'w_cv', 'contained']


def load():
    cols = sorted(set(FEATURES + OLD + META))
    parts = []
    for f in sorted(glob.glob(DUMP + 'pidtrk-*.root')):
        a = uproot.open(f)['trk'].arrays(cols, library='np')
        parts.append(a)
        print(f'  {os.path.basename(f)}: {len(a["label"])} tracks', flush=True)
    d = {c: np.concatenate([p[c] for p in parts]) for c in cols}
    X = np.stack([d[c].astype(np.float32) for c in FEATURES], axis=1)
    # missing per-plane quantities (-9999 from the dumper, +-inf or float-max sentinels from the ntuple) are
    # set to the sentinel -9999, an ordinary value below every physical one, not to NaN: ROOT's RBDT, which
    # evaluates the model in the selection, has no missing-value branches, so XGBoost must not use them
    X[~np.isfinite(X) | (X <= -9998) | (np.abs(X) > 1e30)] = MISSING
    return d, X


def auc(y, s):
    return float(roc_auc_score(y, s)) if 0 < y.sum() < len(y) else float('nan')


def main():
    max_train = int(sys.argv[sys.argv.index('--max-train') + 1]) if '--max-train' in sys.argv else 0
    # --unmatched: tracks with backtracked purity <= 0.5 (mostly cosmic rays of the overlay, 23% of the
    # candidate pool) are kept and labelled "other" instead of being dropped, in training and in testing
    unmatched = '--unmatched' in sys.argv
    tag = ('_unmatched' if unmatched else '') + ('_common' if '--common' in sys.argv else '')
    t0 = time.time()
    d, X = load()
    lab = d['label'].copy()
    if unmatched:
        lab[~(d['purity'] > 0.5)] = 3
        base = (d['ts'] >= 0.5) & (d['len'] >= 5.)
    else:
        base = (d['purity'] > 0.5) & (d['ts'] >= 0.5) & (d['len'] >= 5.)
    tr = base & np.isin(d['period'], TRAIN_PERIODS)
    te = base & np.isin(d['period'], TEST_PERIODS)
    rng = np.random.default_rng(1)
    itr = np.flatnonzero(tr)
    if max_train and len(itr) > max_train: itr = rng.choice(itr, max_train, replace=False)
    # 10% of the training runs' tracks for early stopping
    rng.shuffle(itr); nval = len(itr) // 10; iva, itr = itr[:nval], itr[nval:]
    cnt = np.bincount(lab[itr], minlength=4).astype(float)
    cw = cnt.sum() / (4 * np.maximum(cnt, 1.))                  # balanced classes
    dtr = xgb.DMatrix(X[itr], label=lab[itr], weight=cw[lab[itr]], feature_names=FEATURES, missing=np.nan)
    dva = xgb.DMatrix(X[iva], label=lab[iva], weight=cw[lab[iva]], feature_names=FEATURES, missing=np.nan)
    params = dict(objective='multi:softprob', num_class=4, eval_metric='mlogloss', tree_method='hist',
                  max_depth=6, eta=0.1, subsample=0.8, colsample_bytree=0.8, min_child_weight=5, nthread=8, seed=1)
    print(f'training on {len(itr)} tracks (validation {len(iva)}), class counts {cnt.astype(int).tolist()}', flush=True)
    bst = xgb.train(params, dtr, num_boost_round=2000, evals=[(dva, 'val')], early_stopping_rounds=50, verbose_eval=100)
    bst = bst[: bst.best_iteration + 1]          # keep the trees up to the best iteration only
    bst.save_model(DUMP + f'model_xgb{tag}.json')  # scripts/mp_pid_export_rbdt.py writes the ROOT copy
    ite = np.flatnonzero(te)
    P = bst.predict(xgb.DMatrix(X[ite], feature_names=FEATURES, missing=np.nan))
    yt = lab[ite]
    out = dict(unmatched_as_other=unmatched, train_periods=TRAIN_PERIODS, test_periods=TEST_PERIODS, n_train=int(len(itr)), n_val=int(len(iva)),
               n_test=int(len(ite)), best_iteration=int(bst.num_boosted_rounds()) - 1, class_counts_train=cnt.astype(int).tolist(),
               test_class_counts=np.bincount(yt, minlength=4).tolist())
    cm = confusion_matrix(yt, P.argmax(1), labels=[0, 1, 2, 3])
    out['confusion_test'] = cm.tolist()
    out['auc_one_vs_rest_test'] = {CLASSES[k]: auc(yt == k, P[:, k]) for k in range(4)}
    gain = bst.get_score(importance_type='gain'); tot = sum(gain.values())
    out['importance_gain'] = sorted(((k, v / tot) for k, v in gain.items()), key=lambda x: -x[1])[:25]

    # ---- comparison on the selection's pion-candidate pool, held-out runs ----
    # CC1mu2pi counts a candidate as a pion if it is not the muon candidate, its vertex distance is within
    # pion_vtx_distance_cut() and the split BDT passes pion_bdt_cut_soft()/_hard() (base-class values unless
    # CC1mu2pi overrides them; the optional Bragg and MIP requirements are off)
    import re
    def header_value(name, default):
        val = default
        for h in ('CC1mu1piXp.hh', 'CC1mu2pi.hh'):          # the subclass override wins
            txt = open(os.path.join(REPO, 'xsec_analyzer', 'include', 'XSecAnalyzer', 'Selections', h)).read()
            m = re.search(name + r'\(\)\s*const(?:\s*override)?\s*\{\s*return\s*([-0-9.]+)', txt)
            if m: val = float(m.group(1))
        return val
    cut_soft, cut_hard = header_value('pion_bdt_cut_soft', 0.10), header_value('pion_bdt_cut_hard', -0.20)
    dcut = header_value('pion_vtx_distance_cut', 9.5)
    pool = (d['is_mu_cand'][ite] == 0) & (d['dist'][ite] <= dcut)
    newp = P[:, 1]
    soft = d['len'][ite] < 20.
    old = {'P(pi), new': newp,
           'mp_pion_bdt split (CC1mu2pi)': np.where(soft, d['s_mppi_soft'][ite], d['s_mppi_hard'][ite]),
           'mp_pion_bdt single': d['s_mppi'][ite], 'inclusive pion BDT': d['s_pion'][ite],
           'MIP BDT': d['s_mip'][ite], 'LLR score': d['llr'][ite]}
    hz, ch, cont = d['horn'][ite], d['charge'][ite], d['contained'][ite]
    groups = {'all': pool, 'FHC': pool & (hz == 0), 'RHC': pool & (hz == 1),
              'pi+ vs rest': pool & ((yt != 1) | (ch > 0)), 'pi- vs rest': pool & ((yt != 1) | (ch < 0)),
              'contained': pool & (cont == 1), 'uncontained': pool & (cont == 0),
              'length < 20 cm': pool & soft, 'length >= 20 cm': pool & ~soft}
    comp = {}
    for g, m in groups.items():
        row = {}
        for name, s in old.items():
            row[name] = dict(pi_vs_rest=auc(yt[m] == 1, s[m]),
                             pi_vs_p=auc(yt[m & np.isin(yt, (1, 2))] == 1, s[m & np.isin(yt, (1, 2))]),
                             pi_vs_mu=auc(yt[m & np.isin(yt, (0, 1))] == 1, s[m & np.isin(yt, (0, 1))]))
        comp[g] = dict(n=int(m.sum()), n_pion=int((yt[m] == 1).sum()), auc=row)
    out['comparison'] = comp

    # pion efficiency and purity at the CC1mu2pi working point, and the new P(pi) at the same efficiency
    pass_old = np.where(soft, d['s_mppi_soft'][ite] > cut_soft, d['s_mppi_hard'][ite] > cut_hard) & pool
    is_pi = (yt == 1) & pool
    eff_old = pass_old[is_pi].mean(); pur_old = (yt[pass_old] == 1).mean()
    thr = np.quantile(newp[is_pi], 1. - eff_old)
    pass_new = (newp > thr) & pool
    out['working_point'] = dict(cut_soft=cut_soft, cut_hard=cut_hard, dist_cut=dcut,
                                old=dict(pion_eff=float(eff_old), pion_purity=float(pur_old),
                                         proton_pass=float(pass_old[(yt == 2) & pool].mean()),
                                         muon_pass=float(pass_old[(yt == 0) & pool].mean())),
                                new_same_eff=dict(threshold=float(thr), pion_eff=float(pass_new[is_pi].mean()),
                                                  pion_purity=float((yt[pass_new] == 1).mean()),
                                                  proton_pass=float(pass_new[(yt == 2) & pool].mean()),
                                                  muon_pass=float(pass_new[(yt == 0) & pool].mean())))
    out['seconds'] = round(time.time() - t0)
    os.makedirs(OUTDIR, exist_ok=True)
    json.dump(out, open(os.path.join(OUTDIR, f'phase1_classifier{tag}.json'), 'w'), indent=1)
    write_md(out, tag)
    print(json.dumps({k: out[k] for k in ('n_train', 'n_test', 'best_iteration', 'auc_one_vs_rest_test', 'working_point')}, indent=1))


def write_md(o, tag=''):
    L = ['# Phase 1: multiclass particle classifier (held-out runs)', '',
         f'Trained on run periods {o["train_periods"]} ({o["n_train"]} tracks, balanced classes, {o["best_iteration"] + 1} trees), '
         f'tested on periods {o["test_periods"]} ({o["n_test"]} tracks). Track selection: backtracked purity > 0.5, '
         'track score >= 0.5, length >= 5 cm' + ('; tracks with purity <= 0.5 kept and labelled other.' if o.get('unmatched_as_other') else '.'), '',
         '## One-vs-rest ROC area, test runs', '', '| class | AUC |', '|---|---|']
    L += [f'| {k} | {v:.4f} |' for k, v in o['auc_one_vs_rest_test'].items()]
    L += ['', '## Confusion matrix, test runs (rows true, columns predicted: mu, pi, p, other)', '']
    for k, r in zip(CLASSES, o['confusion_test']):
        tot = max(sum(r), 1); L.append(f'- {k}: ' + ', '.join(f'{x} ({100 * x / tot:.1f}%)' for x in r))
    L += ['', '## Pion identification on the selection\'s candidate pool (no muon candidate), test runs', '']
    names = list(next(iter(o['comparison'].values()))['auc'].keys())
    for metric in ('pi_vs_rest', 'pi_vs_p', 'pi_vs_mu'):
        L += [f'### ROC area, {metric.replace("_", " ")}', '', '| group | n (pions) | ' + ' | '.join(names) + ' |',
              '|---|---|' + '---|' * len(names)]
        for g, c in o['comparison'].items():
            L.append(f'| {g} | {c["n"]} ({c["n_pion"]}) | ' + ' | '.join(f'{c["auc"][n][metric]:.3f}' for n in names) + ' |')
        L.append('')
    w = o['working_point']
    L += ['## Working point', '',
          f'The split mp_pion_bdt at its CC1mu2pi cuts (soft {w["cut_soft"]}, hard {w["cut_hard"]}): pion efficiency '
          f'{100 * w["old"]["pion_eff"]:.1f}%, pion purity {100 * w["old"]["pion_purity"]:.1f}%, protons passing '
          f'{100 * w["old"]["proton_pass"]:.1f}%, muons passing {100 * w["old"]["muon_pass"]:.1f}%.',
          f'The new P(pi) at the same pion efficiency: purity {100 * w["new_same_eff"]["pion_purity"]:.1f}%, protons passing '
          f'{100 * w["new_same_eff"]["proton_pass"]:.1f}%, muons passing {100 * w["new_same_eff"]["muon_pass"]:.1f}%.', '',
          '## Largest feature importances (gain share)', '']
    L += [f'- {k}: {v:.3f}' for k, v in o['importance_gain'][:15]]
    open(os.path.join(OUTDIR, f'phase1_classifier{tag}.md'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
