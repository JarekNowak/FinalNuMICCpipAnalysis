"""cr_llr_detvar.py -- does a detector variation move the LLR of long proton-like tracks as the beam-on data do
(scripts/cr_llr_spike.py)? Same track population in the beam-on data, the overlays and the Run-4 detector-variation
samples: tracks of the two-pion classifier's pool (generation 2, track score >= 0.5, valid LLR score) at least 40 cm
long with both ends in the containment volume, in events of the four opened control regions. Per sample: the median
LLR sum and collection-plane LLR per hit of the proton-like ones (score < -0.5), and the share of them beyond -0.9.
Then the event score of the two-pion domain: the shape of each variation over the CV against data over prediction.

Data and overlay from the control-region skim (pool from the track dump; the emulation below is checked against it on
the overlay); detector variations from /data/uboone/processed/mp2pi (CV and the eight variations, FHC and RHC Run 4).

    python3 scripts/cr_llr_detvar.py  ->  report/planning/multipion/cr_llr_detvar.{json,md}
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'report', 'planning', 'multipion'))
SKIM = '/data/uboone/processed/mp2pi_cr_skim/'
DV = '/data/uboone/processed/mp2pi/'
KNOBS = ['CV', 'LYdown', 'LYrayl', 'Recomb2', 'SCE', 'WMAngleXZ', 'WMAngleYZ', 'WMX', 'WMYZ']
CR = ['CC1mu1piXp_sb_cc0pi', 'CC1mu1piXp_sb_multipi', 'CC1mu1piXp_sb_pi0', 'CC1mu1piXp_sb_cosmic']
CFV = (2.0, 254.35, -113.53, 107.47, 2.1, 1034.9)     # CC1mu1piXp containment volume
TRK = ['pfp_generation_v', 'trk_score_v', 'trk_llr_pid_score_v', 'trk_len_v', 'trk_llr_pid_v', 'trk_llr_pid_y_v',
       'pfnplanehits_Y', 'trk_sce_start_x_v', 'trk_sce_start_y_v', 'trk_sce_start_z_v', 'trk_sce_end_x_v',
       'trk_sce_end_y_v', 'trk_sce_end_z_v']


def inside(x, y, z):
    return (x > CFV[0]) & (x < CFV[1]) & (y > CFV[2]) & (y < CFV[3]) & (z > CFV[4]) & (z < CFV[5])


def emulate(path):
    """pool tracks (long, both ends contained) of control-region events, emulated from the ntuple vectors"""
    t = uproot.open(path)['stv_tree']
    out = {'score': [], 'sum': [], 'perhit': []}
    for a in t.iterate(CR + ['nslice'] + TRK, step_size='300 MB', library='ak'):
        ev = np.zeros(len(a), bool)
        for c in CR: ev |= ak.to_numpy(a[c])
        a = a[ev & (ak.to_numpy(a['nslice']) == 1)]
        if len(a) == 0: continue
        s = a['trk_llr_pid_score_v']
        # the pool's Bragg-likelihood validity requirement is left out (the processed files do not carry the Bragg
        # branches); the overlay row "emulated" against "dump" shows its size
        m = (a['pfp_generation_v'] == 2) & (a['trk_score_v'] >= 0.5) & (s > -1) & (s < 2) & (a['trk_len_v'] >= 40) \
            & inside(a['trk_sce_start_x_v'], a['trk_sce_start_y_v'], a['trk_sce_start_z_v']) \
            & inside(a['trk_sce_end_x_v'], a['trk_sce_end_y_v'], a['trk_sce_end_z_v'])
        out['score'].append(ak.to_numpy(ak.flatten(s[m])))
        out['sum'].append(ak.to_numpy(ak.flatten(a['trk_llr_pid_v'][m])))
        out['perhit'].append(ak.to_numpy(ak.flatten(a['trk_llr_pid_y_v'][m] / np.maximum(a['pfnplanehits_Y'][m], 1))))
    return {k: np.concatenate(v).astype(float) if v else np.zeros(0) for k, v in out.items()}


def dumped(path):
    """the same population from the track dump of the skim"""
    t = uproot.open(path)['stv_tree']
    B = 'CC1mu2piBDTEvt'
    out = {'score': [], 'sum': [], 'perhit': []}
    for a in t.iterate([f'{B}_etrk_{v}' for v in ('idx', 'len', 'llr', 'contained')] + ['trk_llr_pid_v', 'trk_llr_pid_y_v',
                       'pfnplanehits_Y'], step_size='300 MB', library='ak'):
        idx = a[f'{B}_etrk_idx']
        m = (a[f'{B}_etrk_len'] >= 40) & (a[f'{B}_etrk_contained'] == 3)
        out['score'].append(ak.to_numpy(ak.flatten(a[f'{B}_etrk_llr'][m])))
        out['sum'].append(ak.to_numpy(ak.flatten(a['trk_llr_pid_v'][idx][m])))
        out['perhit'].append(ak.to_numpy(ak.flatten(a['trk_llr_pid_y_v'][idx][m] / np.maximum(a['pfnplanehits_Y'][idx][m], 1))))
    return {k: np.concatenate(v).astype(float) if v else np.zeros(0) for k, v in out.items()}


SCORE = np.linspace(0, 0.92, 24)


def score_hist(path):
    """event-score distribution of the two-pion domain (regions 1-4) in control-region events, CV-weighted"""
    t = uproot.open(path)['stv_tree']; h = np.zeros(len(SCORE) - 1)
    br = CR + ['CC1mu2piBDT_evt_region', 'CC1mu2piBDT_evt_score', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
    for a in t.iterate(br, step_size='300 MB', library='np'):
        m = np.zeros(len(a[CR[0]]), bool)
        for c in CR: m |= a[c].astype(bool)
        m &= a['CC1mu2piBDT_evt_region'] >= 1
        w = a['tuned_cv_weight'] * a['ppfx_cv_weight'] * a['normalisation_weight']
        w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)
        h += np.histogram(np.clip(a['CC1mu2piBDT_evt_score'][m], 0, SCORE[-1] - 1e-9), SCORE, weights=w[m])[0]
    return h


def metrics(x):
    p = x['score'] < -0.5
    n = int(p.sum())
    return dict(n_long=int(len(x['score'])), n_protonlike=n, protonlike_share=n / max(len(x['score']), 1),
                median_sum=float(np.median(x['sum'][p])) if n else float('nan'),
                median_perhit_y=float(np.median(x['perhit'][p])) if n else float('nan'),
                beyond_m09=float((x['score'][p] < -0.9).mean()) if n else float('nan'))


def main():
    res = {}
    for lab, f in (('beam-on FHC Run 4', 'beamon_fhc_run4c'), ('beam-on FHC Run 4d', 'beamon_fhc_run4d'),
                   ('beam-on RHC Run 4a', 'beamon_rhc_run4a'), ('beam-on RHC Run 4b', 'beamon_rhc_run4b'),
                   ('overlay FHC Run 4', 'Run4_fhc_new_numi_flux_fhc_pandora_ntuple'),
                   ('overlay RHC Run 4c', 'Run4c_rhc_new_numi_flux_rhc_pandora_ntuple')):
        p = SKIM + f'xsec-ana-{f}.root'
        if os.path.exists(p): res[lab + ' (dump)'] = metrics(dumped(p))
    res['overlay FHC Run 4 (emulated)'] = metrics(emulate(SKIM + 'xsec-ana-Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root'))
    for mode in ('fhc', 'rhc'):
        for k in KNOBS:
            p = DV + f'xsec-ana-detvar_run4{mode}_{k}.root'
            res[f'detector variation {mode.upper()} {k}'] = metrics(emulate(p))
            print(mode, k, res[f'detector variation {mode.upper()} {k}'], flush=True)
    json.dump(res, open(os.path.join(OUTDIR, 'cr_llr_detvar.json'), 'w'), indent=1)
    L = ['# LLR of long proton-like tracks: beam-on data, overlays and detector variations', '',
         'Generated by `scripts/cr_llr_detvar.py`. Pool tracks at least 40 cm long with both ends contained, in events of '
         'the four opened control regions; proton-like means LLR score < -0.5 (LLR sum < -100). Unweighted.', '',
         '| sample | long tracks | proton-like share | median LLR sum | median collection-plane LLR per hit | share beyond -0.9 |',
         '|---|---|---|---|---|---|']
    for k, v in res.items():
        L.append(f"| {k} | {v['n_long']} | {v['protonlike_share']:.3f} | {v['median_sum']:.0f} | {v['median_perhit_y']:.3f} | "
                 f"{v['beyond_m09']:.3f} |")
    # event score: shape of each variation over the CV against the shape of data over prediction (cr_mp2pi_check.json)
    chk = json.load(open(os.path.join(OUTDIR, 'cr_mp2pi_check.json')))['all']['event | score']
    dat, pred = np.array(chk['data_hist']), np.array(chk['pred_hist'])
    rd = (dat / dat.sum()) / (pred / pred.sum())
    sh = {}
    for mode in ('fhc', 'rhc'):
        hcv = score_hist(DV + f'xsec-ana-detvar_run4{mode}_CV.root')
        for k in KNOBS[1:]:
            hv = score_hist(DV + f'xsec-ana-detvar_run4{mode}_{k}.root')
            sh[f'{mode.upper()} {k}'] = (hv / hv.sum()) / (hcv / hcv.sum())
    res['score_shape'] = {'data_over_prediction': list(rd), **{k: list(v) for k, v in sh.items()}}
    L += ['', '## Event score: shape of data over prediction and of each variation over the CV', '',
          'Two-pion domain (regions 1-4) of control-region events, each distribution normalised to one; the correlation '
          'is that of the residuals (ratio - 1) over the 23 bins.', '',
          '| variation | bin 0-0.04 | bin 0.88-0.92 | correlation with data |', '|---|---|---|---|',
          f'| data / prediction | {rd[0]:.3f} | {rd[-1]:.3f} | |']
    for k, v in sh.items():
        L.append(f'| {k} | {v[0]:.3f} | {v[-1]:.3f} | {np.corrcoef(rd - 1, v - 1)[0, 1]:.2f} |')
    json.dump(res, open(os.path.join(OUTDIR, 'cr_llr_detvar.json'), 'w'), indent=1)
    open(os.path.join(OUTDIR, 'cr_llr_detvar.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
