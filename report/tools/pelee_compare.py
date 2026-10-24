"""pelee_compare.py -- compare the PeLEE ntuples used by the analysis with the custom PeLEE ntuples in
/data/uboone/temp/custom_pelee_ntuples (and reweighted_custom_pelee_ntuples), pair by pair (same sample type and run
period). Metadata only; no selection is applied, so the blind signal region is not touched.

Per file: entries of nuselection/NeutrinoSelectionFilter, subruns and summed POT of nuselection/SubRun, branch list,
the length of the weight vectors. Per pair: branches present in one file only, the overlap of the (run, subrun) and
(run, subrun, event) sets, and for the common events the agreement of reconstruction quantities (slice count,
topological score, vertex, track and shower counts) and of the weights (ppfx_cv, weightSplineTimesTune, and the first
entries of weightsPPFX, weightsGenie).

    python3 report/tools/pelee_compare.py [N]                ->  report/planning/pelee_compare.json (summary printed)
    python3 report/tools/pelee_compare.py [N] --reweighted   ->  our overlays against reweighted_custom_pelee_ntuples,
                                                                   report/planning/pelee_compare_reweighted.json
"""
import json, os, sys
from multiprocessing import Pool
import numpy as np
import uproot, awkward as ak

C = '/data/uboone/temp/custom_pelee_ntuples/'
R = '/data/uboone/temp/reweighted_custom_pelee_ntuples/'
OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'planning', 'pelee_compare.json'))
O = '/data/uboone/new_numi_flux/'
PAIRS = [
    ('overlay FHC Run 1', [O + 'Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root'], [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run1_fhc.root']),
    ('overlay RHC Run 1', [O + 'Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root'], [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run1_rhc.root']),
    ('overlay FHC Run 2', [O + 'Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root'], [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run2_fhc.root']),
    ('overlay RHC Run 2', [O + 'Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root'], [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run2_rhc.root']),
    ('overlay RHC Run 3', [O + f'Run3_rhc_new_numi_flux_rhc_pandora_ntuple_a{c}.root' for c in 'abcde'],
     [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run3_rhc.root', C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run3b_rhc_extra.root']),
    ('overlay FHC Run 4', [O + 'Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root'],
     [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run4_fhc.root', C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run4_fhc_extra.root']),
    ('overlay RHC Run 4', [O + f'Run4{c}_rhc_new_numi_flux_rhc_pandora_ntuple.root' for c in 'abc'],
     [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run4_rhc.root', C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run4_rhc_extra.root']),
    ('overlay FHC Run 5', [O + 'reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root'],
     [C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run5_fhc.root', C + 'nu_overlay_withppfxfixes/numi_nu_overlay_pion_ntuples_run5_fhc_extra.root']),
    ('beam-off Run 1', ['/data/uboone/EXT/neutrinoselection_filt_run1_beamoff.root'],
     [C + 'beam_off/numi_beam_off_pion_ntuples_run1_fhc.root', C + 'beam_off/numi_beam_off_pion_ntuples_run1_rhc.root']),
    ('beam-off Run 3b', ['/data/uboone/EXT/neutrinoselection_filt_run3b_beamoff.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run3b_rhc.root']),
    ('beam-off Run 4a', ['/data/uboone/EXT/numi_pelee_ntuple_beam_off_run4a_rhc_ana.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run4a_rhc.root']),
    ('beam-off Run 4b', ['/data/uboone/EXT/numi_pelee_ntuple_beam_off_run4b_rhc_ana.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run4b_rhc.root']),
    ('beam-off Run 4c', ['/data/uboone/EXT/numi_pelee_ntuple_beam_off_run4c_fhc_ana.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run4c_fhc.root']),
    ('beam-off Run 4d', ['/data/uboone/EXT/numi_pelee_ntuple_beam_off_run4d_fhc_ana.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run4d_fhc.root']),
    ('beam-off Run 5', ['/data/uboone/EXT/numi_pelee_ntuple_beam_off_run5_fhc_ana.root'], [C + 'beam_off/numi_beam_off_pion_ntuples_run5_fhc.root']),
    ('beam-on FHC Run 1', ['/data/uboone/beam_on/neutrinoselection_filt_run1_beamon_beamgood.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run1_fhc.root']),
    ('beam-on RHC Run 1', ['/data/uboone/beam_on/rhc/neutrinoselection_filt_run1_RHC_beamon_beamgood.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run1_rhc.root']),
    ('beam-on Run 3b', ['/data/uboone/beam_on/neutrinoselection_filt_run3b_beamon_beamgood.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run3b_rhc.root']),
    ('beam-on Run 4a', ['/data/uboone/beam_on/numi_pelee_ntuple_beam_on_run4a_ana.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run4a_rhc.root']),
    ('beam-on Run 4b', ['/data/uboone/beam_on/numi_pelee_ntuple_beam_on_run4b_ana.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run4b_rhc.root']),
    ('beam-on Run 4c', ['/data/uboone/beam_on/numi_pelee_ntuple_beam_on_run4c_ana.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run4c.root']),
    ('beam-on Run 4d', ['/data/uboone/beam_on/numi_pelee_ntuple_beam_on_run4d_fhc_ana.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run4d_fhc_goodruns.root']),
    ('beam-on Run 5', ['/data/uboone/beam_on/numi_pelee_ntuple_beam_on_run5_ana.root'], [C + 'beam_on/numi_beam_on_pion_ntuples_run5_fhc.root']),
    ('dirt FHC Run 1', ['/data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run1_fhc_slim.root']),
    ('dirt Run 3b', ['/data/uboone/dirt/neutrinoselection_filt_run3b_dirt_overlay.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run3b_rhc_slim.root']),
    ('dirt Run 4a', ['/data/uboone/dirt/run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run4a_rhc_slim.root']),
    ('dirt Run 4b', ['/data/uboone/dirt/numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run4b_rhc_slim.root']),
    ('dirt Run 4c', ['/data/uboone/dirt/numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run4c_fhc_slim.root']),
    ('dirt Run 4d', ['/data/uboone/dirt/numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run4d_fhc_slim.root']),
    ('dirt Run 5', ['/data/uboone/dirt/run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns.root'], [C + 'dirt_overlay/numi_dirt_overlay_pion_ntuples_run5_fhc_slim.root']),
]
RW = R + 'nu_overlay/reweightedPPFX_numi_nu_overlay_pion_ntuples_'
PAIRS_RW = [
    ('overlay FHC Run 1 (rw)', [O + 'Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root'], [RW + 'run1_fhc.root']),
    ('overlay RHC Run 1 (rw)', [O + 'Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root'], [RW + 'run1_rhc.root']),
    ('overlay FHC Run 2 (rw)', [O + 'Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root'], [RW + 'run2_fhc.root']),
    ('overlay RHC Run 2 (rw)', [O + 'Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root'], [RW + 'run2_rhc.root']),
    ('overlay RHC Run 3 (rw)', [O + f'Run3_rhc_new_numi_flux_rhc_pandora_ntuple_a{c}.root' for c in 'abcde'], [RW + 'run3_rhc.root', RW + 'run3b_rhc_extra.root']),
    ('overlay FHC Run 4 (rw)', [O + 'Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root'], [RW + 'run4_fhc.root', RW + 'run4_fhc_extra.root']),
    ('overlay RHC Run 4 (rw)', [O + f'Run4{c}_rhc_new_numi_flux_rhc_pandora_ntuple.root' for c in 'abc'], [RW + 'run4_rhc.root', RW + 'run4_rhc_extra.root']),
    ('overlay FHC Run 5 (rw)', [O + 'reweightedPPFX_numi_nu_overlay_pion_ntuples_run5_fhc.root'], [RW + 'run5_fhc.root', RW + 'run5_fhc_extra.root']),
]
EV = ['run', 'sub', 'evt']
RECO = ['nslice', 'topological_score', 'reco_nu_vtx_x', 'reco_nu_vtx_z', 'n_tracks', 'n_showers', 'CosmicIP']
WGT = ['ppfx_cv', 'weightSplineTimesTune', 'nu_e', 'nu_pdg']


def info(paths):
    out = dict(entries=0, subruns=0, pot=0., branches=None, wlen={}, files=[])
    ids, rs, cols = [], set(), {}
    for p in paths:
        if not os.path.exists(p): out['files'].append(p + ' (MISSING)'); continue
        out['files'].append(p)
        f = uproot.open(p); t = f['nuselection/NeutrinoSelectionFilter']
        out['entries'] += t.num_entries
        br = set(t.keys()); out['branches'] = br if out['branches'] is None else out['branches'] & br
        s = f['nuselection/SubRun']; sa = s.arrays([k for k in ('run', 'subRun', 'pot') if k in s.keys()], library='np')
        out['subruns'] += s.num_entries
        if 'pot' in sa: out['pot'] += float(sa['pot'].sum())
        rs |= set(zip(sa['run'].tolist(), sa['subRun'].tolist()))
        for wb in ('weightsPPFX', 'weightsGenie', 'weightsReint', 'weightsFlux'):
            if wb in br and wb not in out['wlen']:
                a = t[wb].array(entry_stop=200)
                n = ak.to_numpy(ak.num(a)); out['wlen'][wb] = int(n.max()) if len(n) else 0
        want = [b for b in EV + RECO + WGT if b in br]
        a = t.arrays(want, library='np')
        for b in want: cols.setdefault(b, []).append(a[b])
    for b in cols: cols[b] = np.concatenate(cols[b])
    return out, rs, cols


def compare(pair):
    name, ours, cust = pair
    try:
        A, ra, ca = info(ours); B, rb, cb = info(cust)
    except Exception as e:
        return dict(name=name, error=str(e))
    res = dict(name=name, ours={k: v for k, v in A.items() if k != 'branches'}, custom={k: v for k, v in B.items() if k != 'branches'})
    if A['branches'] is not None and B['branches'] is not None:
        res['only_ours'] = sorted(A['branches'] - B['branches']); res['only_custom'] = sorted(B['branches'] - A['branches'])
    res['subrun_common'] = len(ra & rb); res['subrun_only_ours'] = len(ra - rb); res['subrun_only_custom'] = len(rb - ra)
    if all(k in ca and k in cb for k in EV):
        ka = ca['run'].astype(np.int64) * 10 ** 12 + ca['sub'].astype(np.int64) * 10 ** 6 + ca['evt'].astype(np.int64)
        kb = cb['run'].astype(np.int64) * 10 ** 12 + cb['sub'].astype(np.int64) * 10 ** 6 + cb['evt'].astype(np.int64)
        ua, ia = np.unique(ka, return_index=True); ub, ib = np.unique(kb, return_index=True)
        com, xa, xb = np.intersect1d(ua, ub, return_indices=True)
        res.update(events_common=int(len(com)), events_only_ours=int(len(ua) - len(com)), events_only_custom=int(len(ub) - len(com)),
                   duplicates_ours=int(len(ka) - len(ua)), duplicates_custom=int(len(kb) - len(ub)))
        ja, jb = ia[xa], ib[xb]; agree = {}
        for b in RECO + WGT:
            if b in ca and b in cb and len(ja):
                va, vb = ca[b][ja].astype(float), cb[b][jb].astype(float)
                ok = np.isclose(va, vb, rtol=1e-5, atol=1e-6, equal_nan=True)
                agree[b] = dict(identical_fraction=float(ok.mean()),
                                mean_ours=float(np.nanmean(va)), mean_custom=float(np.nanmean(vb)),
                                median_ratio=float(np.nanmedian(np.where(vb != 0, va / np.where(vb != 0, vb, 1), np.nan))))
        res['common_event_agreement'] = agree
    print('done', name, flush=True)
    return res


if __name__ == '__main__':
    rw = '--reweighted' in sys.argv; args = [x for x in sys.argv[1:] if not x.startswith('--')]
    if rw: OUT = OUT.replace('.json', '_reweighted.json')
    with Pool(int(args[0]) if args else 5) as pool:
        res = pool.map(compare, PAIRS_RW if rw else PAIRS, chunksize=1)
    json.dump(res, open(OUT, 'w'), indent=1, default=float)
    for r in res:
        if 'error' in r: print(r['name'], 'ERROR', r['error']); continue
        print(f"{r['name']:20s} ours {r['ours']['entries']:9d} ev {r['ours']['subruns']:7d} sr pot {r['ours']['pot']:.3e} | custom {r['custom']['entries']:9d} ev "
              f"{r['custom']['subruns']:7d} sr pot {r['custom']['pot']:.3e} | common ev {r.get('events_common')} only-ours {r.get('events_only_ours')} "
              f"only-custom {r.get('events_only_custom')} | br only ours {len(r.get('only_ours', []))} only custom {len(r.get('only_custom', []))}")
