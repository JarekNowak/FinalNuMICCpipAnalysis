"""braggfix_gate.py -- staging against live, before the promotion of the trees reprocessed without the
Bragg-pion cut (commit 1fb6be9). For every regular file of processed/{incl_braggfix,w_braggfix,
ext_perrun_braggfix} and its live counterpart:

  - same entry count, same summed_pot, same branch names;
  - event identity and central-value weights bit-identical (the tree has no event numbers: tuned_cv_weight,
    ppfx_cv_weight, normalisation_weight, spline_weight, mc_nu_energy, mc_nu_pdg, topological_score,
    nslice), so the same-seed throws reproduce the same fake-data events;
  - the truth signal definitions (<sel>_MC_Signal) unchanged;
  - no live branch missing (branches added since the live files were made are reported);
  - provenance, where /data/uboone/processed/bragg_eval holds the same input (overlays, beam-off, dirt;
    scripts/bragg_cut_eval.py on multipion/phase1): every <sel>_* branch of the live file equals the
    current code WITH the cut (<sel>_*), and every <sel>_* branch of the staging file equals the current
    code without it (<sel>NoBragg_*), so the cut is the only change;
  - reported, not required: the fraction of events whose <sel>_Selected or control-region flags changed,
    and the selected count before and after.

    python3 slurm/braggfix_gate.py  ->  ../logs/braggfix/gate.tsv; exit status 1 if any required check fails
"""
import os, sys
import numpy as np
import uproot

PROC = '/data/uboone/processed'
PAIRS = [(PROC + '/incl_braggfix', PROC), (PROC + '/w_braggfix', PROC + '/w'), (PROC + '/ext_perrun_braggfix', PROC + '/ext_perrun')]
# the processed tree carries no event numbers: identity = CV weights + selection-independent reco and truth
IDW = ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight', 'spline_weight', 'mc_nu_energy', 'mc_nu_pdg',
       'topological_score', 'nslice']
SELS = ['CC1mu1piXp', 'CC1mu1pi1p']
FLAGS = ['sb_cc0pi', 'sb_pi0', 'sb_multipi', 'sb_cosmic']
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'logs', 'braggfix')
EVAL = PROC + '/bragg_eval/'
# release dirt names that are regular files -> name of the processed raw sample in bragg_eval
DIRT_ALIAS = {
    'xsec-ana-dirt_fhc_run4c.root': 'xsec-ana-numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana.root',
    'xsec-ana-dirt_fhc_run4d.root': 'xsec-ana-numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana.root',
    'xsec-ana-dirt_fhc_run5.root': 'xsec-ana-run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns.root',
    'xsec-ana-dirt_rhc_run3b.root': 'xsec-ana-neutrinoselection_filt_run3b_dirt_overlay.root',
    'xsec-ana-dirt_rhc_run4a.root': 'xsec-ana-run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana.root',
    'xsec-ana-dirt_rhc_run4b.root': 'xsec-ana-numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana.root',
}


def same(a, b):
    if a.shape != b.shape: return False
    return np.array_equal(a, b, equal_nan=True) if a.dtype.kind == 'f' else np.array_equal(a, b)


def check(new, old):
    fn, fo = uproot.open(new), uproot.open(old)
    tn, to = fn['stv_tree'], fo['stv_tree']
    res = dict(file=new, ok=True, notes=[])
    def fail(msg): res['ok'] = False; res['notes'].append(msg)
    if tn.num_entries != to.num_entries: fail(f'entries {tn.num_entries} vs {to.num_entries}')
    pn, po = float(fn['summed_pot'].member('fVal')), float(fo['summed_pot'].member('fVal'))
    if pn != po: fail(f'summed_pot {pn} vs {po}')
    kn, ko = set(tn.keys()), set(to.keys())
    if ko - kn: fail(f'live branches missing: {sorted(ko - kn)[:5]}')
    if kn - ko: res['notes'].append(f'added: {",".join(sorted(kn - ko))}')
    if not res['ok']: return res
    ids = [b for b in IDW if b in kn]
    sels = [s for s in SELS if f'{s}_Selected' in kn]
    want = ids + [f'{s}_{v}' for s in sels for v in ['Selected', 'MC_Signal'] + FLAGS if f'{s}_{v}' in kn]
    acc = {s: np.zeros(4) for s in sels}           # selected old, selected new, changed Selected, changed CR flag
    for an, ao in zip(tn.iterate(want, step_size='300 MB', library='np'), to.iterate(want, step_size='300 MB', library='np')):
        for b in ids:
            if not same(an[b], ao[b]): fail(f'{b} differs')
        for s in sels:
            if not same(an[f'{s}_MC_Signal'], ao[f'{s}_MC_Signal']): fail(f'{s}_MC_Signal differs')
            sn, so = an[f'{s}_Selected'].astype(bool), ao[f'{s}_Selected'].astype(bool)
            crn = np.zeros(len(sn), bool); cro = np.zeros(len(sn), bool)
            for v in FLAGS:
                if f'{s}_{v}' in an: crn |= an[f'{s}_{v}'].astype(bool); cro |= ao[f'{s}_{v}'].astype(bool)
            acc[s] += [so.sum(), sn.sum(), (sn != so).sum(), (crn != cro).sum()]
        if not res['ok']: break
    res['entries'] = tn.num_entries
    res['sel'] = {s: [int(x) for x in v] for s, v in acc.items()}
    ev = EVAL + DIRT_ALIAS.get(os.path.basename(new), os.path.basename(new))
    if res['ok'] and os.path.exists(ev) and not os.path.basename(new).startswith('xsec-ana-fakedata'):
        te = uproot.open(ev)['stv_tree']; ke = set(te.keys())
        if te.num_entries != tn.num_entries: fail(f'bragg_eval entries {te.num_entries}')
        else:
            nb = 0
            for s in sels:
                live_b = [b for b in ko if b.startswith(s + '_') and b in ke]
                stag_b = [b for b in kn if b.startswith(s + '_') and f'{s}NoBragg_{b[len(s) + 1:]}' in ke]
                for b in live_b:
                    if not same(to[b].array(library='np'), te[b].array(library='np')): fail(f'live {b} != current code with the cut')
                for b in stag_b:
                    if not same(tn[b].array(library='np'), te[f'{s}NoBragg_{b[len(s) + 1:]}'].array(library='np')): fail(f'staging {b} != NoBragg')
                nb += len(live_b) + len(stag_b)
            res['notes'].append(f'provenance: {nb} branches checked against bragg_eval')
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, bad = [], 0
    for stage, live in PAIRS:
        for b in sorted(os.listdir(stage)):
            p = os.path.join(stage, b)
            if not b.endswith('.root') or os.path.islink(p): continue
            r = check(p, os.path.join(live, b))
            bad += not r['ok']
            sel = '; '.join(f'{s}: selected {v[0]} -> {v[1]}, changed {v[2]}, CR flag changed {v[3]}' for s, v in r.get('sel', {}).items())
            rows.append(f"{'OK' if r['ok'] else 'FAIL'}\t{os.path.relpath(p, PROC)}\t{r.get('entries', '')}\t{sel}\t{' | '.join(r['notes'])}")
            print(rows[-1], flush=True)
    open(os.path.join(OUT, 'gate.tsv'), 'w').write('\n'.join(rows) + '\n')
    print(f'{len(rows) - bad} / {len(rows)} files pass')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
