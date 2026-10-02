"""mp_pid_eval.py -- event-level comparison of the released pion identification with the particle
classifier (Phase 2 of report/MULTIPION_BNB_ADAPTATION_PLAN.md), for the single-pion selection and
the two multi-pion selections, on the run periods held out of the classifier training (FHC Run 5,
RHC Run 3).

Input: /data/uboone/processed/mp_pid_eval (slurm/slurm_mp_pid_eval.sbatch), one pass with every
selection and its *NewPID variant. Normalisation as the framework (scripts/mp_phase0_baseline.py):
MC and dirt by data POT over the distinct summed_pot, beam-off by 0.98 x triggers over gates, all
from the release file list.

    python3 scripts/mp_pid_eval.py [--no-ext] [--scan]   ->  report/multipion/phase2_eval[_scan][_noext].json, .md
--no-ext leaves out the beam-off samples (a preliminary while their processing runs); --scan reads the
threshold scan (/data/uboone/processed/mp_pid_scan, slurm/slurm_mp_pid_scan.sbatch) with every
"<NewPID selection>_tNN" variant (P(pi) > NN/100) next to its parent
"""
import json, os, sys
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_phase0_baseline as B

SCAN = '--scan' in sys.argv
B.NEW = '/data/uboone/processed/mp_pid_scan/' if SCAN else '/data/uboone/processed/mp_pid_eval/'
OUTDIR = B.OUTDIR
PERIODS = {'FHC Run 5': 5, 'RHC Run 3': 13}
PAIRS = [('CC1mu1piXp', 'CC1mu1piXpNewPID'), ('CC1mu2pi', 'CC1mu2piNewPID'), ('CC1mu3pi', 'CC1mu3piNewPID')]
SELS = [s for p in PAIRS for s in p]
if SCAN:
    SELS += ['CC1mu1piXpNewPID_t20', 'CC1mu1piXpNewPID_t40', 'CC1mu2piNewPID_t10', 'CC1mu2piNewPID_t15',
             'CC1mu2piNewPID_t30', 'CC1mu3piNewPID_t10', 'CC1mu3piNewPID_t30']
    SELS = sorted(SELS, key=lambda x: (x.split('NewPID')[0].replace('CC1mu1piXp', 'CC1mu1pi'), 'NewPID' in x, x))
WEIGHTS = ['spline_weight', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']
ID = ['n_reco_pion_truepi', 'n_reco_pion_proton', 'n_reco_pion_muon', 'n_reco_pion_other']


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def read(path, kind):
    t = uproot.open(path)['stv_tree']
    br = [f'{s}_{v}' for s in SELS for v in ['Selected', 'MC_Signal'] + ID]
    if kind != 'ext':
        br += WEIGHTS
    a = t.arrays([b for b in br if b in t.keys()], library='np')
    w = np.ones(t.num_entries) if kind == 'ext' else safe(np.prod([a[x].astype(float) for x in WEIGHTS], axis=0))
    return a, w


def main():
    no_ext = '--no-ext' in sys.argv
    tag = ('_scan' if SCAN else '') + ('_noext' if no_ext else '')
    on, ext, mc, dirt = B.parse_conf()
    keep = set(PERIODS.values())
    ext, mc, dirt = ({r: v for r, v in tab.items() if r in keep} for tab in (ext, mc, dirt))
    sc, checks = B.scales(on, ext, mc, dirt)
    bad = [c for c in checks if abs(c[3] - c[4]) > 1e-5 * abs(c[4])]
    if bad: sys.exit(f'summed_pot differs from the release copies: {bad}')
    acc = {}       # (period, sel) -> sums
    both = {}      # (period, pair) -> overlap of the selected samples
    for pname, run in PERIODS.items():
        for kind in (('mc', 'dirt') if no_ext else ('mc', 'ext', 'dirt')):
            for path, s in sc[(run, kind)]:
                a, w = read(path, kind)
                ww = w * s
                print(f'  {pname} {kind}: {os.path.basename(path)} scale {s:.4f}', flush=True)
                for sel in SELS:
                    d = acc.setdefault((pname, sel), dict(gen=0., sig=0., bkg=0., ext=0., dirt=0., bkg_by_id=np.zeros(5)))
                    selc = a[f'{sel}_Selected'].astype(bool)
                    if kind == 'mc':
                        sig = a[f'{sel}_MC_Signal'].astype(bool)
                        d['gen'] += ww[sig].sum(); d['sig'] += ww[selc & sig].sum(); d['bkg'] += ww[selc & ~sig].sum()
                        # true identity of the counted pion(s) in selected background: proton, muon, other,
                        # true pion only, none recorded
                        ids = np.stack([a[f'{sel}_{x}'] for x in ID], axis=1)
                        m = selc & ~sig
                        cls = np.where(ids[:, 1] > 0, 0, np.where(ids[:, 2] > 0, 1, np.where(ids[:, 3] > 0, 2, np.where(ids[:, 0] > 0, 3, 4))))
                        for c in range(5): d['bkg_by_id'][c] += ww[m & (cls == c)].sum()
                    elif kind == 'ext':
                        d['ext'] += ww[selc].sum()
                    else:
                        d['dirt'] += ww[selc].sum()
                for old, new in PAIRS:
                    o = both.setdefault((pname, old), dict(both_sig=0., old_only_sig=0., new_only_sig=0.,
                                                           both_bkg=0., old_only_bkg=0., new_only_bkg=0.))
                    so, sn = a[f'{old}_Selected'].astype(bool), a[f'{new}_Selected'].astype(bool)
                    sig = a[f'{old}_MC_Signal'].astype(bool) if kind == 'mc' else np.zeros(len(so), bool)
                    for part, m in (('both', so & sn), ('old_only', so & ~sn), ('new_only', ~so & sn)):
                        o[part + '_sig'] += ww[m & sig].sum(); o[part + '_bkg'] += ww[m & ~sig].sum()
    out = {}
    for pname in list(PERIODS) + ['both periods']:
        for sel in SELS:
            if pname == 'both periods':
                parts = [acc[(p, sel)] for p in PERIODS]
                d = {k: (sum(p[k] for p in parts)) for k in ('gen', 'sig', 'bkg', 'ext', 'dirt')}
                d['bkg_by_id'] = sum(p['bkg_by_id'] for p in parts)
            else:
                d = acc[(pname, sel)]
            tot = d['sig'] + d['bkg'] + d['ext'] + d['dirt']
            out.setdefault(pname, {})[sel] = dict(
                generated=d['gen'], signal=d['sig'], nu_bkg=d['bkg'], ext=d['ext'], dirt=d['dirt'], total=tot,
                efficiency=d['sig'] / d['gen'] if d['gen'] else float('nan'), purity=d['sig'] / tot if tot else float('nan'),
                stat_rel_unc_one_bin=np.sqrt(tot) / d['sig'] if d['sig'] else float('nan'),
                bkg_by_counted_pion_id=dict(zip(['proton', 'muon', 'other', 'true pion only', 'not recorded'],
                                                [float(x) for x in d['bkg_by_id']])))
    out['overlap'] = {f'{p} {o}': v for (p, o), v in both.items()}
    out['beam_off_included'] = not no_ext
    json.dump(out, open(os.path.join(OUTDIR, f'phase2_eval{tag}.json'), 'w'), indent=1, default=float)
    L = ['# Phase 2: released pion identification against the particle classifier, event level', '',
         'Held-out run periods of the classifier (FHC Run 5, RHC Run 3), normalised to their data exposure; '
         + ('PRELIMINARY: beam-off NOT included (its processing was still running); dirt included. ' if no_ext else 'Beam-off and dirt included. ')
         + 'Background split by the true identity of the counted pion candidate(s).', '']
    for pname in list(PERIODS) + ['both periods']:
        L += [f'## {pname}', '', '| Selection | Signal | ν bkg | Beam-off | Dirt | Efficiency | Purity | stat. unc. (one bin) | bkg: p as π | μ as π | other | true π |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for sel in SELS:
            r = out[pname][sel]; b = r['bkg_by_counted_pion_id']
            L.append(f'| {sel} | {r["signal"]:.1f} | {r["nu_bkg"]:.1f} | {r["ext"]:.1f} | {r["dirt"]:.1f} | {100 * r["efficiency"]:.2f}% | '
                     f'{100 * r["purity"]:.1f}% | {100 * r["stat_rel_unc_one_bin"]:.2f}% | {b["proton"]:.1f} | {b["muon"]:.1f} | {b["other"]:.1f} | {b["true pion only"]:.1f} |')
        L.append('')
    L += ['## Overlap of the selected samples (MC signal / MC background + beam-off + dirt)', '']
    for k, v in out['overlap'].items():
        L.append(f'- {k}: both {v["both_sig"]:.1f} / {v["both_bkg"]:.1f}; released only {v["old_only_sig"]:.1f} / {v["old_only_bkg"]:.1f}; '
                 f'new only {v["new_only_sig"]:.1f} / {v["new_only_bkg"]:.1f}')
    open(os.path.join(OUTDIR, f'phase2_eval{tag}.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
