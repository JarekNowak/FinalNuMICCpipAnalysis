"""bragg_cut_eval.py -- event-level effect of the Bragg-pion >= 0.08 cut of the single-pion
selections (report/multipion/PHASE2_SUMMARY.md).

The beam-on, beam-off and dirt ntuples of the older production (FHC Run 1, RHC Runs 1 and 3b) lack
trk_bragg_pion_v, so there the cut passes every track, while every overlay applies it. At unblinding
the data of those periods would be selected without the cut and compared with an overlay prediction
that has it (their beam-off and dirt are already without it). This compares the released selections
with their NoBragg variants in the reprocessed release samples (/data/uboone/processed/bragg_eval,
slurm/slurm_bragg_eval.sbatch), normalised as the framework (scripts/mp_phase0_baseline.py), per run
period and per configuration:

  - the expected data excess of an older-production period: the overlay difference NoBragg - released
    over the released total prediction;
  - the bias of the extracted one-bin cross section of a configuration: the same overlay difference,
    summed over its older-production periods, over the released signal (the background subtraction
    and the efficiency both come from the overlay with the cut);
  - dropping the cut in every sample: efficiency, purity and the one-bin statistical uncertainty.

    python3 scripts/bragg_cut_eval.py  ->  report/multipion/bragg_cut_eval.{json,md}
"""
import json, os, sys
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_phase0_baseline as B

B.NEW = '/data/uboone/processed/bragg_eval/'
OUTDIR = B.OUTDIR
PAIRS = [('CC1mu1piXp', 'CC1mu1piXpNoBragg'), ('CC1mu1pi1p', 'CC1mu1pi1pNoBragg')]
SELS = [s for p in PAIRS for s in p]
OLD = (1, 11, 13)                     # runs whose beam-on data come from the older production
RUN2 = (2, 12)                        # no beam-on data in hand
PNAME = {1: 'FHC Run 1', 2: 'FHC Run 2', 4: 'FHC Run 4', 5: 'FHC Run 5',
         11: 'RHC Run 1', 12: 'RHC Run 2', 13: 'RHC Run 3', 14: 'RHC Run 4'}
WEIGHTS = ['spline_weight', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']


def safe(w):
    return np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.)


def main():
    on, ext, mc, dirt = B.parse_conf()
    sc, checks = B.scales(on, ext, mc, dirt)
    bad = [c for c in checks if abs(c[3] - c[4]) > 1e-5 * abs(c[4])]
    if bad: sys.exit(f'summed_pot differs from the release copies: {bad}')
    acc = {}                          # (run, sel) -> sums
    for (run, kind), files in sorted(sc.items()):
        for path, s in files:
            t = uproot.open(path)['stv_tree']
            br = [f'{x}_{v}' for x in SELS for v in ('Selected', 'MC_Signal')] + ([] if kind == 'ext' else WEIGHTS)
            a = t.arrays(br, library='np')
            w = s * (np.ones(t.num_entries) if kind == 'ext' else safe(np.prod([a[x].astype(float) for x in WEIGHTS], axis=0)))
            for sel in SELS:
                d = acc.setdefault((run, sel), dict(gen=0., sig=0., bkg=0., ext=0., dirt=0.))
                m = a[f'{sel}_Selected'].astype(bool)
                if kind == 'mc':
                    sig = a[f'{sel}_MC_Signal'].astype(bool)
                    d['gen'] += w[sig].sum(); d['sig'] += w[m & sig].sum(); d['bkg'] += w[m & ~sig].sum()
                else:
                    d[kind] += w[m].sum()
            print(f'  run {run} {kind}: {os.path.basename(path)} scale {s:.4f}', flush=True)
    runs = sorted({r for r, _ in acc})
    tot = lambda d: d['sig'] + d['bkg'] + d['ext'] + d['dirt']
    out = {'per_run': {}, 'config': {}}
    for r in runs:
        for rel, nb in PAIRS:
            a, b = acc[(r, rel)], acc[(r, nb)]
            dmc = (b['sig'] + b['bkg']) - (a['sig'] + a['bkg'])          # overlay change only
            out['per_run'][f'{PNAME[r]} {rel}'] = dict(
                released=a, nobragg=b, overlay_change=dmc, excess_over_prediction=dmc / tot(a),
                purity_released=a['sig'] / tot(a), purity_nobragg=b['sig'] / tot(b),
                signal_change=b['sig'] / a['sig'] - 1, total_change=tot(b) / tot(a) - 1)
    for cname, cruns in B.CONFIGS.items():
        for rel, nb in PAIRS:
            S = lambda sel, k: sum(acc[(r, sel)][k] for r in cruns if (r, sel) in acc)
            a = {k: S(rel, k) for k in ('gen', 'sig', 'bkg', 'ext', 'dirt')}
            b = {k: S(nb, k) for k in ('gen', 'sig', 'bkg', 'ext', 'dirt')}
            dold = sum((acc[(r, nb)]['sig'] + acc[(r, nb)]['bkg']) - (acc[(r, rel)]['sig'] + acc[(r, rel)]['bkg'])
                       for r in cruns if r in OLD)
            drun2 = sum((acc[(r, nb)]['sig'] + acc[(r, nb)]['bkg']) - (acc[(r, rel)]['sig'] + acc[(r, rel)]['bkg'])
                        for r in cruns if r in RUN2)
            out['config'][f'{cname} {rel}'] = dict(
                released=a, nobragg=b,
                bias_onebin_old_periods=dold / a['sig'], bias_onebin_if_run2_old=(dold + drun2) / a['sig'],
                eff_released=a['sig'] / a['gen'], eff_nobragg=b['sig'] / b['gen'],
                purity_released=a['sig'] / tot(a), purity_nobragg=b['sig'] / tot(b),
                stat_released=np.sqrt(tot(a)) / a['sig'], stat_nobragg=np.sqrt(tot(b)) / b['sig'],
                bs_released=(tot(a) - a['sig']) / a['sig'], bs_nobragg=(tot(b) - b['sig']) / b['sig'])
    json.dump(out, open(os.path.join(OUTDIR, 'bragg_cut_eval.json'), 'w'), indent=1, default=float)
    L = ['# The Bragg-pion cut of the single-pion selections at event level', '',
         'Release samples at data exposure (overlay, beam-off, dirt). "NoBragg" drops Bragg-pion >= 0.08 in every sample. '
         'Older-production periods (beam-on, beam-off and dirt without the branch): FHC Run 1, RHC Runs 1 and 3; '
         'Run 2 has no beam-on data in hand.', '',
         '## Per period', '',
         '| Period | Selection | Signal released / NoBragg | Total released / NoBragg | Purity released / NoBragg | Overlay change / prediction |',
         '|---|---|---|---|---|---|']
    for k, v in out['per_run'].items():
        a, b = v['released'], v['nobragg']
        L.append(f'| {k.rsplit(" ", 1)[0]} | {k.rsplit(" ", 1)[1]} | {a["sig"]:.1f} / {b["sig"]:.1f} | {tot(a):.1f} / {tot(b):.1f} | '
                 f'{100 * v["purity_released"]:.1f} / {100 * v["purity_nobragg"]:.1f}% | {100 * v["excess_over_prediction"]:+.1f}% |')
    L += ['', '## Per configuration', '',
          '| Configuration | Selection | Efficiency released / NoBragg | Purity released / NoBragg | B/S released / NoBragg | Stat. unc. released / NoBragg | One-bin bias, older-production periods (with Run 2) |',
          '|---|---|---|---|---|---|---|']
    for k, v in out['config'].items():
        L.append(f'| {k.split(" ")[0]} | {k.split(" ")[1]} | {100 * v["eff_released"]:.2f} / {100 * v["eff_nobragg"]:.2f}% | '
                 f'{100 * v["purity_released"]:.1f} / {100 * v["purity_nobragg"]:.1f}% | {v["bs_released"]:.2f} / {v["bs_nobragg"]:.2f} | '
                 f'{100 * v["stat_released"]:.2f} / {100 * v["stat_nobragg"]:.2f}% | {100 * v["bias_onebin_old_periods"]:+.1f}% '
                 f'({100 * v["bias_onebin_if_run2_old"]:+.1f}%) |')
    open(os.path.join(OUTDIR, 'bragg_cut_eval.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
