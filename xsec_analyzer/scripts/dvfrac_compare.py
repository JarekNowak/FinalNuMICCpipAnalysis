#!/usr/bin/env python3
"""dvfrac_compare.py -- tables of scripts/dvfrac_impact.sh: for every released single-pion product the integrated
cross section, its detector and total uncertainty, and the per-bin detector and total uncertainties, as released
("base", detVar as DV) and with the fractional treatment ("dvfrac"); plus the check that "base" reproduces the
release (SYSTDUMP of logs/systdump, or of the harvest logs for the products outside regen_fdfix.sh).

    python3 scripts/dvfrac_compare.py [OUT]       (default /data/uboone/processed/dvfrac_study)
Writes OUT/dvfrac_impact.md and OUT/dvfrac_impact.tsv.
"""
import os, re, sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else '/data/uboone/processed/dvfrac_study'
MAIN = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer'
CFG = {'FHC5': 'fhc5', 'RHCFULL': 'rhcfull', 'COMB': 'comb'}


def systdump(path):
    d = {}
    if not os.path.exists(path): return None
    for line in open(path):
        m = re.match(r'\[SYSTDUMP\]\s+(\S+)\s+(\S+)', line)
        if m: d[m.group(1)] = float(m.group(2))
    return d or None


def release_dump(tag):
    fam, t, obs = tag.split('_', 2)
    cfg = CFG[t]
    o = 'ppi' if obs == 'ppi3bin' else obs
    name = f'{cfg}_{o}' if fam == 'incl' else f'1p_{cfg}_{o}'
    d = systdump(f'{MAIN}/logs/systdump/{name}.dump')
    return d if d else systdump(f'{MAIN}/logs/fdfix/harvest_proposed_{tag}.log')


def table(path):
    rows = [l.split() for l in open(path)]
    n = int(rows[0][1])
    if rows[1][0] == 'numYbins':
        m = np.zeros((n, n))
        for r in rows[3:]:
            if len(r) == 3: m[int(r[0]), int(r[1])] = float(r[2])
        return m
    v = np.zeros(n)
    for r in rows[1:]:
        if len(r) == 2: v[int(r[0])] = float(r[1])
    return v


def closure(d, tag):
    """truth chi2 / bins of the framework (fake-data truth through A_C) and the trace of A_C."""
    chi, nb = np.nan, 0
    for line in open(f'{d}/run.log'):
        m = re.match(r'truth: \S+ = (\S+)/(\d+)', line)
        if m: chi, nb = float(m.group(1)), int(m.group(2)); break
    import uproot
    tr = np.trace(uproot.open(f'{d}/closure_hists_xsec_{tag}.root')['h_A_C'].values())
    return chi, nb, tr


def perbin(d):
    x = table(f'{d}/unfolded_signal.txt')
    det = np.sqrt(np.diag(table(f'{d}/cov_detVar_total.txt'))) / np.abs(x)
    tot = np.sqrt(np.diag(table(f'{d}/cov_total.txt'))) / np.abs(x)
    return x, det, tot


def main():
    tags = [l.split('|')[0] for l in open(f'{OUT}/products.txt') if l.strip()]
    rows, bad = [], []
    for tag in tags:
        b, f = systdump(f'{OUT}/base/{tag}/systdump.txt'), systdump(f'{OUT}/dvfrac/{tag}/systdump.txt')
        if not b or not f:
            bad.append(f'{tag}: missing output'); continue
        rel = release_dump(tag)
        same = 'n/a' if rel is None else ('yes' if all(abs(rel[k] - b.get(k, np.nan)) <= 1e-3 * abs(rel[k]) + 1e-9
                                                        for k in rel) else 'NO')
        xb, db, tb = perbin(f'{OUT}/base/{tag}')
        xf, df, tf = perbin(f'{OUT}/dvfrac/{tag}')
        cb, cf = closure(f'{OUT}/base/{tag}', tag), closure(f'{OUT}/dvfrac/{tag}', tag)
        rows.append(dict(tag=tag, same=same, nb=len(xb), chi_b=cb[0], chi_f=cf[0], ndf=cb[1], tr_b=cb[2], tr_f=cf[2], sig_b=b['sigma_int'], sig_f=f['sigma_int'],
                         det_b=b['detVar_total'], det_f=f['detVar_total'], tot_b=b['total'], tot_f=f['total'],
                         bdet_b=db, bdet_f=df, btot_b=tb, btot_f=tf, x_ratio=xf / xb))
    with open(f'{OUT}/dvfrac_impact.tsv', 'w') as fp:
        fp.write('tag\tbase_reproduces_release\tnbins\tsigma_base\tsigma_dvfrac\tdet_base\tdet_dvfrac\ttot_base\t'
                 'tot_dvfrac\tbin_det_base\tbin_det_dvfrac\tbin_tot_base\tbin_tot_dvfrac\tbin_xsec_ratio\n')
        for r in rows:
            fmt = lambda a: ','.join(f'{v:.4f}' for v in a)
            fp.write(f"{r['tag']}\t{r['same']}\t{r['nb']}\t{r['sig_b']:.5g}\t{r['sig_f']:.5g}\t{r['det_b']:.3f}\t"
                     f"{r['det_f']:.3f}\t{r['tot_b']:.3f}\t{r['tot_f']:.3f}\t{fmt(r['bdet_b'])}\t{fmt(r['bdet_f'])}\t"
                     f"{fmt(r['btot_b'])}\t{fmt(r['btot_f'])}\t{fmt(r['x_ratio'])}\n")
    L = ['# Fractional detector variations: impact on the released single-pion results', '',
         f'{len(rows)} products ({OUT}). Integrated cross section and its detector and total uncertainty from the '
         'SYSTDUMP lines (percent); per bin, the detector and total uncertainty relative to the unfolded value, as '
         'the ratio dvfrac / base (median and range over bins), and the ratio of the unfolded values; the framework\'s '
         'chi2 of the unfolded pseudo-data against its truth through A_C, and the trace of A_C. "Base" is the '
         'release configuration run with this checkout ("Release" = it reproduces the released SYSTDUMP).', '',
         '| Product | Release | sigma dvfrac/base | Detector % base | dvfrac | Total % base | dvfrac | '
         'Bin detector ratio | Bin total ratio | Bin xsec ratio | Truth chi2 base | dvfrac | tr A_C base | dvfrac |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    rng = lambda a: f'{np.median(a):.2f} ({a.min():.2f}-{a.max():.2f})' if len(a) > 1 else f'{a[0]:.2f}'
    for r in rows:
        L.append(f"| {r['tag']} | {r['same']} | {r['sig_f'] / r['sig_b']:.3f} | {r['det_b']:.1f} | {r['det_f']:.1f} | "
                 f"{r['tot_b']:.1f} | {r['tot_f']:.1f} | {rng(r['bdet_f'] / r['bdet_b'])} | "
                 f"{rng(r['btot_f'] / r['btot_b'])} | {rng(r['x_ratio'])} | "
                 + (f"{r['chi_b']:.2f}/{r['ndf']} | {r['chi_f']:.2f}/{r['ndf']} | " if r['ndf'] else '- | - | ')
                 + f"{r['tr_b']:.2f} | {r['tr_f']:.2f} |")
    if rows:
        alldet = np.concatenate([r['bdet_f'] / r['bdet_b'] for r in rows])
        alltot = np.concatenate([r['btot_f'] / r['btot_b'] for r in rows])
        allx = np.concatenate([r['x_ratio'] for r in rows])
        L += ['', f'All {len(alldet)} bins: detector ratio median {np.median(alldet):.2f} '
              f'(range {alldet.min():.2f}-{alldet.max():.2f}); total ratio median {np.median(alltot):.2f} '
              f'({alltot.min():.2f}-{alltot.max():.2f}); unfolded value ratio median {np.median(allx):.3f} '
              f'({allx.min():.3f}-{allx.max():.3f}).']
        for t in CFG:
            sel = [r for r in rows if f'_{t}_' in r['tag']]
            d = np.concatenate([r['bdet_f'] / r['bdet_b'] for r in sel])
            L.append(f'{t}: detector ratio median {np.median(d):.2f} over {len(d)} bins; integrated detector '
                     f'{np.median([r["det_b"] for r in sel]):.1f}% -> {np.median([r["det_f"] for r in sel]):.1f}% '
                     f'(medians), total {np.median([r["tot_b"] for r in sel]):.1f}% -> '
                     f'{np.median([r["tot_f"] for r in sel]):.1f}%.')
    if bad: L += ['', 'Problems: ' + '; '.join(bad)]
    open(f'{OUT}/dvfrac_impact.md', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
