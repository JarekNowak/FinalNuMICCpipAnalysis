"""ext_gates_table.py -- beam-off (EXT) normalisation per beam-on run period, from the live configuration.

Reads the onBNB (beam-on triggers, POT) and extBNB (beam-off gates) entries of the live FHC and RHC file
lists, counts the selected beam-off events in each processed beam-off file (inclusive
CC1mu1piXp_Selected; proton-tagged CC1mu1pi1p_Selected from the w/ files), and writes
data_release/ext_gates.tsv and figures/ext_perrun_validation.pdf. Files of one run are pooled, as the
framework does: scale = 0.98 * beam-on triggers / sum(beam-off gates).

    python3 report/tools/ext_gates_table.py
"""
import os, collections, uproot, numpy as np
R = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/'
CFG = R + 'xsec_analyzer/configs/'
OCC = 0.98


def entries(fp):
    on, ext = {}, collections.defaultdict(list)
    for line in open(fp):
        p = line.split()
        if len(p) < 4 or p[0].startswith('#'):
            continue
        if p[2] == 'onBNB':
            on[int(p[1])] = (float(p[3]), float(p[4]))
        elif p[2] == 'extBNB':
            ext[int(p[1])].append((p[0], float(p[3])))
    return on, ext


def nsel(path, flag):
    t = uproot.open(path)['stv_tree']
    return int(np.sum(t[flag].array(library='np')))


rows = []
for mode, fp, fpw in [('FHC', 'file_properties_numi_fhc5.txt', 'file_properties_numi_fhc5_w.txt'),
                      ('RHC', 'file_properties_numi_rhcfull.txt', 'file_properties_numi_rhcfull_w.txt')]:
    on, ext = entries(CFG + fp)
    _, extw = entries(CFG + fpw)
    for run in sorted(on):
        trig, pot = on[run]
        files = ext[run]
        gates = sum(g for _, g in files)
        assert abs(gates - sum(g for _, g in extw[run])) < 1e-3, (mode, run)
        scale = OCC * trig / gates
        raw = sum(nsel(f, 'CC1mu1piXp_Selected') for f, _ in files)
        raw1p = sum(nsel(f, 'CC1mu1pi1p_Selected') for f, _ in extw[run])
        # processed names are xsec-ana-ext_<sample>_run<beam-on run>.root; <sample> names the period recorded
        samples = [os.path.basename(f)[len('xsec-ana-ext_'):].rsplit('_run', 1)[0] for f, _ in files]
        tags = '+'.join(samples)
        if not samples[0].startswith(f'run{run}'):
            tags += f' (stand-in for run {run})'
        rows.append(dict(mode=mode, run=run, pot=pot, trig=trig, tags=tags, gates=gates, scale=scale,
                         raw=raw, scaled=raw * scale, raw1p=raw1p, scaled1p=raw1p * scale, rate=raw / gates * 1e6))

out = R + 'report/data_release/ext_gates.tsv'
with open(out, 'w') as o:
    o.write('# Beam-off (EXT) normalisation per beam-on run period, from the live FHC and RHC file lists '
            '(tools/ext_gates_table.py).\n# Beam-off data are cosmic-only and are matched by run period irrespective of '
            'horn polarity; files within a run are pooled. scale = 0.98 * bnb_gates / sum(ext_gates).\n')
    o.write('mode\trun\tbnb_pot\tbnb_gates\text_files\text_gates_sum\tscale\text_sel_incl_raw\text_sel_incl_scaled'
            '\text_sel_1p_raw\text_sel_1p_scaled\trate_per_Mgate_incl\n')
    for r in rows:
        o.write(f"{r['mode']}\t{r['run']}\t{r['pot']:.4e}\t{r['trig']:.0f}\t{r['tags']}\t{r['gates']:.2f}\t{r['scale']:.5f}"
                f"\t{r['raw']}\t{r['scaled']:.2f}\t{r['raw1p']}\t{r['scaled1p']:.2f}\t{r['rate']:.2f}\n")
    tot = {m: sum(r['scaled'] for r in rows if r['mode'] == m) for m in ('FHC', 'RHC')}
    tot1 = {m: sum(r['scaled1p'] for r in rows if r['mode'] == m) for m in ('FHC', 'RHC')}
    o.write(f"# totals at the final cut, inclusive: FHC {tot['FHC']:.1f}, RHC {tot['RHC']:.1f}, combined "
            f"{tot['FHC'] + tot['RHC']:.1f}; proton-tagged: FHC {tot1['FHC']:.1f}, RHC {tot1['RHC']:.1f}, combined "
            f"{tot1['FHC'] + tot1['RHC']:.1f}\n")
print(open(out).read())

import sys
sys.path.insert(0, R + 'report/tools')
import xsec_figs_style as S
import matplotlib.pyplot as plt
S.apply()
lab = [f"{r['mode']} {r['run']}" for r in rows]
x = np.arange(len(rows))
fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
rate = np.array([r['rate'] for r in rows])
rerr = np.array([np.sqrt(r['raw']) / r['gates'] * 1e6 for r in rows])
ax[0].errorbar(x, rate, yerr=rerr, fmt='o', color='black', ms=4, capsize=2)
ax[0].set_ylabel(r'selected beam-off events per $10^6$ gates')
ax[1].bar(x, [r['scale'] for r in rows], color='#999999', edgecolor='black', lw=0.6)
ax[1].set_ylabel(r'scale $=0.98\,N_{\rm gates}^{\rm on}/N_{\rm gates}^{\rm off}$')
ev = np.array([r['scaled'] for r in rows])
everr = np.array([np.sqrt(r['raw']) * r['scale'] for r in rows])
ax[2].bar(x, ev, yerr=everr, color='#999999', edgecolor='black', lw=0.6, capsize=2)
ax[2].set_ylabel('beam-off events at the final cut')
for a in ax:
    a.set_xticks(x)
    a.set_xticklabels(lab, rotation=45, ha='right', fontsize=8)
    a.axvline(3.5, color='black', lw=0.5, ls=':')
    a.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(R + 'report/figures/ext_perrun_validation.pdf')
print('wrote figures/ext_perrun_validation.pdf')
