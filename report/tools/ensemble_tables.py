"""ensemble_tables.py -- evaluate the six configuration-aware fake-data ensembles of version 1.7 (2026-10-05,
selection without the Bragg-pion requirement, xsec_analyzer/ens_braggfix.sh: FHC p_mu, FHC and RHC cos theta_mu, COMB p_pi
(three regions), FHC and RHC one-bin totals; 100 members each) with ensemble_stat_pulls.evaluate, write
data_release/ensemble_2026-10-05.tsv and
patch the coverage tables of the analysis note (statistical reference) and the technical supplement
(statistical + full). Run slurm/ens_statcov_all.sh first so every member has its DataStats table.

    python3 report/tools/ensemble_tables.py
"""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ensemble_stat_pulls import evaluate
R=os.path.abspath(os.path.join(HERE,'..'))+'/'
TAGS=[('fhc5_pmu',r'$d\sigma/dp_\mu$, FHC'),('fhc5_costhmu',r'$d\sigma/d\cos\theta_\mu$, FHC'),('rhcfull_costhmu',r'$d\sigma/d\cos\theta_\mu$, RHC'),
      ('comb_ppi3bin',r'$d\sigma/dp_\pi$ (three regions), combined'),('fhc5_total','one-bin total, FHC'),('rhcfull_total','one-bin total, RHC')]
res={}
for tag,lab in TAGS:
    try: res[tag]=evaluate(tag, 100, quiet=True)
    except Exception as e: print(tag,'failed:',e)
with open(R+'data_release/ensemble_2026-10-05.tsv','w') as o:
    o.write('# Fake-data ensembles of the version 1.7 release (2026-10-05; selection without the Bragg-pion requirement): Poisson throws of the central-value\n'
            '# prediction of each configuration (slurm_ensemble_cfg.sbatch), each carried through univmake + UnfolderNuMI; pulls against\n'
            '# the fixed central-value truth smeared by each member\'s own A_C. Statistical reference = DataStats covariance of the\n'
            '# member; full = total covariance. Offsets are of the ensemble-mean integral from the reference integral.\n')
    o.write('extraction\tmembers\tbins\tpull_mean_stat\tpull_width_stat\tcov68_stat\tcov95_stat\tpull_mean_full\tpull_width_full\tcov68_full\tcov95_full\tmean_integral\treference_integral\toffset_pct\toffset_sigma_of_mean\tthrow_to_throw_sd_pct\n')
    for tag,_ in TAGS:
        if tag not in res: continue
        r=res[tag]; s=r['statistical only']; f=r['full']
        o.write(f"{tag}\t{r['members']}\t{r['bins']}\t{s['mean']:.2f}\t{s['width']:.2f}\t{100*s['cov68']:.1f}\t{100*s['cov95']:.1f}\t{f['mean']:.2f}\t{f['width']:.2f}\t{100*f['cov68']:.1f}\t{100*f['cov95']:.1f}\t{r['mean_integral']:.4f}\t{r['reference_integral']:.4f}\t{r['offset_pct']:.2f}\t{r['offset_sigma']:.1f}\t{r['sd_pct']:.1f}\n")
print('wrote data_release/ensemble_2026-10-05.tsv')
def row(tag,lab,key):
    r=res[tag][key]; return f"{lab:<40} & ${r['mean']:+.2f}$ & ${r['width']:.2f}$ & ${100*r['cov68']:.1f}\\%$ / ${100*r['cov95']:.1f}\\%$ \\\\"
def patch(path, header_line, blocks):
    t=open(path).read(); i=t.index(header_line); a=t.rindex('\\begin{tabular}',0,i); b=t.index('\\end{tabular}',i)+len('\\end{tabular}')
    body='\\begin{tabular}{lccc}\n\\toprule\n'+header_line+'\n\\midrule\n'+'\n'.join(blocks)+'\n\\midrule\nexpected                & $0$     & $1$    & $68\\%$ / $95\\%$ \\\\\n\\bottomrule\n\\end{tabular}'
    open(path,'w').write(t[:a]+body+t[b:]); print('patched',os.path.basename(path))
stat=[row(t,l,'statistical only') for t,l in TAGS if t in res]; full=[row(t,l,'full') for t,l in TAGS if t in res]
patch(R+'notes/analysis_note.tex','Reference: statistical covariance & pull mean & pull width & 68\\% / 95\\% coverage \\\\', stat)
patch(R+'notes/technical_supplement.tex','Reference covariance & pull mean & pull width & 68\\% / 95\\% coverage \\\\', ['\\multicolumn{4}{l}{\\emph{statistical only}}\\\\']+stat+['\\midrule','\\multicolumn{4}{l}{\\emph{full}}\\\\']+full)
for tag,_ in TAGS:
    if tag in res: r=res[tag]; print(f"{tag:16s} n={r['members']:3d} stat width {r['statistical only']['width']:.2f} mean {r['statistical only']['mean']:+.2f} cov68 {100*r['statistical only']['cov68']:.1f}%  offset {r['offset_pct']:+.2f}% ({r['offset_sigma']:+.1f} sigma)")
