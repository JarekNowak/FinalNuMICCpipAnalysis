# extract_eval.py -- A_C row sums, A_C conditioning s_min/s_1, data-statistical term, per-bin total
# uncertainty and closure for every extraction in a slurm_extract.sbatch manifest. For multi-slice
# (2D) configs the closure/A_C numbers come from the full "bin number" slice
# (closure_hists_all_*), otherwise from slice 0.
#   python3 report/tools/extract_eval.py xsec_analyzer/slurm/newobs_manifest.list
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..', 'xsec_analyzer'))
sys.path.insert(0, HERE)
from binning_eval_c50 import dump, side
man = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, 'slurm', 'newobs_manifest.list')
print(f"{'extraction':34s} {'nb':>2s} {'A_C rows':>11s} {'smin/s1':>7s} {'DataStat':>8s} {'bin unc %':>11s} {'closure':>11s} {'chi2/ndf p':>14s}")
for line in open(man):
    xc = line.split()[0]
    U = next(l.split()[1] for l in open(os.path.join(REPO, xc)) if l.startswith('UnivFile'))
    D, B = os.path.dirname(U), os.path.basename(U).replace('_univmake.root', '')
    log = f'{D}/unfold_{B}.log'
    side_all, side0 = f'{D}/closure_hists_all_xsec_{B}.root', f'{D}/closure_hists_xsec_{B}.root'
    sidef = side_all if os.path.exists(side_all) else side0
    if not (os.path.exists(log) and os.path.exists(sidef)):
        print(f'{B:34s} pending'); continue
    try:
        d, chi = dump(log); nb, rows, cond, rel, clo = side(sidef)
    except Exception as e:
        print(f'{B:34s} unreadable ({e})'); continue
    cl = f"{min(clo):.2f}-{max(clo):.2f}" if clo else 'n/a'
    ch = f"{chi[0]:.2f}/{chi[1]} {chi[2]:.2f}" if chi else 'n/a'
    print(f"{B:34s} {nb:2d} {min(rows):.2f}-{max(rows):.2f} {cond:7.3f} {d.get('DataStats', float('nan')):7.1f}% "
          f"{min(rel):4.0f}-{max(rel):4.0f}   {cl:>11s} {ch:>14s}")
