"""ensemble_tables_v1.py -- evaluate the pseudo-data ensembles of prerequisite V1 (first run 2026-09-28, re-run for version 1.7 on
2026-10-05): the proton angles
theta_p and theta_pi_p in FHC, RHC and combined, and the two combined-only two-dimensional results,
theta_p x delta p_T (proton-tagged) and cos theta_pi x cos theta_mu (inclusive); 100 members each
(xsec_analyzer/ens_braggfix.sh). Writes data_release/ensemble_v1_2026-10-05.tsv with the per-bin statistical pull
means and widths as well, since a two-dimensional cell can behave differently from the extraction as a whole.
The coverage table of the proton-tagged note (tables/ensemble_1p.tex) reads this file through
ensemble_tables_1p.py; the inclusive two-dimensional row is quoted in the analysis note. Run
xsec_analyzer/slurm/ens_statcov.sh <cfg> <obs> first, so every member has its DataStats table, and
ensemble_entry_scan.py on every tag.

    python3 report/tools/ensemble_tables_v1.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ensemble_stat_pulls import evaluate
R = os.path.abspath(os.path.join(HERE, '..')) + '/'
OUT = R + 'data_release/ensemble_v1_2026-10-05.tsv'
TAGS = ['1p_fhc5_thetap', '1p_rhcfull_thetap', '1p_comb_thetap',
        '1p_fhc5_thpipr', '1p_rhcfull_thpipr', '1p_comb_thpipr',
        '1p_comb_thetap_dpt_2d', 'comb_costhpi_costhmu_2d']


def main():
    res = {}
    for tag in TAGS:
        try: res[tag] = evaluate(tag, 100, quiet=True)
        except Exception as e: print(tag, 'failed:', e)
    with open(OUT, 'w') as o:
        o.write('# Pseudo-data ensembles of prerequisite V1: Poisson throws of the central-value prediction, run by run, each\n'
                '# carried through univmake + UnfolderNuMI (slurm_ensemble_1p.sbatch, slurm_ensemble_cfg.sbatch; ens_braggfix.sh); pulls\n'
                '# against the fixed central-value truth smeared by each member\'s own A_C. Statistical reference = DataStats\n'
                '# covariance of the member; full = total covariance. Offsets are of the ensemble-mean integral from the reference\n'
                '# integral. Two-dimensional bins are ordered outer slice by outer slice, inner variable fastest.\n')
        o.write('extraction\tmembers\tbins\tpull_mean_stat\tpull_width_stat\tcov68_stat\tcov95_stat\tpull_mean_full\tpull_width_full'
                '\tcov68_full\tcov95_full\tmean_integral\treference_integral\toffset_pct\toffset_sigma_of_mean\tthrow_to_throw_sd_pct'
                '\tbin_pull_means_stat\tbin_pull_widths_stat\n')
        for tag in TAGS:
            if tag not in res: continue
            r = res[tag]; s = r['statistical only']; f = r['full']; pb = r['statistical only per-bin']
            o.write(f"{tag}\t{r['members']}\t{r['bins']}\t{s['mean']:.2f}\t{s['width']:.2f}\t{100*s['cov68']:.1f}\t{100*s['cov95']:.1f}"
                    f"\t{f['mean']:.2f}\t{f['width']:.2f}\t{100*f['cov68']:.1f}\t{100*f['cov95']:.1f}\t{r['mean_integral']:.4f}"
                    f"\t{r['reference_integral']:.4f}\t{r['offset_pct']:.2f}\t{r['offset_sigma']:.1f}\t{r['sd_pct']:.1f}"
                    f"\t{','.join(f'{b[0]:.2f}' for b in pb)}\t{','.join(f'{b[1]:.2f}' for b in pb)}\n")
    print('wrote', os.path.relpath(OUT, R))
    for tag in TAGS:
        if tag not in res: continue
        r = res[tag]; s = r['statistical only']; f = r['full']; pb = r['statistical only per-bin']
        print(f"{tag:24s} n={r['members']:3d} bins={r['bins']} stat mean {s['mean']:+.2f} width {s['width']:.2f} "
              f"cov {100*s['cov68']:.1f}/{100*s['cov95']:.1f}  full width {f['width']:.2f}  offset {r['offset_pct']:+.2f}% "
              f"({r['offset_sigma']:+.1f} sigma)  bins: " + ' '.join(f'{b[0]:+.2f}/{b[1]:.2f}' for b in pb))


if __name__ == '__main__':
    main()
