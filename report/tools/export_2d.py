"""export_2d.py -- release the two-dimensional results proposed for approval (2026-09-27): the curves
(unfolded result and every A_C-smeared model, per analysis bin) and the additional-smearing matrix, in
the same conventions as export_curves.C / export_matrices.C, appended to index_curves.tsv and
index_A_C.tsv. Run AFTER the two ROOT exporters, which rewrite the indices.

The flattened measurement is ordered outer slice by outer slice (inner variable fastest); each entry of
the closure file is the cross section integrated over the analysis bin, so d2sigma/dXdY = value/area.
    python3 report/tools/export_2d.py
"""
import os, sys
import numpy as np, uproot
HERE = os.path.dirname(os.path.abspath(__file__)); REL = os.path.join(os.path.dirname(HERE), 'data_release')
RB = '/data/uboone/processed/rebuild_2d'
# tag, closure file, inner (name, edges), outer (name, edges), open-top note
PAIRS = [('1p_COMB_thetap_dpt2d', 'closure_hists_all_xsec_ccpi1p_COMB_thetap_dpt.root',
          ('theta_p [rad]', [0, 0.8, np.pi]), ('delta_pT [GeV/c]', [0, 0.3, 2.5]),
          'the upper delta p_T edge (2.5 GeV/c) is a plotting boundary: the truth bin is open'),
         ('incl_COMB_costhpi_costhmu2d', 'closure_hists_all_xsec_ccpi_COMB_costhpi_costhmu.root',
          ('cos_theta_pi', [-1, 0.35, 1]), ('cos_theta_mu', [-1, 0.65, 0.85, 1]), '')]
KEYS = [('h_unfolded_nuwro', 'unfolded_data'), ('h_fakedata_truth', 'truth_smeared'), ('h_genie_tune', 'tune_smeared'),
        ('h_gen_GENIE', 'GENIE_smeared'), ('h_gen_GiBUU', 'GiBUU_smeared'), ('h_gen_NEUT', 'NEUT_smeared'), ('h_gen_NuWro', 'NuWro_smeared')]
for tag, fn, (xin, xe), (yout, ye), note in PAIRS:
    f = uproot.open(os.path.join(RB, fn)); nx, ny = len(xe) - 1, len(ye) - 1
    vals = {lab: f[k].values() for k, lab in KEYS}; err = f['h_unfolded_nuwro'].errors()
    assert len(vals['unfolded_data']) == nx * ny, (tag, len(vals['unfolded_data']))
    with open(os.path.join(REL, f'curves_{tag}.tsv'), 'w') as o:
        o.write(f'# Measurement-space curves for {tag} (two-dimensional, combined configuration).\n'
                '# Analysis bins ordered outer slice by outer slice, inner variable fastest. All model columns are\n'
                '# ALREADY smeared by A_C. Values are d2sigma/dXdY in 1e-38 cm^2/(unit X)/(unit Y)/Ar; multiply by\n'
                f'# area to integrate. Inner X = {xin}, outer Y = {yout}.' + (f' Note: {note}.' if note else '') + '\n')
        o.write('bin\tinner_low\tinner_high\touter_low\touter_high\tarea\tunfolded_data\tstat_err\t' + '\t'.join(l for _, l in KEYS[1:]) + '\n')
        for k in range(nx * ny):
            j, i = divmod(k, nx); a = (xe[i + 1] - xe[i]) * (ye[j + 1] - ye[j])
            o.write(f'{k}\t{xe[i]:.6g}\t{xe[i+1]:.6g}\t{ye[j]:.6g}\t{ye[j+1]:.6g}\t{a:.6g}\t{vals["unfolded_data"][k]/a:.8g}\t{err[k]/a:.8g}\t'
                    + '\t'.join(f'{vals[l][k]/a:.8g}' for _, l in KEYS[1:]) + '\n')
    A = f['h_A_C'].values()          # uproot: [x (true), y (smeared)]
    n = A.shape[0]
    with open(os.path.join(REL, f'A_C_{tag}.tsv'), 'w') as o:
        o.write(f'# Additional smearing matrix A_C for {tag}\n# Compare a prediction p to the published result as (A_C . p), summing over\n'
                '# the TRUE index. Row = smeared bin i, column = true bin j, value = A_C[i][j].\n'
                f'# {n} smeared bins x {n} true bins; analysis-bin order as in curves_{tag}.tsv.\nsmeared_bin' + ''.join(f'\ttrue_{j+1}' for j in range(n)) + '\n')
        rs = []
        for i in range(n):
            o.write(f'{i+1}' + ''.join(f'\t{A[j, i]:.10g}' for j in range(n)) + '\n'); rs.append(A[:, i].sum())
    for idx, row in [('index_curves.tsv', f'{tag}\t{n}\n'), ('index_A_C.tsv', f'{tag}\t{n}\t{n}\t{min(rs):.6g}\t{max(rs):.6g}\t{fn}\n')]:
        p = os.path.join(REL, idx); lines = [l for l in open(p) if not l.startswith(tag + '\t')]
        open(p, 'w').write(''.join(lines) + row)
    print(f'wrote curves_{tag}.tsv and A_C_{tag}.tsv ({n} bins; A_C row sums {min(rs):.3f}-{max(rs):.3f})')
