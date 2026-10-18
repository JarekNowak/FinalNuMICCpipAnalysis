"""mp2pi_kin_configs.py -- framework configurations of the two-pion differential measurement, one variable at a time,
from the binnings of scripts/mp2pi_kinematics.py (report/planning/multipion/mp2pi_kinematics.json: edges per variable,
chosen under the 0.50 criterion with at least 40 predicted signal events per bin).

Per variable v:
  configs/mp2pi_kin_<v>_bin_config.txt    true bins of the signal in the truth variable (the first open below, the last
                                          open above, so that every signal event has a bin) and one background bin;
                                          reco bins of the signal region in the reco variable (defined values only)
                                          and the four sidebands of the one-bin total as sideband bins of the signal
                                          block, so that the constraint works as for the total
  configs/mp2pi_kin_<v>_slice_config.txt  one slice in the variable
  configs/mp2pi_xsec_config_numi_kin_<v>_<cfg>[_constr].txt   Wiener-SVD (filter on, second derivative) as the
                                          single-pion differential results; detector variations fractional; the
                                          universes of slurm/slurm_univmake_mp2pi.sbatch with MP2PI_SUFFIX=_kin_fluxw
                                          and MP2PI_OUTTAG=_<v>

    python3 scripts/mp2pi_kin_configs.py <mp2pi_kinematics.json> [variables] [--prod 2f]
        (--prod P: extraction configs mp2pi_xsec_config_numi_kinP_<v>_<cfg>[_constr].txt on the production of
        /data/uboone/processed/mp2pi_P, universes mp2pi_<CFG>_total_P_fluxw_<v>_univmake.root; v = total uses the bin
        and slice configs of the one-bin total, copied to mp2pi_kin_total_{bin,slice}_config.txt, with Wiener-SVD off)
"""
import json, os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
E = 'CC1mu2piBDT_'
# variable: reco branch, truth branch, axis title, units, upper display edge of an open last bin
VAR = {
    'pmu':    ('kin_pmu', 'kin_true_pmu', 'p_{#mu}', ' (GeV/c)', 2.5),
    'cmu':    ('kin_cmu', 'kin_true_cmu', 'cos#theta_{#mu}', '', 1.0),
    'cfwd':   ('kin_cpi_fwd', 'kin_true_cpi_fwd', 'cos#theta_{#pi}^{fwd}', '', 1.0),
    'cbwd':   ('kin_cpi_bwd', 'kin_true_cpi_bwd', 'cos#theta_{#pi}^{bwd}', '', 1.0),
    'cpipi':  ('kin_cpipi', 'kin_true_cpipi', 'cos#theta_{#pi#pi}', '', 1.0),
    'cmupi':  ('kin_cmupi', 'kin_true_cmupi', 'cos#theta_{#mu#pi}', '', 1.0),
    'ppi_lo': ('kin_ppi_lo', 'kin_true_ppi_lo', 'p_{#pi}^{low}', ' (GeV/c)', 0.6),
    'ppi_hi': ('kin_ppi_hi', 'kin_true_ppi_hi', 'p_{#pi}^{high}', ' (GeV/c)', 0.8),
    'np':     ('kin_np', 'mc_n_proton_sub', 'N_{p}', '', 1.5),
}
SIDEBANDS = [f'1 0 "{E}evt_region == {r} && !CC1mu1piXp_Selected"' for r in (1, 2, 3, 4)]
TAG = {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}
FLUX = {'fhc5': '6.81159e-10', 'rhcfull': '6.44646e-10', 'comb': '6.596906e-10'}


def cut(br, lo, hi, first, last, reco):
    c = []
    if reco: c.append(f'{E}{br} > -998')
    if not first: c.append(f'{E}{br} >= {lo:g}')
    if not last: c.append(f'{E}{br} < {hi:g}')
    return ' && '.join(c)


PROD = sys.argv[sys.argv.index('--prod') + 1] if '--prod' in sys.argv else ''


def xsec_configs(v, n):
    pk = PROD or 'kin'
    for cfg in TAG:
        for mode in ('', '_constr'):
            lines = [f'# Two-pion {"one-bin total" if v == "total" else "differential measurement, variable " + v} ({cfg}), '
                     f'{n} bins; production {pk} (MCS-calibrated classifier inputs' + (', two-fold event classifier' if PROD == '2f' else '')
                     + '), detector variations fractional and weighted to the current flux' + ('; with the sideband constraint' if mode else '') + '.',
                     f'UnivFile /data/uboone/processed/mp2pi_{TAG[cfg]}_total_{pk}_fluxw_{v}_univmake.root',
                     'SystFile configs/ccpi_systcalc_numi_dvfrac.conf',
                     f'FPFile configs/file_properties_numi_{cfg}_mp2pi_{pk}_fluxw.txt',
                     'Unfold WienerSVD 0' if v == 'total' else 'Unfold WienerSVD 1 second-deriv', f'Flux {FLUX[cfg]}'] + \
                    (['Calculator Constrained'] if mode else []) + \
                    ['Prediction uBTune "MicroBooNE Tune" univ CV', 'Prediction FakeData "Fakedata" univ FakeData']
            name = f'mp2pi_xsec_config_numi_kin{PROD}_{v}_{cfg}{mode}.txt'
            with open(f'{REPO}/configs/{name}', 'w') as f: f.write('\n'.join(lines) + '\n')


def write(v, edges):
    rb, tb, title, unit, top = VAR[v]
    n = len(edges) - 1
    true = [f'0 0 "{E}MC_Signal' + (' && ' + cut(tb, edges[i], edges[i + 1], i == 0, i == n - 1, False) if n > 1 else '') + '"'
            for i in range(n)] + [f'1 -1 "!{E}MC_Signal"']
    reco = [f'0 0 "{E}evt_region == 0 && ' + cut(rb, edges[i], edges[i + 1], i == 0, i == n - 1, True) + '"'
            for i in range(n)] + SIDEBANDS
    with open(f'{REPO}/configs/mp2pi_kin_{v}_bin_config.txt', 'w') as f:
        f.write('\n'.join([f'mp2pi_CC1mu2piBDT_kin_{v}_1D', 'stv_tree', 'CC1mu2piBDT', str(len(true))] + true
                          + [str(len(reco))] + reco) + '\n')
    disp = list(edges[:-1]) + [min(edges[-1], top) if edges[-1] > top else edges[-1]]
    if disp[-1] <= disp[-2]: disp[-1] = top
    sl = ['2', f'"{title}" "{unit}" "{title}" "{unit}"', '"bin number" "" "bin number" ""', '2',
          '"events"', '1', f'0 {n + 1} ' + ' '.join(f'{x:g}' for x in disp), '0', str(n)] + \
         [f'{i} 1 {i + 1}' for i in range(n)] + \
         ['"events"', '1', f'1 {n + 1} ' + ' '.join(str(i) for i in range(n + 1)), '0', str(n)] + \
         [f'{i} 1 {i + 1}' for i in range(n)]
    with open(f'{REPO}/configs/mp2pi_kin_{v}_slice_config.txt', 'w') as f: f.write('\n'.join(sl) + '\n')
    xsec_configs(v, n)
    print(f'{v}: {n} bins, edges {edges}')


if __name__ == '__main__':
    res = json.load(open(sys.argv[1]))
    args = [a for a in sys.argv[2:] if not a.startswith('--') and a != PROD]
    for v in (args[0].split(',') if args else res):
        if v == 'total': xsec_configs('total', 1); print('total: one bin'); continue
        if PROD: xsec_configs(v, len(res[v]['edges']) - 1); print(f'{v}: extraction configs for production {PROD}'); continue
        write(v, res[v]['edges'])
