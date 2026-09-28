"""check_1pinc_bins.py -- the edges of the proton-tagged muon/pion generator histograms (obs_1pinc.h, written to
gen2d/<gen>_1p_inc_<flavour>.root) against the TRUE bins of the adopted proton-tagged bin configs: interior
edges equal, the first edge equal to the first lower bound, the last bin open exactly where the config's is.
    python3 gen2d/check_1pinc_bins.py      (from generator_predictions/)"""
import re, sys, numpy as np, uproot
CF = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/configs/'
MAP = {'pmu1p': 'ccpi1p_pmu_bin_config.txt', 'ppi1p': 'ccpi1p_ppi_bin_config_2bin.txt',
       'costhmu1p': 'ccpi1p_costhmu_bin_config.txt', 'costhpi1p': 'ccpi1p_costhpi_bin_config.txt',
       'thmupi1p': 'ccpi1p_thmupi_bin_config.txt'}
OPEN = {'pmu1p': True, 'ppi1p': True, 'costhmu1p': False, 'costhpi1p': True, 'thmupi1p': True}   # obs_1pinc.h topOpen

def true_bins(f):
    out = []
    for l in open(CF + f):
        l = l.strip()
        if l.startswith('0 0 "') and '_MC_Signal' in l:
            s = l.split('"')[1]
            lo = re.search(r'>=\s*(-?[0-9.]+)', s); hi = re.search(r'<\s*(-?[0-9.]+)', s)
            out.append((float(lo.group(1)) if lo else None, float(hi.group(1)) if hi else None))
    return out

bad = 0
for g in ('genie', 'gibuu', 'neut', 'nuwro'):
    for fl in ('numu', 'numubar'):
        f = uproot.open(f'gen2d/{g}_1p_inc_{fl}.root')
        for h, cfg in MAP.items():
            e = f[h].axis().edges(); tb = true_bins(cfg); n = len(e) - 1; msg = []
            if n != len(tb): msg.append(f'{n} bins vs {len(tb)}')
            else:
                for i, (lo, hi) in enumerate(tb):
                    if lo is None or abs(lo - e[i]) > 1e-9: msg.append(f'bin {i} low {e[i]} vs {lo}')
                    if i < n - 1 and (hi is None or abs(hi - e[i + 1]) > 1e-9): msg.append(f'bin {i} high {e[i+1]} vs {hi}')
                if (tb[-1][1] is None) != OPEN[h]: msg.append(f'last bin open={tb[-1][1] is None} in config, {OPEN[h]} in obs_1pinc.h')
                if tb[-1][1] is not None and abs(tb[-1][1] - e[-1]) > 1e-9: msg.append(f'last edge {e[-1]} vs {tb[-1][1]}')
            bad += bool(msg)
            if g == 'genie' and fl == 'numu' or msg:
                print(f'  {g:6s} {fl:8s} {h:10s} {list(np.round(e, 4))} vs {cfg}: ' + ('; '.join(msg) if msg else 'OK'))
print('bin check:', 'ALL OK (4 generators x 2 flavours x 5 observables)' if not bad else f'{bad} MISMATCHES')
sys.exit(1 if bad else 0)
