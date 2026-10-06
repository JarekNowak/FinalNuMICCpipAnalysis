"""make_braggfix_univ_manifest.py -- the 65 universe files of the release, each with the EXACT bin
configuration it was built from, for the rebuild after the Bragg-pion cut was dropped (commit 1fb6be9).

A univmake file stores, per directory, the reco and true bin specifications, the selection, the ntuple
name and the file-properties list it was built from (directory total_configs+<list>). This writes each
directory back as a bin-configuration file in configs/braggfix/ (checked byte-identical against
configs/ccpi_pmu_bin_config_opt.txt) and one manifest line per universe file:

    <live universe file>|<file-properties list>|<bin config>[;<bin config>...]

    python3 slurm/make_braggfix_univ_manifest.py  ->  slurm/braggfix_univ_manifest.list, configs/braggfix/*.txt
"""
import os
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CFG = ('fhc5', 'rhcfull', 'comb')
INCL = ['pmu', 'ppi3bin', 'costhmu', 'costhpi', 'thmupi', 'total']
ONEP = ['pmu', 'ppi2bin', 'costhmu', 'costhpi', 'thmupi', 'Wpipr', 'Whad', 'dpt2bin', 'dphit3bin', 'dalphat2bin',
        'pn2bin', 'pp', 'total', 'thetap', 'thpipr']
D2 = ['d2/ccpi_xsec_config_numi_costhpi_costhmu_2d_comb.txt', 'd2/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb.txt']


def univfile(xsec_config):
    for line in open(os.path.join(REPO, 'configs', xsec_config)):
        tok = line.split()
        if tok and tok[0] == 'UnivFile': return tok[1]
    raise SystemExit(f'no UnivFile in {xsec_config}')


def main():
    xcs = [f'ccpi_xsec_config_numi_{o}_{c}.txt' for o in INCL for c in CFG] + \
          [f'ccpi1p_xsec_config_numi_{o}_{c}.txt' for o in ONEP for c in CFG] + D2
    out = os.path.join(REPO, 'configs', 'braggfix'); os.makedirs(out, exist_ok=True)
    lines, seen = [], set()
    for xc in xcs:
        u = univfile(xc)
        if u in seen: continue
        seen.add(u)
        f = uproot.open(u)
        bins, fpm = [], set()
        for name in sorted({k.split(';')[0] for k in f.keys(recursive=False)}):
            d = f[name]
            sub = {k.split(';')[0] for k in d.keys(recursive=False)}
            if 'reco_bin_spec' not in sub: continue            # unfolder caches, not a univmake directory
            t = d['true_bin_spec'].strip('\n').split('\n'); r = d['reco_bin_spec'].strip('\n').split('\n')
            txt = '\n'.join([name, d['ntuple_name'], d['sel_for_categ'], str(len(t))] + t + [str(len(r))] + r) + '\n'
            b = os.path.join(out, f'{os.path.basename(u)[:-len("_univmake.root")]}__{name}.txt')
            open(b, 'w').write(txt); bins.append(os.path.relpath(b, REPO))
            fpm |= {s[len('total_configs+'):] for s in sub if s.startswith('total_configs+')}
        if len(fpm) != 1 or not bins: raise SystemExit(f'{u}: file-properties lists {fpm}, {len(bins)} directories')
        lines.append(f'{u}|configs/{fpm.pop()}|{";".join(bins)}')
    open(os.path.join(HERE, 'braggfix_univ_manifest.list'), 'w').write('\n'.join(lines) + '\n')
    print(len(lines), 'universe files')


if __name__ == '__main__':
    main()
