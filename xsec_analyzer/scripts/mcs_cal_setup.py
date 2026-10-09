"""mcs_cal_setup.py -- inputs of the MCS-calibrated build (scripts/mcs_cal_derive.py, macros/mcs_cal_apply.C) of the
released results that use the muon momentum.

  slurm/mcscal_files.list          in|out|map: every simulated file of the inclusive and proton-tagged configurations
                                   (overlay, dirt, detector variations, pseudo-data) with the nominal map, and the
                                   pseudo-data with every variant map (the calibration uncertainty)
  configs/mcscal/file_properties_numi_<cfg>[_w].txt   the release lists with the calibrated files (beam-off unchanged)
  configs/mcscal/<xsec config>     the release extraction configs with these lists, their own universe files and the
                                   version-1.8 systematics (detector variations fractional)
  slurm/mcscal_manifest.list       "xsec bin slice" for slurm/slurm_extract.sbatch

    python3 scripts/mcs_cal_setup.py
"""
import os, re

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
P = '/data/uboone/processed/'
CAL = P + 'mcscal/'
MAPS = CAL + 'maps/'
VARIANTS = ['shrink', 'coarse', 'fine', 'old_production', 'new_production', 'FHC', 'RHC']
CFGS = ['fhc5', 'rhcfull', 'comb']
SYST = 'configs/ccpi_systcalc_numi_dvfrac.conf'          # identical to version 1.8 (inclusive and proton-tagged)
PRODUCTS = [l.strip().split('|') for l in open(os.path.join(REPO, 'slurm', 'mcsmap_products.list')) if l.strip()]
PRODUCTS += [['1p_COMB_thetap_dpt2d', 'configs/d2/ccpi1p_xsec_config_numi_thetap_dpt_2d_comb.txt',
              'configs/d2/ccpi1p_thetap_dpt_slice_config_2d.txt', 'configs/d2/ccpi1p_thetap_dpt_bin_config_2d.txt', '1p']]


def cal_path(path, variant='nominal'):
    """calibrated location of a processed file (inclusive: processed/, proton-tagged: processed/w/)"""
    if path.startswith(P + 'w/'): return CAL + variant + '/w/' + os.path.basename(path)
    if os.path.dirname(path) + '/' == P: return CAL + variant + '/' + os.path.basename(path)
    raise ValueError(path)


def rows(fp):
    for l in open(os.path.join(REPO, fp)):
        f = l.split()
        if l.lstrip().startswith('#') or len(f) < 3: continue
        yield l, f


def main():
    os.makedirs(os.path.join(REPO, 'configs', 'mcscal'), exist_ok=True)
    sim, fake = set(), set()
    for cfg in CFGS:
        for w in ('', '_w'):
            fp = f'configs/file_properties_numi_{cfg}{w}.txt'
            out = []
            for l in open(os.path.join(REPO, fp)):
                f = l.split()
                if not l.lstrip().startswith('#') and len(f) >= 3 and f[2] != 'extBNB':
                    sim.add(f[0])
                    if f[2] == 'onBNB': fake.add(f[0])
                    l = l.replace(f[0], cal_path(f[0]), 1)
                out.append(l)
            open(os.path.join(REPO, 'configs', 'mcscal', os.path.basename(fp)), 'w').write(''.join(out))
    lines = [f'{s}|{cal_path(s)}|{MAPS}nominal.txt' for s in sorted(sim)]
    lines += [f'{s}|{cal_path(s, v)}|{MAPS}{v}.txt' for v in VARIANTS for s in sorted(fake)]
    open(os.path.join(REPO, 'slurm', 'mcscal_files.list'), 'w').write('\n'.join(lines) + '\n')
    man = []
    for tag, xc, sc, bc, fam in PRODUCTS:
        txt = open(os.path.join(REPO, xc)).read()
        fp = re.search(r'^FPFile\s+(\S+)', txt, re.M).group(1)
        txt = re.sub(r'^UnivFile\s+\S+', f'UnivFile {CAL}univ/{tag}_univmake.root', txt, flags=re.M)
        txt = re.sub(r'^SystFile\s+\S+', f'SystFile {SYST}', txt, flags=re.M)
        txt = re.sub(r'^FPFile\s+\S+', f'FPFile configs/mcscal/{os.path.basename(fp)}', txt, flags=re.M)
        txt = '# MCS-calibrated build of ' + xc + ' (scripts/mcs_cal_setup.py)\n' + txt
        out = f'configs/mcscal/{tag}.txt'
        open(os.path.join(REPO, out), 'w').write(txt)
        man.append(f'{out} {bc} {sc}')
    open(os.path.join(REPO, 'slurm', 'mcscal_manifest.list'), 'w').write('\n'.join(man) + '\n')
    print(f'{len(sim)} simulated files ({len(fake)} pseudo-data), {len(lines)} calibration jobs, {len(man)} products')


if __name__ == '__main__':
    main()
