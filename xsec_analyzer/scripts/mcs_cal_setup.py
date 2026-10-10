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


ALTMAN = 'slurm/braggfix_altmodel_manifest.list'
MUMOM = ('pmu', 'Whad', 'dpt', 'dalphat', 'pn')
# release ensembles of the calibrated results: member name in /data/uboone/processed/ens, calibrated tag
ENS = [('fhc5_pmu', 'incl_FHC5_pmu'), ('1p_fhc5_Whad', '1p_FHC5_Whad'), ('1p_fhc5_dpt2bin', '1p_FHC5_dpt2bin'),
       ('1p_comb_pn2bin', '1p_COMB_pn2bin'), ('1p_comb_thetap_dpt_2d', '1p_COMB_thetap_dpt2d')]


def validation():
    """version 1.9 validation inputs of the calibrated results:
      slurm/mcscal_valid_files.list    in|out|map: the alternative-model pseudo-data and the ensemble throws, nominal map
      configs/mcscal/file_properties_numi_<cfg>[_w]_alt<x>.txt, configs/mcscal/alt/<config>
                                       the alternative-model extractions on the calibrated simulation
      slurm/mcscal_alt_manifest.list   "alt config, bin, slice, calibrated tag, version-1.8 output" (slurm/slurm_mcscal_alt.sbatch)
      slurm/mcscal_ens_manifest.list   "member name, calibrated tag, release config, bin, slice" (slurm/slurm_mcscal_ens.sbatch)"""
    import glob
    os.makedirs(os.path.join(REPO, 'configs', 'mcscal', 'alt'), exist_ok=True)
    prod = {p[0]: p for p in PRODUCTS}
    files, alt = {}, []
    for l in open(os.path.join(REPO, ALTMAN)):
        if not l.strip(): continue
        xc, bc, sc = l.split()[:3]
        if not any(m in os.path.basename(xc) for m in MUMOM): continue
        txt = open(os.path.join(REPO, xc)).read()
        U = re.search(r'^UnivFile\s+(\S+)', txt, re.M).group(1); fp = re.search(r'^FPFile\s+(\S+)', txt, re.M).group(1)
        stem = os.path.basename(U)[:-len('_univmake.root')]                       # e.g. ccpi1p_FHC5_Whad_altgenie
        pfx, cfg, rest = stem.split('_', 2)
        obs = rest.rsplit('_', 1)[0]
        tag = f'{"incl" if pfx == "ccpi" else "1p"}_{cfg}_{"thetap_dpt2d" if obs == "thetap_dpt" else obs}'
        assert tag in prod, tag
        out = []
        for row in open(os.path.join(REPO, fp)):
            f = row.split()
            if not row.lstrip().startswith('#') and len(f) >= 3 and f[2] != 'extBNB':
                if f[2] == 'onBNB':
                    new = CAL + 'alt/' + ('w/' if f[0].startswith(P + 'w/') else '') + os.path.basename(f[0])
                    files[f[0]] = f'{f[0]}|{new}|{MAPS}nominal.txt'
                else:
                    new = cal_path(f[0])
                row = row.replace(f[0], new, 1)
            out.append(row)
        open(os.path.join(REPO, 'configs', 'mcscal', os.path.basename(fp)), 'w').write(''.join(out))
        txt = re.sub(r'^UnivFile\s+\S+', f'UnivFile {CAL}alt/univ/{os.path.basename(U)}', txt, flags=re.M)
        txt = re.sub(r'^FPFile\s+\S+', f'FPFile configs/mcscal/{os.path.basename(fp)}', txt, flags=re.M)
        open(os.path.join(REPO, 'configs', 'mcscal', 'alt', os.path.basename(xc)), 'w').write(
            '# MCS-calibrated alternative-model extraction (scripts/mcs_cal_setup.py validation)\n' + txt)
        alt.append(f'configs/mcscal/alt/{os.path.basename(xc)} {bc} {sc} {tag} {os.path.dirname(U)}/xsec_{stem}.root')
    lines = sorted(files.values())
    for d, sub in ((P + 'ens/throws/', 'ens/throws/'), (P + 'w/ens/throws/', 'ens/w/throws/')):
        for f in sorted(glob.glob(d + 'fakedata_*_t*.root')):
            if sub == 'ens/throws/' and '_rhc_' in f: continue                    # the inclusive ensemble is FHC only
            lines.append(f'{f}|{CAL}{sub}{os.path.basename(f)}|{MAPS}nominal.txt')
    open(os.path.join(REPO, 'slurm', 'mcscal_valid_files.list'), 'w').write('\n'.join(lines) + '\n')
    open(os.path.join(REPO, 'slurm', 'mcscal_alt_manifest.list'), 'w').write('\n'.join(alt) + '\n')
    ens = [f'{n} {t} {prod[t][1]} {prod[t][3]} {prod[t][2]}' for n, t in ENS]
    open(os.path.join(REPO, 'slurm', 'mcscal_ens_manifest.list'), 'w').write('\n'.join(ens) + '\n')
    print(f'{len(files)} alternative-model pseudo-data files, {len(lines) - len(files)} ensemble throws, '
          f'{len(alt)} alternative-model extractions, {len(ens)} ensembles')


if __name__ == '__main__':
    import sys
    validation() if sys.argv[1:] == ['validation'] else main()
