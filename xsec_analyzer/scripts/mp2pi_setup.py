#!/usr/bin/env python3
"""mp2pi_setup.py -- framework inputs of the two-pion one-bin total (Phase 4, item 3, of
report/planning/MULTIPION_BNB_ADAPTATION_PLAN.md), from the production in /data/uboone/processed/mp2pi
(slurm/slurm_mp2pi.sbatch, CC1mu2piBDT and CC1mu1piXp).

1. File lists configs/file_properties_numi_{fhc5,rhcfull,comb}_mp2pi.txt: the release lists
   (configs/file_properties_numi_*_w.txt) with every path moved to the mp2pi directory under the same file name;
   per-run names that the release resolves through links (beam-off and dirt shared by several run periods) are made
   links to the reprocessed sample there.
2. Pseudo-data plan /data/uboone/processed/mp2pi/throw_plan.txt for macros/throw_mp2pi.C: per beam-on entry of the
   release lists, the overlays AND the dirt of that run period, each with the scale data POT / sum of the distinct
   summed_pot of its type (the framework's normalisation), so the pseudo-data carry dirt as the prediction does;
   the framework adds the beam-off to fake data itself.

    python3 scripts/mp2pi_setup.py
"""
import os, sys
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mp_phase0_baseline as B

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = '/data/uboone/processed/mp2pi/'
CFGS = {'fhc5': 'file_properties_numi_fhc5_w.txt', 'rhcfull': 'file_properties_numi_rhcfull_w.txt',
        'comb': 'file_properties_numi_comb_w.txt'}


def target(path, ftype):
    """File in the mp2pi production that a release entry stands for."""
    real = os.path.basename(os.path.realpath(path))
    if ftype == 'dirtMC':
        return B.DIRT_ALIAS.get(real, real)
    if ftype == 'onBNB':
        return None                      # thrown by macros/throw_mp2pi.C
    return real if ftype == 'extBNB' else os.path.basename(path)


def summed_pot(f):
    return float(uproot.open(f)['summed_pot'].member('fVal'))


def main():
    links, plan = {}, {}
    for cfg, name in CFGS.items():
        lines_out = []
        per_run = {}
        for raw in open(os.path.join(REPO, 'configs', name)):
            line = raw.rstrip('\n')
            body = line.split('#')[0].strip()
            if not body:
                lines_out.append(line); continue
            f = body.split()
            path, run, ftype = f[0], int(f[1]), f[2]
            base = os.path.basename(path)
            new = OUT + base
            tgt = target(path, ftype)
            if tgt is not None and tgt != base:
                links[base] = tgt
            if tgt is not None and not os.path.exists(OUT + tgt):
                sys.exit(f'missing production output {OUT + tgt} for {path}')
            lines_out.append(' '.join([new] + f[1:]))
            per_run.setdefault(run, {}).setdefault(ftype, []).append((new, f))
        open(os.path.join(REPO, 'configs', f'file_properties_numi_{cfg}_mp2pi.txt'), 'w').write(
            f'# two-pion one-bin total ({cfg}): {name} with the files of /data/uboone/processed/mp2pi (scripts/mp2pi_setup.py)\n'
            + '\n'.join(lines_out) + '\n')
        for run, d in per_run.items():
            for new, f in d.get('onBNB', []):
                out = os.path.basename(new)
                if out in plan: continue
                dpot = float(f[4])
                ins = []
                for typ in ('numuMC', 'dirtMC'):
                    files = [OUT + target(f_[0], typ) for _, f_ in d.get(typ, [])]
                    distinct = []
                    for fl in files:
                        sp = summed_pot(fl)
                        if not any(abs(sp - s) <= 1e-6 * abs(s) for s in distinct): distinct.append(sp)
                    for fl in files: ins.append((fl, dpot / sum(distinct)))
                plan[out] = (dpot, ins)
    for base, tgt in links.items():
        lp = OUT + base
        if os.path.islink(lp) or os.path.exists(lp):
            if os.path.realpath(lp) != os.path.realpath(OUT + tgt): sys.exit(f'{lp} exists and is not the expected link')
            continue
        os.symlink(tgt, lp)
    with open(OUT + 'throw_plan.txt', 'w') as fp:
        for k, (out, (dpot, ins)) in enumerate(sorted(plan.items())):
            fp.write(f'OUT {out} {dpot:.6e} {k + 1}\n')
            for fl, sc in ins: fp.write(f'IN {fl} {sc:.8e}\n')
    print(f'{len(links)} links, {len(plan)} pseudo-data files planned:')
    for out, (dpot, ins) in sorted(plan.items()):
        print(f'  {out}: data POT {dpot:.4g}, ' + ', '.join(f'{os.path.basename(f)} x {s:.4g}' for f, s in ins))


if __name__ == '__main__':
    main()
