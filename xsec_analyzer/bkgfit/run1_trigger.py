"""run1_trigger.py -- is the FHC Run-1 beam-on sample the open-trigger dataset or not?

RESOLVED: the beam bookkeeping confirms 2.192e20 POT and 5 748 692 triggers (no open trigger),
which is the pairing this test selected. The configs carry it since commit 619e5e8; this script
remains as the record of the evidence.

Two bookkeeping values exist for FHC Run 1:

    without open trigger   2.192e20 POT,  5 748 692 triggers
    with    open trigger   3.283e20 POT,  9 846 635 triggers

Every configuration in this analysis assumes the second. The beam-on skim carries no POT of its own
(summed_pot is zero in the skim and the source ntuple's SubRun tree holds only run/subRun), so the
question cannot be answered by reading a file -- but it can be answered by the data, because the
assumption enters the prediction three ways at once and all three scale with it:

    neutrino MC   scales with the assumed data POT,
    beam-off      scales with the assumed beam-on trigger (gate) count,
    dirt          scales with the assumed data POT.

If the wrong pair is assumed, every control region in that period is off by the same factor, which
is exactly the signature reported for FHC Run 1. This recomputes the four regions under both
hypotheses. The one that puts the ratios at unity is the exposure the data actually has.

    python3 -m bkgfit.run1_trigger
"""
import os
import numpy as np

from .farsb import read, klass, CLASSES, P, BO, G1, OCCX, DIRT

MC_R1 = P + 'xsec-ana-Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root'
EXT_R1 = P + 'xsec-ana-neutrinoselection_filt_run1_beamoff.root'
DATA_R1 = BO + 'xsec-ana-beamon_fhc_run1.root'

# the two bookkeeping pairs, and the scales the analysis derives from them
# The POT and the trigger count enter different parts of the prediction, so the two mixed
# hypotheses are tested as well: they separate "the exposure is wrong" from "the gate count is
# wrong", which a single combined test cannot do.
HYP = {
    'assumed: open trigger POT + open trigger gates': dict(pot=3.283, trig=9846635., mc=0.14101),
    'neither: 2.192e20 POT + 5.75M gates':            dict(pot=2.192, trig=5748692., mc=0.14101 * 2.192 / 3.283),
    'mixed: open POT + non-open gates':               dict(pot=3.283, trig=5748692., mc=0.14101),
    'mixed: non-open POT + open gates':               dict(pot=2.192, trig=9846635., mc=0.14101 * 2.192 / 3.283),
    # A third pairing exists in the repository for THIS beam-on file: configs/file_properties_numi_3283.txt
    # records 3.283e20 POT with 7809962 triggers for neutrinoselection_filt_run1_beamon_beamgood.root,
    # and run_dv_fullrebuild.sh derives its trigger counts from the same 7809962 / 3.283e20. The live
    # 9846635 appears only in the scaling code, never next to a POT for this file.
    'repo pairing: 3.283e20 + 7.81M gates':           dict(pot=3.283, trig=7809962., mc=0.14101),
    'non-open POT + 7.81M gates':                     dict(pot=2.192, trig=7809962., mc=0.14101 * 2.192 / 3.283),
}
REGIONS = [('CC0pi', 'cc0pi'), ('pi0', 'pi0'), ('multi-pi', 'multipi'), ('cosmic', 'cosmic')]


def main():
    for f in (MC_R1, EXT_R1, DATA_R1):
        if not os.path.exists(f):
            print('MISSING', f); return
    mc, ext, dat = read(MC_R1, True), read(EXT_R1, False), read(DATA_R1, False)
    dirt = read(DIRT, True)

    print('%-24s %8s %9s %9s %9s %8s' % ('', 'region', 'nu MC', 'beam-off', 'dirt', 'data'))
    for name, h in HYP.items():
        print(f'\n--- {name}: {h["pot"]}e20 POT, {h["trig"]:.0f} triggers')
        sd = 0.092402 * 0.65 * h['pot'] / 8.857          # dirt, as bkgfit_templates.C scales it
        se = OCCX * h['trig'] / G1                        # beam-off, by gate ratio
        for lab, key in REGIONS:
            m = mc[key]; e = ext[key]; d = dat[key]; dt = dirt[key]
            nu = (mc['w'][m] * h['mc']).sum()
            eo = e.sum() * se
            di = (dirt['w'][dt] * sd).sum()
            pred = nu + eo + di
            obs = d.sum()
            print('%-24s %8s %9.1f %9.1f %9.1f %8d   data/pred = %.3f'
                  % ('', lab, nu, eo, di, obs, obs / pred if pred else 0))


if __name__ == '__main__':
    main()
