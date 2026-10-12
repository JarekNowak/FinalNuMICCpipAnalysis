#!/usr/bin/env python3
"""detvar_fluxw.py -- per-event new / old NuMI flux weight of one Run-4 detector-variation sample, and its reweighted
copy (report/planning/multipion/PHASE4_SUMMARY.md, "Detector variations").

The detector-variation samples were simulated with the old NuMI flux, the nominal overlays with the g4.10.4 one. The
weight is the CV ratio of the two in true energy and angle to the beamline,
/data/uboone/NuMIFlux/NewFluxFiles/ratios_g4_10_4_cvonly_fornuwrofakedata.root ("New Flux / Old (bugged) Flux"), per
flavour and horn mode; the angle is that of the true neutrino direction in beam coordinates, with the rotation of
src/app/NuMI/addBeamlineGeometryWeightsToMap.cpp (GetNuMIAngle, "beam"). The true direction exists only in the
unprocessed ntuples; their entries are those of stv_tree in the same order (checked on the flavour and energy).
Events outside the populated ratio map (0.3%) keep weight 1.

    python3 scripts/detvar_fluxw.py <line of slurm/mp2pi_detvar_manifest.list, from 0>
writes /data/uboone/processed/mp2pi/detvar_fluxw/<name>.root through macros/add_fluxw_branch.C.
"""
import os, subprocess, sys
import numpy as np
import uproot

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
# another production: MP2PI_DETVAR_MANIFEST (e.g. slurm/mp2pi_mcscal_detvar_manifest.list) and MP2PI_DIR
MAN = os.environ.get('MP2PI_DETVAR_MANIFEST', os.path.join(REPO, 'slurm', 'mp2pi_detvar_manifest.list'))
OUT = os.path.join(os.environ.get('MP2PI_DIR', '/data/uboone/processed/mp2pi'), 'detvar_fluxw', '')
RATIO = '/data/uboone/NuMIFlux/NewFluxFiles/ratios_g4_10_4_cvonly_fornuwrofakedata.root'
# detector -> beam rotation; TRotation::RotateAxes takes the three vectors as the columns
ROT = np.column_stack([[0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021],
                       [4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359],
                       [-0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291]])
FLAV = {14: 'numu', -14: 'numubar', 12: 'nue', -12: 'nuebar'}


def weights(raw_files, proc_file, mode):
    br = ['nu_pdg', 'nu_e', 'true_nu_px', 'true_nu_py', 'true_nu_pz']
    raw = {b: [] for b in br}
    for fn in raw_files:
        a = uproot.open(fn)['nuselection/NeutrinoSelectionFilter'].arrays(br, library='np')
        for b in br: raw[b].append(a[b])
    raw = {b: np.concatenate(v) for b, v in raw.items()}
    p = uproot.open(proc_file)['stv_tree'].arrays(['mc_nu_pdg', 'mc_nu_energy'], library='np')
    if len(p['mc_nu_pdg']) != len(raw['nu_pdg']) or not np.array_equal(p['mc_nu_pdg'], raw['nu_pdg']) \
            or not np.allclose(p['mc_nu_energy'], raw['nu_e']):
        sys.exit(f'{proc_file}: entries do not match the unprocessed ntuples')
    beam = ROT @ np.vstack([raw['true_nu_px'], raw['true_nu_py'], raw['true_nu_pz']])
    ang = np.degrees(np.arccos(np.clip(beam[2] / np.maximum(np.linalg.norm(beam, axis=0), 1e-12), -1., 1.)))
    w = np.ones(len(ang), dtype=np.float32)
    rf = uproot.open(RATIO)
    for pdg, fl in FLAV.items():
        h = rf[f'ratio_cv_{mode}_enu_angle_{fl}']
        v, ex, ey = h.values(), h.axis(0).edges(), h.axis(1).edges()
        m = raw['nu_pdg'] == pdg
        ix = np.searchsorted(ex, raw['nu_e'][m], side='right') - 1
        iy = np.searchsorted(ey, ang[m], side='right') - 1
        inside = (ix >= 0) & (ix < len(ex) - 1) & (iy >= 0) & (iy < len(ey) - 1)
        r = np.ones(m.sum())
        r[inside] = v[ix[inside], iy[inside]]
        r[~(r > 0)] = 1.
        w[m] = r
    return w


def main():
    k = int(sys.argv[1])
    raws, var, _, proc = [l.strip().split('|') for l in open(MAN) if l.strip()][k]
    mode = 'fhc' if 'run4fhc' in proc else 'rhc'
    os.makedirs(OUT, exist_ok=True)
    out = OUT + os.path.basename(proc)
    w = weights(raws.split('+'), proc, mode)
    wfile = out + '.fluxw.f32'
    w.tofile(wfile)
    print(f'{var} {mode}: {len(w)} events, mean weight {w.mean():.4f}, at 1 (outside the map or other flavours) '
          f'{np.mean(w == 1.):.4f}', flush=True)
    rc = subprocess.call(['root', '-l', '-b', '-q', f'{REPO}/macros/add_fluxw_branch.C("{proc}","{wfile}","{out}")'])
    sys.exit(rc)


if __name__ == '__main__':
    main()
