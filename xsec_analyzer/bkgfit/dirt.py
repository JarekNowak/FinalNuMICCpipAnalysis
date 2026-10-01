"""dirt.py -- per-run dirt normalisation, the Python counterpart of macros/dirt_perrun.h (2026-10-01).

Each run's dirt sample(s) are scaled by that run's data POT over the summed_pot of the sample(s)
(the sum of the DISTINCT per-file values), as SystematicsCalculator scales dirtMC. Runs, data POT and
dirt samples are read from the release file-properties list. The scale does NOT contain the 0.65
dirt normalisation: ProcessNTuples stores it in normalisation_weight, so it enters once, with the
central-value weight.

Until 2026-10-01 farsb.py and run1_trigger.py scaled the Run-1 dirt sample alone by the data POT
over the summed OVERLAY POT (0.092402 x POT/8.857 in FHC, 0.071666 x POT/11.082 in RHC) and applied
the 0.65 a second time.
"""
import os
import uproot

XA = os.environ.get('XSEC_ANALYZER_DIR', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
FP = os.path.join(XA, 'configs', 'file_properties_numi_comb_w.txt')

# release alias (xsec-ana-<alias>.root) -> name of the processed file it stands for, so that
# directories holding only the processed names (e.g. a fresh reprocess) can be read as well
ORIG = {
    'dirt_fhc_run1': 'prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot',
    'dirt_fhc_run2': 'prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot',
    'dirt_fhc_run4c': 'numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana',
    'dirt_fhc_run4d': 'numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana',
    'dirt_fhc_run5': 'run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns',
    'dirt_rhc_run1': 'neutrinoselection_filt_run3b_dirt_overlay',
    'dirt_rhc_run2': 'neutrinoselection_filt_run3b_dirt_overlay',
    'dirt_rhc_run3': 'neutrinoselection_filt_run3b_dirt_overlay',
    'dirt_rhc_run4a': 'run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana',
    'dirt_rhc_run4b': 'numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana',
}


def read_fp(fp=FP):
    """({run: data POT}, {run: [dirt file basenames]}) from the file-properties list."""
    pot, files = {}, {}
    for line in open(fp):
        tok = line.split('#')[0].split()
        if len(tok) < 3: continue
        path, run, typ = tok[0], int(tok[1]), tok[2]
        if typ == 'onBNB': pot[run] = float(tok[4])
        elif typ == 'dirtMC': files.setdefault(run, []).append(os.path.basename(path))
    return pot, files


def dirt_perrun(mode='comb', where=None, skip=(), fp=FP):
    """[(run, path, scale)] for mode 'fhc' (runs 1-5), 'rhc' (runs 11-14) or 'comb'.

    where(run) gives the directory to read that run's samples from (default /data/uboone/processed/);
    the alias name is tried first, then the processed-file name. The POT is read from the file read.
    skip: runs to leave out (e.g. (2, 12): no Run-2 beam-on data)."""
    pot, files = read_fp(fp)
    out = []
    for run in sorted(files):
        fhc = run < 10
        if (mode == 'fhc' and not fhc) or (mode == 'rhc' and fhc) or run in skip: continue
        d = where(run) if where else '/data/uboone/processed/'
        paths = []
        for alias in files[run]:
            stem = alias[len('xsec-ana-'):-len('.root')]
            cands = [os.path.join(d, alias), os.path.join(d, 'xsec-ana-' + ORIG.get(stem, stem) + '.root')]
            p = next((c for c in cands if os.path.exists(c)), None)
            if p is None: raise FileNotFoundError(f'dirt run {run}: none of {cands}')
            paths.append(p)
        seen = []
        for p in paths:
            v = float(uproot.open(p)['summed_pot'].member('fVal'))
            if not any(abs(s - v) <= 1e-6 * abs(s) for s in seen): seen.append(v)
        if run not in pot or sum(seen) <= 0: raise ValueError(f'dirt run {run}: no data POT or dirt POT')
        out += [(run, p, pot[run] / sum(seen)) for p in paths]
    return out


if __name__ == '__main__':
    for mode in ('fhc', 'rhc'):
        for run, p, s in dirt_perrun(mode):
            print(f'{mode} run {run:2d}  scale {s:.5f}  {os.path.basename(p)}')
