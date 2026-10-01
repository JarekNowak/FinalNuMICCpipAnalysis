#!/usr/bin/env python3
"""Phase 0 baseline of the NuMI multi-pion selections (report/MULTIPION_BNB_ADAPTATION_PLAN.md).

Reads the single-pass reprocessing in /data/uboone/processed/mp_phase0, in which the
inclusive (CC1mu1piXp), two-pion (CC1mu2pi) and three-pion (CC1mu3pi) selections were run
into one tree, and normalises it as SystematicsCalculator does:

  numuMC, dirtMC : per run, data POT / sum of the DISTINCT per-file summed_pot of that type,
                   event weight safe(spline * tune * ppfx_cv * normalisation)
  extBNB         : per run, 0.98 * beam-on triggers / summed beam-off gates, unweighted

Exposures, triggers and gate counts are read from the release file-properties list
(configs/file_properties_numi_comb_w.txt); nothing is hard-coded. Every reprocessed MC and
dirt file is checked against the summed_pot of the release copy it replaces.

Outputs (report/multipion/): phase0_baseline.json (all sums), phase0_cutflow.tsv.
Usage: python3 scripts/mp_phase0_baseline.py
"""
import json, os, sys
import numpy as np
import uproot
import awkward as ak

REPO = '/home/t2k/nowak/MicroBooNE/working_xsec_analyzer'
CONF = REPO + '/xsec_analyzer/configs/file_properties_numi_comb_w.txt'
NEW = '/data/uboone/processed/mp_phase0/'
OUTDIR = REPO + '/report/multipion/'
FHC_RUNS, RHC_RUNS = (1, 2, 4, 5), (11, 12, 13, 14)
CONFIGS = {'FHC': FHC_RUNS, 'RHC': RHC_RUNS, 'COMB': FHC_RUNS + RHC_RUNS}
SELS = {'CC1mu2pi': 2, 'CC1mu3pi': 3}
# release dirt files that are regular files under an alias name -> raw sample they came from
DIRT_ALIAS = {
    'xsec-ana-dirt_fhc_run4c.root': 'xsec-ana-numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_ana.root',
    'xsec-ana-dirt_fhc_run4d.root': 'xsec-ana-numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_ana.root',
    'xsec-ana-dirt_fhc_run5.root': 'xsec-ana-run5_numi_fhc_dirt_overlay_pandora_ntuple_v08_00_00_67_slim_run5_ana_nonzerolifetime_goodruns.root',
    'xsec-ana-dirt_rhc_run3b.root': 'xsec-ana-neutrinoselection_filt_run3b_dirt_overlay.root',
    'xsec-ana-dirt_rhc_run4a.root': 'xsec-ana-run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_ana.root',
    'xsec-ana-dirt_rhc_run4b.root': 'xsec-ana-numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_ana.root',
}
# cumulative cut-flow of the multi-pion selections as applied (bits of <sel>_cutflow_bits);
# stages 6-7 are diagnostics, not part of the selection
STAGES = [('None', 0), ('Vertex in FV', (1 << 0) | (1 << 1)), ('Topological', 1 << 2),
          ('MuonCandidate', 1 << 3), ('N pions', 1 << 4), ('Opening angle', 1 << 8),
          ('Final (selected)', 1 << 9),
          ('diag: + shower veto', 1 << 7), ('diag: + 3-plane mu, pi', (1 << 5) | (1 << 6))]
N_SEL_STAGES = 7   # stages 0..6 make up the selection; the rest are diagnostics
TOPO = ['signal', 'CC0pi', 'CC1pi', 'CC2pi', 'CC3pi+', 'CCpi0', 'CCother', 'NC', 'nueCC', 'OOFV']


def summed_pot(path):
    return float(uproot.open(path)['summed_pot'].member('fVal'))


def parse_conf():
    on, ext, mc, dirt = {}, {}, {}, {}
    for raw in open(CONF):
        line = raw.split('#')[0].strip()
        if not line:
            continue
        f = line.split()
        path, run, ft = f[0], int(f[1]), f[2]
        if ft == 'onBNB':
            on[run] = (float(f[3]), float(f[4]))          # triggers, POT
        elif ft == 'extBNB':
            ext.setdefault(run, []).append((NEW + os.path.basename(os.path.realpath(path)), float(f[3])))
        elif ft == 'numuMC':
            mc.setdefault(run, []).append((NEW + os.path.basename(path), path))
        elif ft == 'dirtMC':
            real = os.path.realpath(path)
            b = os.path.basename(real)
            dirt.setdefault(run, []).append((NEW + DIRT_ALIAS.get(b, b), real))
    return on, ext, mc, dirt


def scales(on, ext, mc, dirt):
    """(run, kind) -> list of (new file, scale); checks summed_pot against the release copies."""
    out, checks = {}, []
    for kind, table in (('mc', mc), ('dirt', dirt)):
        for run, files in table.items():
            seen = []
            for new, old in files:
                pn, po = summed_pot(new), summed_pot(old)
                checks.append((kind, run, os.path.basename(new), pn, po))
                if not any(abs(s - pn) <= 1e-6 * abs(s) for s in seen):
                    seen.append(pn)
            denom = sum(seen)
            out[(run, kind)] = [(new, on[run][1] / denom) for new, _ in files]
    for run, files in ext.items():
        gates = sum(g for _, g in files)
        out[(run, 'ext')] = [(new, 0.98 * on[run][0] / gates) for new, _ in files]
    return out, checks


def branches():
    b = ['spline_weight', 'tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight',
         'CC1mu1piXp_Selected', 'CC1mu1piXp_MC_Signal', 'CC1mu1piXp_sb_multipi']
    per = ['Selected', 'MC_Signal', 'MC_Signal_1p', 'mc_topology', 'cutflow_bits',
           'n_contained_pion', 'n_uncontained_pion', 'pion_number_reco', 'pic_true_pdg',
           'pic_true_mom', 'pic_idx', 'CandidatePionIndex', 'mc_pi_mom',
           'mc_n_threshold_pionpm', 'mc_pionpm_min_mom', 'mc_n_threshold_pion0', 'mc_n_kaons',
           'mc_n_heaviermeson', 'mc_n_threshold_muon', 'candidate_muon_mom_true',
           'true_mu_pi_opening_angle', 'sig_truevertex_in_fv', 'sig_ccnc', 'sig_is_numu',
           'mc_n_pipm_thr0175']
    for s in SELS:
        b += [s + '_' + p for p in per]
    return b


def lead_match(true_pdg, true_mom, lead_mom, pos):
    """Candidate at position pos (per event, -1 = none) is the true leading pion (1% momentum)."""
    ok = np.zeros(len(pos), dtype=bool)
    for i, p in enumerate(pos):
        if p < 0 or p >= len(true_pdg[i]) or len(lead_mom[i]) == 0:
            continue
        lm = lead_mom[i][0]
        ok[i] = abs(true_pdg[i][p]) == 211 and abs(true_mom[i][p] - lm) <= 0.01 * lm
    return ok


def process(path, kind):
    """Unscaled weighted sums for one file. kind: mc | ext | dirt."""
    acc = {}

    def add(key, w):
        acc[key] = acc.get(key, 0.) + float(w)

    for d in uproot.iterate(path + ':stv_tree', branches(), step_size='400 MB', library='ak'):
        if kind == 'ext':
            w = np.ones(len(d['spline_weight']))
        else:
            w = (np.asarray(d['spline_weight'], dtype=np.float64) * np.asarray(d['tuned_cv_weight'], dtype=np.float64)
                 * np.asarray(d['ppfx_cv_weight'], dtype=np.float64) * np.asarray(d['normalisation_weight'], dtype=np.float64))
            w = np.where(np.isfinite(w) & (w >= 0.) & (w <= 30.), w, 1.0)
        sel1 = np.asarray(d['CC1mu1piXp_Selected'], bool)
        sig1 = np.asarray(d['CC1mu1piXp_MC_Signal'], bool) & (kind == 'mc')
        cr1 = np.asarray(d['CC1mu1piXp_sb_multipi'], bool)
        # 1pi validation against the released cut-flow (final stage = Selected)
        for lab, m in (('sig', sel1 & sig1), ('bkg', sel1 & ~sig1)) if kind == 'mc' else ((kind, sel1),):
            add(('v1pi', lab), w[m].sum())
        add(('v1pi', 'cr_multipi', kind), w[cr1].sum())
        sels = {}
        for s, n in SELS.items():
            g = lambda k: np.asarray(d[s + '_' + k])
            sel = g('Selected').astype(bool)
            sels[s] = sel
            sig = g('MC_Signal').astype(bool) & (kind == 'mc')
            topo = g('mc_topology').astype(int)
            bits = g('cutflow_bits').astype(int)
            cat = topo if kind == 'mc' else None
            # cumulative cut-flow by category
            cum = np.ones(len(w), bool)
            for si, (_, mask) in enumerate(STAGES):
                cum = cum & ((bits & mask) == mask)
                if kind == 'mc':
                    for c in range(10):
                        m = cum & (cat == c)
                        if m.any():
                            add(('cf', s, si, TOPO[c]), w[m].sum())
                else:
                    add(('cf', s, si, kind), w[cum].sum())
            # final stage must reproduce Selected
            fin = np.ones(len(w), bool)
            for _, mask in STAGES[:N_SEL_STAGES]:
                fin &= ((bits & mask) == mask)
            assert np.array_equal(fin, sel), (path, s, 'cut-flow final != Selected')
            if kind != 'mc':
                add(('ovl', s, 'sel1pi', kind), w[sel & sel1].sum())
                add(('ovl', s, 'cr_multipi', kind), w[sel & cr1].sum())
                add(('unc', s, kind), w[sel & (g('n_uncontained_pion') > 0)].sum())
                continue
            sig1p = g('MC_Signal_1p').astype(bool)
            add(('sig1p', s, 'gen'), w[sig1p].sum()); add(('sig1p', s, 'sel'), w[sig1p & sel].sum())
            all175 = sig & (g('mc_n_pipm_thr0175') == n)
            add(('thr0175', s, 'gen'), w[all175].sum()); add(('thr0175', s, 'sel'), w[all175 & sel].sum())
            theta_only = (g('sig_truevertex_in_fv').astype(bool) & g('sig_ccnc').astype(bool) & g('sig_is_numu').astype(bool)
                          & (g('mc_n_threshold_muon') == 1) & (g('mc_n_threshold_pionpm') == n) & (g('mc_n_threshold_pion0') == 0)
                          & (g('mc_n_kaons') == 0) & (g('mc_n_heaviermeson') == 0) & (g('candidate_muon_mom_true') > 0.15)
                          & (g('mc_pionpm_min_mom') > 0.10) & (g('true_mu_pi_opening_angle') >= 2.6))
            assert not (theta_only & sig).any()
            add(('thetaonly', s, 'gen'), w[theta_only].sum()); add(('thetaonly', s, 'sel'), w[theta_only & sel].sum())
            for lab, m in (('sig', sig), ('bkg', ~sig)):
                add(('ovl', s, 'sel1pi', lab), w[sel & sel1 & m].sum())
                add(('ovl', s, 'cr_multipi', lab), w[sel & cr1 & m].sum())
                add(('unc', s, lab), w[sel & m & (g('n_uncontained_pion') > 0)].sum())
            # truth exclusivity with the inclusive signal
            add(('excl', s, 'sig1pi'), w[sig & sig1].sum())
            # leading-pion convention, selected signal events
            ss = sel & sig
            if ss.any():
                idx = np.flatnonzero(ss)
                tpdg = ak.to_list(d[s + '_pic_true_pdg'][idx]); tmom = ak.to_list(d[s + '_pic_true_mom'][idx])
                pidx = ak.to_list(d[s + '_pic_idx'][idx]); lmom = ak.to_list(d[s + '_mc_pi_mom'][idx])
                cpi = g('CandidatePionIndex')[idx]
                pos_llr = np.array([pl.index(c) if c in pl else -1 for pl, c in zip(pidx, cpi)])
                ww = w[idx]
                add(('lead', s, 'n'), ww.sum())
                add(('lead', s, 'longest'), ww[lead_match(tpdg, tmom, lmom, np.zeros(len(idx), int))].sum())
                add(('lead', s, 'llrcand'), ww[lead_match(tpdg, tmom, lmom, pos_llr)].sum())
                npi = np.array([sum(1 for x in tp if abs(x) == 211) for tp in tpdg])
                add(('lead', s, 'allNtruepi'), ww[npi >= n].sum())
        add(('ovl', '2pi3pi', kind), w[sels['CC1mu2pi'] & sels['CC1mu3pi']].sum())
        if kind == 'mc':
            s2 = np.asarray(d['CC1mu2pi_MC_Signal'], bool); s3 = np.asarray(d['CC1mu3pi_MC_Signal'], bool)
            add(('excl', '2pi3pi'), w[s2 & s3].sum())
    return acc


def main():
    on, ext, mc, dirt = parse_conf()
    sc, checks = scales(on, ext, mc, dirt)
    # summed_pot is a float accumulated over subruns, so copies of one sample written in
    # different passes agree only to a few parts per million (Run-1 dirt: 2e-6)
    bad = [c for c in checks if abs(c[3] - c[4]) > 1e-5 * abs(c[4])]
    for c in checks:
        print('POT check %-5s run %2d %-90s new %.6e release %.6e %s' % (c[0], c[1], c[2], c[3], c[4], 'OK' if c not in bad else 'MISMATCH'))
    if bad:
        sys.exit('summed_pot mismatch: the reprocessed files do not match the release inputs')
    files = {}
    for (run, kind), lst in sc.items():
        for new, _ in lst:
            files[new] = kind
    per_file = {}
    for f, kind in sorted(files.items()):
        print('reading', kind, os.path.basename(f), flush=True)
        per_file[f] = process(f, kind)
    totals = {}
    for cfg, runs in CONFIGS.items():
        tot = {}
        for (run, kind), lst in sc.items():
            if run not in runs:
                continue
            for new, scale in lst:
                for k, v in per_file[new].items():
                    tot[k] = tot.get(k, 0.) + scale * v
        totals[cfg] = tot
    os.makedirs(OUTDIR, exist_ok=True)
    json.dump({cfg: {'|'.join(map(str, k)): v for k, v in t.items()} for cfg, t in totals.items()},
              open(OUTDIR + 'phase0_baseline.json', 'w'), indent=1, sort_keys=True)
    with open(OUTDIR + 'phase0_cutflow.tsv', 'w') as out:
        out.write('config\tselection\tstage\tsignal\tnu_bkg\text\tdirt\tpred\teff_pct\tpur_pct\n')
        for cfg, t in totals.items():
            for s in SELS:
                gen = t.get(('cf', s, 0, 'signal'), 0.)
                for si, (name, _) in enumerate(STAGES):
                    sig = t.get(('cf', s, si, 'signal'), 0.)
                    bkg = sum(t.get(('cf', s, si, c), 0.) for c in TOPO[1:])
                    ex, di = t.get(('cf', s, si, 'ext'), 0.), t.get(('cf', s, si, 'dirt'), 0.)
                    pred = sig + bkg + ex + di
                    out.write('%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%.2f\n' % (
                        cfg, s, name, sig, bkg, ex, di, pred, 100 * sig / gen if gen else 0, 100 * sig / pred if pred else 0))
    print('wrote', OUTDIR + 'phase0_baseline.json', OUTDIR + 'phase0_cutflow.tsv')


if __name__ == '__main__':
    main()
