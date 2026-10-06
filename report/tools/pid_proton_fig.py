"""pid_proton_fig.py -- proton rejection of the pion candidate (inclusive note, fig:pid_proton), from the
pion-candidate track dumps of the 14 release overlays (/data/uboone/processed/mp_pid/pidtrk-*.root,
macros/mp_pid/dump_pid_tracks.C on multipion/phase1).

Population: the tracks the pion identification acts on: generation 2, track score >= 0.5, not the muon
candidate, end contained, start within 4 cm of the vertex, longer than 20 cm, in events with a muon
candidate; truth-matched (purity > 0.5) pions and protons, central-value weighted, area-normalised.
  figures/pid_pibdt.png   pion-BDT score, with the applied threshold (> -0.10)
  figures/pid_bpion.png   Bragg-pion score, not applied (absent from the older-production ntuples)

    python3 report/tools/pid_proton_fig.py
"""
import glob, os
import numpy as np
import uproot
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(os.path.dirname(HERE), 'figures')
B = ['pdg', 'purity', 'is_mu_cand', 'contained', 'dist', 'len', 'ts', 's_pion', 'bragg_pion', 'w_cv', 'entry']


def load():
    pi, pr = {'s_pion': [], 'bragg_pion': [], 'w': []}, {'s_pion': [], 'bragg_pion': [], 'w': []}
    for f in sorted(glob.glob('/data/uboone/processed/mp_pid/pidtrk-*.root')):
        a = uproot.open(f)['trk'].arrays(B, library='np')
        key = a['entry']                                   # event = entry of the raw tree (one sample per file)
        mu_events = np.unique(key[a['is_mu_cand'] == 1])
        m = np.isin(key, mu_events) & (a['is_mu_cand'] == 0) & (a['ts'] >= 0.5) & (a['contained'] == 1) \
            & (a['dist'] <= 4) & (a['len'] > 20) & (a['purity'] > 0.5)
        w = np.where(np.isfinite(a['w_cv']) & (a['w_cv'] >= 0) & (a['w_cv'] <= 30), a['w_cv'], 1.)
        for d, sel in ((pi, m & (np.abs(a['pdg']) == 211)), (pr, m & (a['pdg'] == 2212))):
            for v in ('s_pion', 'bragg_pion'): d[v].append(a[v][sel])
            d['w'].append(w[sel])
    return ({k: np.concatenate(v) for k, v in pi.items()}, {k: np.concatenate(v) for k, v in pr.items()})


def draw(pi, pr, var, xlab, lo, hi, nb, out, thr=None):
    hp = ROOT.TH1F('hpi_' + var, f';{xlab};Area-normalised', nb, lo, hi)
    hr = ROOT.TH1F('hpr_' + var, f';{xlab};Area-normalised', nb, lo, hi)
    for h, d in ((hp, pi), (hr, pr)):
        x = np.clip(d[var], lo, hi - 1e-9)
        h.FillN(len(x), x.astype(float), d['w'].astype(float))
        h.Scale(1. / h.Integral())
    hp.SetLineColor(ROOT.kRed + 1); hp.SetFillColorAlpha(ROOT.kRed - 9, 0.5); hp.SetLineWidth(2)
    hr.SetLineColor(ROOT.kAzure + 2); hr.SetFillColorAlpha(ROOT.kAzure - 9, 0.5); hr.SetLineWidth(2)
    c = ROOT.TCanvas('c_' + var, '', 720, 540)
    ymax = max(hp.GetMaximum(), hr.GetMaximum()) * 1.25
    hp.SetMaximum(ymax); hp.SetMinimum(0); hp.Draw('hist'); hr.Draw('hist same')
    l = ROOT.TLegend(0.62, 0.72, 0.89, 0.89)
    l.AddEntry(hp, 'true #pi^{#pm}', 'f'); l.AddEntry(hr, 'true proton', 'f')
    l.SetBorderSize(0); l.SetFillStyle(0); l.Draw()
    if thr is not None:
        ln = ROOT.TLine(thr, 0, thr, ymax * 0.85); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw()
        tx = ROOT.TLatex(); tx.SetTextSize(0.034); tx.DrawLatex(thr, ymax * 0.88, 'cut')
    c.SaveAs(os.path.join(FIG, out))
    return hp, hr


def main():
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
    pi, pr = load()
    f = lambda d, cut: float(d['w'][cut].sum() / d['w'].sum())
    print(f"tracks: {len(pi['w'])} pions, {len(pr['w'])} protons")
    print(f"pion BDT > -0.10: pions {100 * f(pi, pi['s_pion'] > -0.1):.1f}%, protons {100 * f(pr, pr['s_pion'] > -0.1):.1f}%")
    print(f"Bragg-pion >= 0.08: pions {100 * f(pi, pi['bragg_pion'] >= 0.08):.1f}%, protons {100 * f(pr, pr['bragg_pion'] >= 0.08):.1f}%")
    draw(pi, pr, 's_pion', 'Pion BDT score', -0.4, 0.3, 35, 'pid_pibdt.png', thr=-0.1)
    draw(pi, pr, 'bragg_pion', 'Bragg-pion score', 0., 1., 50, 'pid_bpion.png')


if __name__ == '__main__':
    main()
