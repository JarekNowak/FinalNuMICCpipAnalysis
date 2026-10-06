"""final_topo_fig.py -- the topological variables at the final cut (inclusive note, fig:final_topo): the
muon-pion opening angle and the number of primary tracks of the selected events, the per-run
pseudo-data against the simulation, FHC at its full exposure (7.766e20 POT).

Samples and normalisation are those of the release (configs/file_properties_numi_fhc5.txt), as the
framework applies them: overlay and dirt per run by data POT over the summed_pot of that run's files
(distinct values), beam-off by 0.98 x triggers over gates, central-value weight
safe(tune x ppfx x normalisation). The pseudo-data are the release's per-run Poisson throws.
  figures/final_oa.{pdf,png}, figures/final_ntracks.{pdf,png}

    python3 report/tools/final_topo_fig.py
"""
import os
import numpy as np
import uproot
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
REP = os.path.dirname(HERE)
FIG = os.path.join(REP, 'figures')
LIST = os.path.join(REP, '..', 'xsec_analyzer', 'configs', 'file_properties_numi_fhc5.txt')
S = 'CC1mu1piXp'
VARS = {'oa': (f'{S}_mu_pi_opening_angle', '#mu#minus#pi opening angle [rad]', 32, 0., 3.2, 'final_oa'),
        'nt': (f'{S}_sb_nprimtrk', 'Number of primary tracks', 4, 1.5, 5.5, 'final_ntracks')}
W = ['tuned_cv_weight', 'ppfx_cv_weight', 'normalisation_weight']


def parse():
    on, rows = {}, []
    for line in open(LIST):
        t = line.split('#')[0].split()
        if len(t) < 3: continue
        path, run, typ = t[0], int(t[1]), t[2]
        if typ == 'onBNB': on[run] = (float(t[3]), float(t[4]))
        rows.append((path, run, typ, float(t[3]) if typ == 'extBNB' else 0.))
    return on, rows


def scales(on, rows):
    pot = lambda p: float(uproot.open(p)['summed_pot'].member('fVal'))
    out = []
    for run in sorted(on):
        for typ in ('numuMC', 'dirtMC'):
            fs = [p for p, r, ty, _ in rows if r == run and ty == typ]
            distinct = []
            for p in fs:
                v = pot(p)
                if not any(abs(v - d) <= 1e-6 * abs(d) for d in distinct): distinct.append(v)
            out += [(p, run, typ, on[run][1] / sum(distinct)) for p in fs]
        ext = [(p, g) for p, r, ty, g in rows if r == run and ty == 'extBNB']
        out += [(p, run, 'extBNB', 0.98 * on[run][0] / sum(g for _, g in ext)) for p, _ in ext]
        out += [(p, run, 'onBNB', 1.) for p, r, ty, _ in rows if r == run and ty == 'onBNB']
    return out


def main():
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
    on, rows = parse()
    H = {}
    for key, (br, xl, nb, lo, hi, out) in VARS.items():
        for c in ('sig', 'bkg', 'oofv', 'dirt', 'ext', 'data'):
            H[(key, c)] = ROOT.TH1D(f'h_{key}_{c}', f';{xl};Events (7.766#times10^{{20}} POT)', nb, lo, hi)
            H[(key, c)].Sumw2()
    for path, run, typ, sc in scales(on, rows):
        t = uproot.open(path)['stv_tree']
        br = [f'{S}_Selected', f'{S}_MC_Signal', f'{S}_EventCategory'] + [v[0] for v in VARS.values()]
        a = t.arrays(br + (W if typ in ('numuMC', 'dirtMC') else []), library='np')
        sel = a[f'{S}_Selected'].astype(bool)
        if typ in ('numuMC', 'dirtMC'):
            w = np.prod([a[x].astype(float) for x in W], axis=0)
            w = np.where(np.isfinite(w) & (w >= 0) & (w <= 30), w, 1.) * sc
        else:
            w = np.full(len(sel), sc)
        if typ == 'numuMC':
            sig = a[f'{S}_MC_Signal'].astype(bool); oofv = a[f'{S}_EventCategory'] == 2
            cats = {'sig': sel & sig, 'oofv': sel & ~sig & oofv, 'bkg': sel & ~sig & ~oofv}
        else:
            cats = {{'dirtMC': 'dirt', 'extBNB': 'ext', 'onBNB': 'data'}[typ]: sel}
        for key, (b, *_r) in VARS.items():
            x = a[b].astype(float)
            for c, m in cats.items():
                if m.any(): H[(key, c)].FillN(int(m.sum()), x[m], w[m])
        print(f'  run {run:2d} {typ:7s} scale {sc:.4f} selected {sel.sum():6d} {os.path.basename(path)}', flush=True)
    col = {'sig': ROOT.kAzure + 1, 'bkg': ROOT.kOrange + 1, 'oofv': ROOT.kGreen + 2, 'dirt': ROOT.kViolet - 4, 'ext': ROOT.kGray + 1}
    lab = {'sig': 'Signal (CC1#pi^{#pm})', 'bkg': 'Beam background', 'oofv': 'Out of FV', 'dirt': 'Dirt', 'ext': 'Beam-off'}
    for key, (b, xl, nb, lo, hi, out) in VARS.items():
        hs = ROOT.THStack(f'hs_{key}', f';{xl};Events (7.766#times10^{{20}} POT)')
        for c in ('sig', 'bkg', 'oofv', 'dirt', 'ext'):
            h = H[(key, c)]; h.SetFillColor(col[c]); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1); hs.Add(h)
        d = H[(key, 'data')]; d.SetMarkerStyle(20); d.SetMarkerSize(0.85); d.SetLineColor(ROOT.kBlack)
        c = ROOT.TCanvas(f'c_{key}', '', 720, 540)
        hs.SetMaximum(max(hs.GetMaximum(), d.GetMaximum() + d.GetBinError(d.GetMaximumBin())) * 1.35)
        hs.Draw('hist'); d.Draw('E1 same')
        l = ROOT.TLegend(0.58, 0.58, 0.89, 0.89); l.SetBorderSize(0); l.SetFillStyle(0)
        l.AddEntry(d, 'Pseudo-data', 'lep')
        for k in ('sig', 'bkg', 'oofv', 'dirt', 'ext'): l.AddEntry(H[(key, k)], lab[k], 'f')
        l.Draw()
        for ext in ('pdf', 'png'): c.SaveAs(os.path.join(FIG, f'{out}.{ext}'))
        tot = sum(H[(key, k)].Integral() for k in ('sig', 'bkg', 'oofv', 'dirt', 'ext'))
        print(f'{out}: prediction {tot:.1f}, pseudo-data {d.Integral():.0f}, signal {H[(key, "sig")].Integral():.1f}')


if __name__ == '__main__':
    main()
