# binning_eval_c50.py -- secondary gate of the 0.50-criterion study: for every candidate built by
# slurm/slurm_c50.sbatch into processed/rebuild_c50, the quantities of tab:binning_cands --
# A_C row sums, A_C conditioning s_min/s_1, bin-averaged data-statistical term, per-bin total
# uncertainty range and closure ratio range (the minimum diagonals are in logs/c50/screen_*.tsv).
# Candidates whose outputs are missing are listed as pending.
#   python3 report/tools/binning_eval_c50.py
import re, sys, os, ROOT
ROOT.gROOT.SetBatch(True)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', 'xsec_analyzer'))
RB = os.environ.get('C50_RB', '/data/uboone/processed/rebuild_c50/')
sys.path.insert(0, os.path.join(REPO, 'scripts'))

def dump(path):
    d = {}; chi = None
    for l in open(path):
        if l.startswith('[SYSTDUMP]'):
            p = l.split()
            if len(p) >= 3:
                try: d[p[1]] = float(p[2])
                except ValueError: pass
        m = re.search(r'truth: .* = ([0-9.e+-]+)/(\d+) bins?, p-value = ([0-9.e+-]+)', l)
        if m and chi is None: chi = (float(m.group(1)), int(m.group(2)), float(m.group(3)))
    return d, chi

def side(path):
    f = ROOT.TFile.Open(path); u = f.Get('h_unfolded_nuwro'); a = f.Get('h_A_C'); nb = u.GetNbinsX()
    rows = [sum(a.GetBinContent(j, i) for j in range(1, nb + 1)) for i in range(1, nb + 1)]
    M = ROOT.TMatrixD(nb, nb)
    for i in range(nb):
        for j in range(nb): M[i][j] = a.GetBinContent(j + 1, i + 1)
    sv = ROOT.TDecompSVD(M); sv.Decompose(); s = sv.GetSig()
    cond = min(s[k] for k in range(nb)) / s[0] if s[0] > 0 else float('nan')
    rel = [100 * u.GetBinError(i) / abs(u.GetBinContent(i)) if u.GetBinContent(i) else float('nan') for i in range(1, nb + 1)]
    clo = []
    t = f.Get('h_fakedata_truth')  # already A_C-smeared: CrossSectionExtractor smears pred_map in place
    if t:
        clo = [u.GetBinContent(i) / t.GetBinContent(i) for i in range(1, nb + 1) if t.GetBinContent(i)]
    f.Close(); return nb, rows, cond, rel, clo

def main():
    from c50_candidates import CANDS
    print(f"{'candidate':12s} {'cfg':8s} {'nb':>2s} {'A_C rows':>11s} {'smin/s1':>7s} {'DataStat':>8s} {'bin unc %':>11s} {'closure':>11s} {'chi2/ndf p':>14s}")
    for c in CANDS:
        # candidates may declare their own configurations (the RHC theta_p schemes add comb)
        want = c.get('cfgs', ('fhc5', 'rhcfull'))
        for cfg, T in [(x, {'fhc5': 'FHC5', 'rhcfull': 'RHCFULL', 'comb': 'COMB'}[x]) for x in want]:
            base = f"{c['pfx']}_{T}_{c['name']}"
            log, out = RB + f'unfold_{base}.log', RB + f'closure_hists_xsec_{base}.root'
            if not (os.path.exists(log) and os.path.exists(out)):
                print(f"{c['name']:12s} {T:8s} pending"); continue
            try:
                d, chi = dump(log); nb, rows, cond, rel, clo = side(out)
            except Exception as e:
                print(f"{c['name']:12s} {T:8s} unreadable ({e})"); continue
            cl = f"{min(clo):.2f}-{max(clo):.2f}" if clo else 'n/a'
            ch = f"{chi[0]:.2f}/{chi[1]} {chi[2]:.2f}" if chi else 'n/a'
            print(f"{c['name']:12s} {T:8s} {nb:2d} {min(rows):.2f}-{max(rows):.2f} {cond:7.3f} {d.get('DataStats', float('nan')):7.1f}% "
                  f"{min(rel):4.0f}-{max(rel):4.0f}   {cl:>11s} {ch:>14s}")

if __name__ == '__main__':
    main()
