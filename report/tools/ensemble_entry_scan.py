"""ensemble_entry_scan.py TAG [NMAX] -- the entry-count scan of an ensemble (2026-09-24 finding, now a tool):
a member whose universe file was built from a truncated throw shows fewer entries in the throw's
unweighted_0_true histogram (all bins, flow included) than the throw file holds on disk. Compares
every throw of every finished member of TAG (<cfg>_<obs> or 1p_<cfg>_<obs>) with its file in
ens/throws/ (w/ens/throws/ for the proton-tagged family) and reports any mismatch. A truncated read is
short by hundreds or thousands of entries; a shortfall of one or two is the boundary case (five overlay
signal events have cos theta_mu exactly 1.0 and fall outside a top true bin written as "< 1.000";
2026-09-26) and is reported as such. Exit status 1 only for a truncated read.
    python3 report/tools/ensemble_entry_scan.py fhc5_pmu [100]
"""
import os, sys, uproot
P = '/data/uboone/processed/'
tag = sys.argv[1]; NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 100
TH = P + ('w/ens/throws/' if tag.startswith('1p_') else 'ens/throws/')
ok = short = 0; bad = []
ntree = {}
for t in range(1, NMAX + 1):
    u = P + f'ens/univ_{tag}_t{t}.root'
    if not os.path.exists(u) or not os.path.exists(P + f'ens/closure_hists_xsec_{tag}_t{t}.root'): continue
    f = uproot.open(u); d = f[f.keys()[0].split(';')[0]]
    throws = [k.split(';')[0] for k in d.keys() if '+throws+fakedata_' in k and '/' not in k]
    for k in throws:
        fn = k.replace('+', '/')
        if fn not in ntree: ntree[fn] = uproot.open(fn)['stv_tree'].num_entries
        h = d[k + '/unweighted_0_true']; n = int(round(h.values(flow=True).sum()))
        if n == ntree[fn]: ok += 1
        else: short += 1; bad.append((t, os.path.basename(fn), n, ntree[fn]))
trunc = [b for b in bad if b[3] - b[2] > 2]
print(f'{tag}: {ok} throws match, {len(bad) - len(trunc)} short by one or two entries (boundary case), {len(trunc)} truncated')
for t, fn, n, m in bad: print(f'  member {t} {fn}: read {n} of {m} ({100*n/m:.1f}%){"  TRUNCATED" if m - n > 2 else ""}')
sys.exit(1 if trunc else 0)
