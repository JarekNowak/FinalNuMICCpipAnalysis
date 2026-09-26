"""binreview_tables.py -- the binning review of 2026-09-26: the binning each observable would use under the
nominal 0.68 and the relaxed 0.50 migration-diagonal criteria, chosen by a stated rule from every binning
evaluated in logs/binreview/eval.tsv, and the tables of the technical supplement.

Rule, applied separately for each criterion:
  1. admissible: worst column-normalised diagonal (the unfolder's [DIAGDUMP]) at or above the criterion in
     BOTH FHC and RHC, and every true bin holding at least the event floor of the family in selected FHC
     signal at data exposure (70 inclusive, the smallest bin of the released inclusive schemes; 50
     proton-tagged);
  2. information: the degrees of freedom for signal, the trace of A_C (the sum of the Wiener filter
     factors), averaged over FHC and RHC;
  3. choice: among admissible binnings, the fewest bins reaching 90% of the largest averaged trace; ties go
     to the larger worst diagonal. With no admissible binning of two or more bins the observable has no
     shape measurement under that criterion (one bin, i.e. the integrated cross section).
The choice is flagged as rule-sensitive when 85% or 95% instead of 90%, or averaging over FHC, RHC and
combined instead of FHC and RHC, would pick a different binning.

Writes tables/binreview_incl.tex, tables/binreview_1p.tex, tables/binreview_edges.tex and
logs/binreview/choices.tsv.
    python3 report/tools/binreview_tables.py
"""
import os, sys, csv, re, math, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import maximin_binning as MB
REP = os.path.dirname(HERE); TOP = os.path.dirname(REP)
EVAL = os.path.join(TOP, 'logs', 'binreview', 'eval.tsv')
FLOOR = {'incl': 70, '1p': 50}
CRIT = (0.68, 0.50)
# released edges (analysis note, tab:binning); theta_p and theta_pip are the proposed study binnings
PI = math.pi
REL_EDGES = {('incl', 'pmu'): [0.15, 0.35, 0.55, 0.75, 0.95, 1.25, 1.75, 3.0], ('incl', 'ppi'): [0.175, 0.205, 1.0],
             ('incl', 'costhmu'): [-1, 0.45, 0.65, 0.8, 0.9, 1], ('incl', 'thetamu'): [0, 0.43, 0.62, 0.84, 1.17, 3.15],
             ('incl', 'costhpi'): [-1, -0.1, 0.35, 0.6, 0.8, 1], ('incl', 'thmupi'): [0, 0.6, 0.85, 1.3, 1.85, 2.6],
             ('1p', 'Whad'): [0, 0.78, 1.08, 1.21, 1.31, 1.44, 2.74], ('1p', 'Wpipr'): [1.08, 1.19, 1.23, 1.27, 1.34, 1.47, 2.9],
             ('1p', 'dpt'): [0, 0.3, 2.5], ('1p', 'dalphat'): [0, 120, 180], ('1p', 'dphit'): [0, 50, 180],
             ('1p', 'pn'): [0, 0.375, 2.0], ('1p', 'thetap'): [0, 0.659, 1.068, 1.539, PI], ('1p', 'thpipr'): [0, 1.193, 2.010, PI]}
ORDER = {'incl': ['pmu', 'ppi', 'costhmu', 'thetamu', 'costhpi', 'thmupi'],
         '1p': ['Whad', 'Wpipr', 'dpt', 'dalphat', 'dphit', 'pn', 'thetap', 'thpipr']}
LAB = {'pmu': r'$p_\mu$', 'ppi': r'$p_\pi$', 'costhmu': r'$\cos\theta_\mu$', 'thetamu': r'$\theta_\mu$',
       'costhpi': r'$\cos\theta_\pi$', 'thmupi': r'$\theta_{\mu\pi}$', 'Whad': r'$W_{\rm had}$', 'Wpipr': r'$W_{\pi p}$',
       'dpt': r'$\delta p_T$', 'dalphat': r'$\delta\alpha_T$', 'dphit': r'$\delta\phi_T$', 'pn': r'$p_n$',
       'thetap': r'$\theta_p$', 'thpipr': r'$\theta_{\pi p}$'}
UNIT = {'pmu': 'GeV/$c$', 'ppi': 'GeV/$c$', 'thetamu': 'rad', 'thmupi': 'rad', 'Whad': 'GeV/$c^2$', 'Wpipr': 'GeV/$c^2$',
        'dpt': 'GeV/$c$', 'dalphat': 'deg', 'dphit': 'deg', 'pn': 'GeV/$c$', 'thetap': 'rad', 'thpipr': 'rad'}
CF = ('fhc5', 'rhcfull', 'comb')

rows = list(csv.DictReader(open(EVAL), delimiter='\t'))
G = collections.defaultdict(lambda: collections.defaultdict(dict))
for r in rows: G[(r['family'], r['observable'])][r['binning']][r['config']] = r

_mats = {}
def min_events(fam, obs, edges):
    """smallest selected-FHC-signal content of any true bin (data exposure), from the scan's fine matrices"""
    if fam not in _mats: _mats[fam] = {m: MB.load(fam, m) for m in ('fhc', 'rhc')}
    return min(MB.evaluate(fam, obs, edges, _mats[fam])['fhc'][1])

def edges_of(fam, obs, b, d):
    if b == 'released': return REL_EDGES[(fam, obs)]
    return [float(x) for x in d['fhc5']['edges'].split()]

def floor_of(fam, obs, b, d):
    if b == 'released': return min_events(fam, obs, REL_EDGES[(fam, obs)])
    m = re.search(r'floor (\d+)', d['fhc5'].get('note', '')); return float(m.group(1)) if m else float(FLOOR[fam])

def cands(fam, obs):
    out = []
    for b, d in G[(fam, obs)].items():
        if not all(c in d for c in CF): continue
        f = lambda c, k: float(d[c][k])
        out.append(dict(b=b, d=d, K=int(d['fhc5']['bins']), w=min(f('fhc5', 'min_diag'), f('rhcfull', 'min_diag')),
                        fl=floor_of(fam, obs, b, d), dFR=(f('fhc5', 'dfs') + f('rhcfull', 'dfs')) / 2,
                        dFRC=sum(f(c, 'dfs') for c in CF) / 3))
    return out

def choose(fam, obs, thr, tol=0.90, key='dFR'):
    adm = [c for c in cands(fam, obs) if c['w'] >= thr and c['fl'] >= FLOOR[fam] and c['K'] >= 2]
    if not adm: return None, 1
    dmax = max(c[key] for c in adm)
    ok = sorted((c for c in adm if c[key] >= tol * dmax - 1e-9), key=lambda c: (c['K'], -c['w']))
    return ok[0], max(c['K'] for c in adm)

def fmt(d, k, dig, cfgs=CF, pct=False):
    out = []
    for c in cfgs:
        v = float(d[c][k]); out.append(f'{v:.0f}' if pct else (f'{v:.{dig}f}' if dig else f'{int(v)}'))
    return '/'.join(out)

res = {}
for fam in ORDER:
    for obs in ORDER[fam]:
        rr = {}
        for thr in CRIT:
            c, kmax = choose(fam, obs, thr)
            alt = {choose(fam, obs, thr, t, k)[0]['b'] if choose(fam, obs, thr, t, k)[0] else 'one bin'
                   for t in (0.85, 0.90, 0.95) for k in ('dFR', 'dFRC')}
            rr[thr] = dict(c=c, kmax=kmax, sensitive=len(alt) > 1, alts=sorted(alt))
        rel = next((c for c in cands(fam, obs) if c['b'] == 'released'), None)
        res[(fam, obs)] = dict(rel=rel, **{str(t): rr[t] for t in CRIT})

# ---- machine-readable record: the working copy in logs/ and the dated copies kept with the report
REL = os.path.join(REP, 'data_release')
HDR = ('# Binning review of 2026-09-26 (0.68 versus 0.50 migration-diagonal criteria), fake data on the corrected\n'
       '# exposure; technical supplement sec:binreview. Study only: no released binning is changed.\n')
with open(EVAL) as f, open(os.path.join(REL, 'binreview_eval_2026-09-26.tsv'), 'w') as o:
    o.write(HDR + '# One row per binning and configuration (report/tools/binreview_eval.py); binning = released, study or\n'
            '# <obs>K<n> (maximin, n bins) / costhmuT<n> (theta_mu partition in cos). dfs = trace of A_C, n50/n20 = Wiener\n'
            '# filter factors (eigenvalues of A_C) >= 0.5 / >= 0.2, eig = the factors, min_diag = worst diagonal (DIAGDUMP).\n')
    o.write(f.read())
CHO = [os.path.join(TOP, 'logs', 'binreview', 'choices.tsv'), os.path.join(REL, 'binreview_choices_2026-09-26.tsv')]
with open(CHO[0], 'w') as o:
    o.write('family\tobservable\tcriterion\tchoice\tbins\tfinest_admissible\tworst_diag_FHC_RHC\ttrace_AC_FHC_RHC_COMB\tfactors_ge_0.5\t'
            'max_bin_unc_pct\trule_sensitive\talternatives\tedges\n')
    for (fam, obs), v in res.items():
        for thr in CRIT:
            r = v[str(thr)]; c = r['c']
            if c is None:
                o.write(f"{fam}\t{obs}\t{thr}\tone bin\t1\t{r['kmax']}\t\t\t\t\t{r['sensitive']}\t{','.join(r['alts'])}\t\n"); continue
            d = c['d']
            o.write(f"{fam}\t{obs}\t{thr}\t{c['b']}\t{c['K']}\t{r['kmax']}\t{fmt(d, 'min_diag', 3, CF[:2])}\t{fmt(d, 'dfs', 2)}\t"
                    f"{fmt(d, 'n50', 0)}\t{fmt(d, 'bin_unc_max', 0, pct=True)}\t{r['sensitive']}\t{','.join(r['alts'])}\t"
                    f"{' '.join(f'{round(e, 3):g}' for e in edges_of(fam, obs, c['b'], d))}\n")
with open(CHO[0]) as f, open(CHO[1], 'w') as o:
    o.write(HDR + '# Binning selected per observable and criterion by the rule of report/tools/binreview_tables.py.\n'); o.write(f.read())

# ---- LaTeX tables
def line(label, K, d, flag=''):
    if d is None:
        return f"    {label} & 1 & \\multicolumn{{4}}{{l}}{{no binning of two or more bins is admissible}} \\\\"
    return (f"    {label}{flag} & {K} & {fmt(d, 'min_diag', 2, CF[:2])} & {fmt(d, 'dfs', 1)} & {fmt(d, 'n50', 0)} & "
            f"{fmt(d, 'bin_unc_max', 0, pct=True)} \\\\")

def table(fam, fn, caption, label):
    L = [r'% generated by report/tools/binreview_tables.py -- do not edit', r'\begin{table}[H]', r'  \centering',
         r'  \caption{' + caption + '}', r'  \label{' + label + '}', r'  \footnotesize\setlength{\tabcolsep}{4pt}',
         r'  \begin{tabular}{llcccc}', r'    \toprule',
         r'    & bins & worst diag. & trace $A_C$ & factors $\geq0.5$ & largest bin unc.\ [\%] \\',
         r'    & & FHC/RHC & FHC/RHC/COMB & FHC/RHC/COMB & FHC/RHC/COMB \\', r'    \midrule']
    for i, obs in enumerate(ORDER[fam]):
        v = res[(fam, obs)]
        if i: L.append(r'    \midrule')
        L.append(r'    \multicolumn{6}{l}{' + LAB[obs] + r'} \\')
        rel = v['rel']
        tag = 'released' if not (fam == '1p' and obs in ('thetap', 'thpipr')) else 'study'
        if rel:
            note = '' if rel['fl'] >= FLOOR[fam] else f"$^{{\\ast}}$"
            L.append(line(f'\\quad {tag}', rel['K'], rel['d'], note))
        for thr in CRIT:
            r = v[str(thr)]; c = r['c']
            flag = r'$^{\dagger}$' if r['sensitive'] else ''
            same = c is not None and rel is not None and c['b'] == 'released'
            lab = f"\\quad ${thr:.2f}$" + (' (= ' + tag + ')' if same else '')
            L.append(line(lab, c['K'] if c else 1, c['d'] if c else None, flag))
    L += [r'    \bottomrule', r'  \end{tabular}', r'\end{table}']
    open(os.path.join(REP, 'tables', fn), 'w').write('\n'.join(L) + '\n'); print('wrote tables/' + fn)

table('incl', 'binreview_incl.tex',
      r"Inclusive binnings under the two criteria, FHC, RHC and combined. For each observable: the released "
      r"binning and the binning the rule of \S\ref{sec:binreview} selects under the $0.68$ and the $0.50$ "
      r"criterion. ``worst diag.'' is the smallest column-normalised migration diagonal over the true bins; "
      r"``trace $A_C$'' the degrees of freedom for signal; ``factors $\geq0.5$'' the number of Wiener filter "
      r"factors of at least $0.5$, the shape components the data determine with signal-to-noise of at least one; "
      r"the last column the largest per-bin total uncertainty. $^{\dagger}$ The choice depends on the details of "
      r"the rule (\S\ref{sec:binreview}). Fake data.", 'tab:binreview_incl')
table('1p', 'binreview_1p.tex',
      r"Proton-tagged binnings under the two criteria, as Table~\ref{tab:binreview_incl}, with an event floor "
      r"of $50$ selected FHC signal events per bin. ``study'' marks the proposed $\theta_p$ and $\theta_{\pi p}$ "
      r"binnings, which are not in the released result set. $^{\ast}$ Below the event floor: the released "
      r"$W_{\rm had}$ and $W_{\pi p}$ hold $3$ and $45$ events in their sparsest bin, the study $\theta_p$ and "
      r"$\theta_{\pi p}$ $21$ and $42$. $^{\dagger}$ The choice depends on the details of the rule. Fake data.",
      'tab:binreview_1p')

# edges of every selected binning, as written into the bin configurations (three decimals)
def ef(e): return r'$\pi$' if abs(e - PI) < 1e-3 else f'{round(e, 3):g}'

E = [r'% generated by report/tools/binreview_tables.py -- do not edit', r'\begin{table}[H]', r'  \centering',
     r"  \caption{Bin edges of the binnings selected under each criterion (Tables~\ref{tab:binreview_incl} and "
     r"\ref{tab:binreview_1p}); the last bin of $p_\mu$, $p_\pi$, $W_{\rm had}$, $W_{\pi p}$, $\delta p_T$ and $p_n$ is open "
     r"above, as in the release. ``released'' where the rule keeps the released edges.}",
     r'  \label{tab:binreview_edges}', r'  \footnotesize', r'  \begin{tabular}{L{2.6cm}cp{10.4cm}}', r'    \toprule',
     r'    Observable & criterion & edges \\', r'    \midrule']
for fam in ORDER:
    for obs in ORDER[fam]:
        v = res[(fam, obs)]; lab = LAB[obs] + (f' ({UNIT[obs]})' if obs in UNIT else '')
        for j, thr in enumerate(CRIT):
            c = v[str(thr)]['c']
            if c is None: txt = 'one bin'
            elif c['b'] == 'released': txt = 'released'
            else: txt = ', '.join(ef(e) for e in edges_of(fam, obs, c['b'], c['d']))
            E.append(f"    {lab if j == 0 else ''} & ${thr:.2f}$ & {txt} \\\\")
    if fam == 'incl': E.append(r'    \midrule')
E += [r'    \bottomrule', r'  \end{tabular}', r'\end{table}']
open(os.path.join(REP, 'tables', 'binreview_edges.tex'), 'w').write('\n'.join(E) + '\n'); print('wrote tables/binreview_edges.tex')

for (fam, obs), v in res.items():
    s = []
    for thr in CRIT:
        r = v[str(thr)]; c = r['c']
        s.append(f"{thr}: {c['b'] if c else 'one bin'} (K{c['K'] if c else 1}, finest {r['kmax']}){' SENSITIVE ' + '/'.join(r['alts']) if r['sensitive'] else ''}")
    print(f"{fam:4s} {obs:8s} " + ' | '.join(s))
