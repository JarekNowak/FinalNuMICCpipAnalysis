#!/usr/bin/env python3
"""check_tables.py -- refuse a note whose summary tables disagree with the release.
Compares, for every configuration and observable, the values printed in the note's result
tables (tab:fhc/rhc/comb, tab:sigint_all, tab:wtki_*, tab:chi2_incl) with
report/current_results.tsv (written from the release systematics dumps and closure
sidecars), and checks that every table sits under its own label with the column count its
header declares. Exit 1 on any mismatch. Added after a table patcher had rewritten the
PRECEDING table for captions longer than 400 characters (Draft 1.2 review).
  python3 check_tables.py
"""
import re,sys,os,csv
R=os.path.dirname(os.path.abspath(__file__))
DOCS={f:open(os.path.join(R,f)).read() for f in ('analysis_note.tex','technical_supplement.tex','proton_tagged_note.tex')}
def doc_of(label):
    hits=[t for t in DOCS.values() if '\\label{%s}'%label in t]
    assert len(hits)==1, label
    return hits[0]
res={(r['family'],r['config'],r['observable']):r for r in csv.DictReader(open(os.path.join(R,'current_results.tsv')),delimiter='\t')}
# 2026-09-26: 0.50-criterion binnings -- inclusive p_pi keyed ppi3bin, theta_mu dropped, delta phi_T keyed dphit3bin
LAB={'pmu':r'$p_\mu$','ppi3bin':r'$p_\pi$','costhmu':r'$\cos\theta_\mu$','costhpi':r'$\cos\theta_\pi$','thmupi':r'$\theta_{\mu\pi}$',
     'Whad':r'$W_\mathrm{had}$','dpt2bin':r'$\delta p_T$','dalphat2bin':r'$\delta\alpha_T$','dphit3bin':r'$\delta\phi_T$','pn2bin':r'$p_n$','pp':r'$p_p$'}
INV={v:k for k,v in LAB.items()}
# the proton-tagged tables list p_pi (two regions) under the same symbol
INV1=dict(INV); INV1[r'$p_\pi$']='ppi2bin'
bad=[]
def table(label):
    note=doc_of(label)
    i=note.index('\\label{%s}'%label); a=note.index('\\begin{tabular}',i)
    ends=[x for x in (note.find('\\end{table}',i),note.find('\\end{center}',i)) if x>0]
    if a>min(ends): bad.append(f"{label}: no tabular between label and end of environment"); return None
    body=note[a:note.index('\\end{tabular}',a)]
    spec=re.match(r'\\begin\{tabular\}\{(.*?)\}\s*$',body.split('\n')[0]).group(1); ncol=len(re.findall(r'[lcrp]\{[^}]*\}|[lcr]|[LR]\{[^}]*\}',spec.replace('|','')))
    rows=[l for l in body.split('\n') if '&' in l]
    for l in rows:
        width=l.count('&')+1+sum(int(k)-1 for k in re.findall(r'\\multicolumn\{(\d+)\}',l))   # spanned columns count
        if width!=ncol: bad.append(f"{label}: row has {width} columns, header declares {ncol}: {l.strip()[:60]}")
    return rows
def num(cell): 
    m=re.search(r'-?[0-9]+\.[0-9]+',cell); return float(m.group()) if m else None
for lab,cfg in [('tab:fhc','fhc5'),('tab:rhc','rhcfull'),('tab:comb','comb')]:
    for l in table(lab) or []:
        c=[x.strip() for x in l.split('&')]
        if c[0] in INV:
            o=INV[c[0]]; r=res[('incl',cfg,o)]
            for j,key,tol in [(1,'sigma_int',6e-4),(2,'flux_pct',0.06),(3,'detVar_pct',0.06),(4,'PredTotal_pct',0.06)]:
                v=num(c[j]); 
                if v is None or abs(v-float(r[key]))>tol: bad.append(f"{lab} {o} col{j}: note {c[j]} vs release {r[key]}")
            if abs(num(c[5])-float(r['chi2_truth']))>0.006: bad.append(f"{lab} {o} closure chi2: {c[5]} vs {r['chi2_truth']}")
for l in table('tab:sigint_all') or []:
    c=[x.strip() for x in l.split('&')]
    if c[0] in INV:
        for j,cfg in [(1,'fhc5'),(2,'rhcfull'),(3,'comb')]:
            if abs(num(c[j])-float(res[('incl',cfg,INV[c[0]])]['sigma_int']))>6e-4: bad.append(f"tab:sigint_all {INV[c[0]]} {cfg}: {c[j]}")
for l in table('tab:wtki') or []:          # 2026-09-28: FHC, RHC, combined side by side
    c=[x.strip() for x in l.split('&')]
    if c[0] in INV1:
        for k,cfg in [(1,'fhc5'),(4,'rhcfull'),(7,'comb')]:
            r=res[('1p',cfg,INV1[c[0]])]
            if abs(num(c[k])-float(r['sigma_int']))>6e-4 or abs(num(c[k+2])-float(r['chi2_truth']))>0.006: bad.append(f"tab:wtki {cfg} {INV1[c[0]]}: {c[k]} {c[k+2]} vs {r['sigma_int']} {r['chi2_truth']}")
beam=None
for l in table('tab:chi2_incl') or []:
    c=[x.strip() for x in l.split('&')]
    if c[0] in ('FHC','RHC','comb'): beam={'FHC':'fhc5','RHC':'rhcfull','comb':'comb'}[c[0]]
    if len(c)>2 and c[1] in INV and beam:
        if abs(num(c[2])-float(res[('incl',beam,INV[c[1]])]['chi2_truth']))>0.006: bad.append(f"tab:chi2_incl {beam} {INV[c[1]]}: {c[2]}")
for lab in ['tab:systbreak','tab:systematics','tab:cutcount','tab:cutflow','tab:ppi_partial','tab:ppi_partial_main','tab:ac_rowsums']: table(lab)
if bad:
    print("TABLE CONSISTENCY FAILURES (%d):"%len(bad)); [print("  "+b) for b in bad]; sys.exit(1)
print("all note tables consistent with the release (%d extractions checked)"%len(res))
