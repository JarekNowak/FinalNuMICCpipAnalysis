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
note=open(os.path.join(R,'analysis_note.tex')).read()
res={(r['family'],r['config'],r['observable']):r for r in csv.DictReader(open(os.path.join(R,'current_results.tsv')),delimiter='\t')}
LAB={'pmu':r'$p_\mu$','ppi2bin':r'$p_\pi$','costhmu':r'$\cos\theta_\mu$','costhpi':r'$\cos\theta_\pi$','thmupi':r'$\theta_{\mu\pi}$','thetamu':r'$\theta_\mu$',
     'Whad':r'$W_\mathrm{had}$','dpt2bin':r'$\delta p_T$','dalphat2bin':r'$\delta\alpha_T$','dphit2bin':r'$\delta\phi_T$','pn2bin':r'$p_n$'}
INV={v:k for k,v in LAB.items()}
bad=[]
def table(label):
    i=note.index('\\label{%s}'%label); a=note.index('\\begin{tabular}',i)
    ends=[x for x in (note.find('\\end{table}',i),note.find('\\end{center}',i)) if x>0]
    if a>min(ends): bad.append(f"{label}: no tabular between label and end of environment"); return None
    body=note[a:note.index('\\end{tabular}',a)]
    spec=re.match(r'\\begin\{tabular\}\{(.*?)\}\s*$',body.split('\n')[0]).group(1); ncol=len(re.findall(r'[lcrp]\{[^}]*\}|[lcr]|[LR]\{[^}]*\}',spec.replace('|','')))
    rows=[l for l in body.split('\n') if '&' in l]
    for l in rows:
        if l.count('&')+1!=ncol: bad.append(f"{label}: row has {l.count('&')+1} columns, header declares {ncol}: {l.strip()[:60]}")
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
for lab,cfg in [('tab:wtki_fhc','fhc5'),('tab:wtki_rhc','rhcfull'),('tab:wtki_comb','comb')]:
    for l in table(lab) or []:
        c=[x.strip() for x in l.split('&')]
        if c[0] in INV:
            r=res[('1p',cfg,INV[c[0]])]
            if abs(num(c[1])-float(r['sigma_int']))>6e-4 or abs(num(c[3])-float(r['chi2_truth']))>0.006: bad.append(f"{lab} {INV[c[0]]}: {c[1]} {c[3]} vs {r['sigma_int']} {r['chi2_truth']}")
beam=None
for l in table('tab:chi2_incl') or []:
    c=[x.strip() for x in l.split('&')]
    if c[0] in ('FHC','RHC','comb'): beam={'FHC':'fhc5','RHC':'rhcfull','comb':'comb'}[c[0]]
    if len(c)>2 and c[1] in INV and beam:
        if abs(num(c[2])-float(res[('incl',beam,INV[c[1]])]['chi2_truth']))>0.006: bad.append(f"tab:chi2_incl {beam} {INV[c[1]]}: {c[2]}")
for lab in ['tab:systbreak','tab:systematics','tab:cutcount','tab:cutflow','tab:ppi_partial','tab:chi2_theta']: table(lab)
if bad:
    print("TABLE CONSISTENCY FAILURES (%d):"%len(bad)); [print("  "+b) for b in bad]; sys.exit(1)
print("all note tables consistent with the release (%d extractions checked)"%len(res))
