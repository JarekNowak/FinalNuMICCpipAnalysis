# cutflow_1p_table.py -- regenerate the supplement's proton-tagged cut-flow table (tab:cutflow_1p)
# from logs/fdfix/cutflow_1p_{fhc,rhc,comb}.log (macros/cutflow_yields_1p.C on the cf_1p reprocess),
# and print the derived quantities the note quotes (final purity, efficiency, background rejection
# and signal retention of the proton-identification stage).
import re
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; L='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
def fmt(x):
    if x>=1000: return f"{int(round(x)):,}".replace(',', '\\,')
    return f"{x:.1f}" if x<10 else f"{int(round(x))}"
blocks={'fhc':'FHC (Runs 1,2,4,5; $8.857\\times10^{20}$~POT)','rhc':'RHC (Runs 1,2,3,4a,4b,4c; $1.108\\times10^{21}$~POT)','comb':'Combined FHC$+$RHC ($1.994\\times10^{21}$~POT)'}
rows=[]; derived={}
for m,title in blocks.items():
    lines=[l for l in open(L+f'cutflow_1p_{m}.log') if re.match(r'^(None|InFiducialVol|Topological|MuonCandidate|ContainedPion|MuonIn3Planes|PionIn3Planes|ShowerCut|OpeningAngle|Nonprotons|IdentProton)\s',l)]
    assert len(lines)==11, (m,len(lines))
    rows.append(f'    \\multicolumn{{8}}{{l}}{{\\textit{{{title}}}}} \\\\')
    gen=None; st={}
    for l in lines:
        c,sig,bkg,ext,dirt,pred=l.split(); sig,bkg,ext,dirt,pred=map(float,(sig,bkg,ext,dirt,pred))
        if gen is None: gen=sig
        eff=100*sig/gen; pur=100*sig/pred; st[c]=(sig,bkg,ext,dirt,pred,eff,pur)
        cells=[fmt(sig),fmt(bkg),fmt(ext),fmt(dirt),fmt(pred),f'{eff:.1f}',f'{pur:.1f}' if pur>=0.1 else f'{pur:.2f}']
        if c=='IdentProton': rows.append('    \\textbf{IdentProton} & '+' & '.join(f'\\textbf{{{x}}}' for x in cells)+' \\\\')
        else: rows.append(f'    {c:<16} & '+' & '.join(cells)+' \\\\')
    s=st['IdentProton']; p=st['Nonprotons']
    derived[m]=dict(purity=s[6],eff=s[5],pur_nonp=p[6],bkg_rej=100*(1-s[1]/p[1]),sig_ret=100*s[0]/p[0],ext_ratio=p[2]/s[2] if s[2]>0 else float('nan'),nsel=s[4],sig=s[0])
body='\\begin{tabular}{lrrrrrrr}\n    \\toprule\n    Cut & Signal & Bkg MC & EXT & Dirt & Pred. & Eff. [\\%] & Pur. [\\%] \\\\\n    \\midrule\n'+'\n'.join(rows)+'\n    \\bottomrule\n  \\end{tabular}'
t=open(R+'technical_supplement.tex').read()
i=t.index('\\label{tab:cutflow_1p}'); a=t.index('\\begin{tabular}',i); b=t.index('\\end{tabular}',a)+len('\\end{tabular}')
assert '\\end{table}' in t[b:b+400]
t=t[:a]+body+t[b:]; open(R+'technical_supplement.tex','w').write(t)
for m,d in derived.items(): print(m, {k:round(v,1) for k,v in d.items()})
