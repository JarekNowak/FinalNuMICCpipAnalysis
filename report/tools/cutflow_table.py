# cutflow_table.py -- regenerate the note's inclusive cut-flow table (tab:cutflow) from
# logs/fdfix/cutflow_incl_{fhc,rhc,comb}.log (macros/cutflow_yields.C on the cf/ reprocess, full CV
# weight, dirt normalisation in the weight), and print the final-stage purity, efficiency and EXT
# share the note quotes. Companion of cutflow_1p_table.py.
import re
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; L='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
def fmt(x):
    if x>=1000: return f"{int(round(x)):,}".replace(',', '\\,')
    return f"{x:.1f}" if x<10 else f"{int(round(x))}"
blocks={'fhc':'FHC (Runs 1,2,4,5; $7.766\\times10^{20}$~POT)','rhc':'RHC (Runs 1,2,3,4a,4b,4c; $1.108\\times10^{21}$~POT)','comb':'Combined FHC$+$RHC ($1.885\\times10^{21}$~POT)'}
STAGES='None InFiducialVol Topological MuonCandidate ContainedPion MuonIn3Planes PionIn3Planes ShowerCut OpeningAngle Nonprotons'.split()
rows=[]; derived={}
for m,title in blocks.items():
    lines=[l for l in open(L+f'cutflow_incl_{m}.log') if re.match(r'^(%s)\s+[0-9.]+\s'%'|'.join(STAGES),l)]
    assert len(lines)==10, (m,len(lines))
    rows.append(f'    \\multicolumn{{8}}{{l}}{{\\textit{{{title}}}}} \\\\')
    gen=None; st={}
    for l in lines:
        c,sig,bkg,ext,dirt,pred=l.split(); sig,bkg,ext,dirt,pred=map(float,(sig,bkg,ext,dirt,pred))
        if gen is None: gen=sig
        eff=100*sig/gen; pur=100*sig/pred; st[c]=(sig,bkg,ext,dirt,pred,eff,pur)
        cells=[fmt(sig),fmt(bkg),fmt(ext),fmt(dirt),fmt(pred),f'{eff:.1f}',f'{pur:.1f}']
        if c=='Nonprotons': rows.append('    \\textbf{Nonprotons} & '+' & '.join(f'\\textbf{{{x}}}' for x in cells)+' \\\\')
        else: rows.append(f'    {c} & '+' & '.join(cells)+' \\\\')
    if m!='comb': rows.append('    \\midrule')
    s=st['Nonprotons']; derived[m]=dict(sig=s[0],bkg=s[1],ext=s[2],dirt=s[3],pred=s[4],eff=s[5],pur=s[6],ext_share=100*s[2]/s[4])
body='\\begin{tabular}{lrrrrrrr}\n    \\toprule\n    Cut & Signal & Bkg MC & EXT & Dirt & Pred. & Eff. [\\%] & Pur. [\\%] \\\\\n    \\midrule\n'+'\n'.join(rows)+'\n    \\bottomrule\n  \\end{tabular}'
t=open(R+'analysis_note.tex').read()
i=t.index('\\label{tab:cutflow}'); a=t.index('\\begin{tabular}',i); b=t.index('\\end{tabular}',a)+len('\\end{tabular}')
assert '\\end{table}' in t[b:b+400]
t=t[:a]+body+t[b:]; open(R+'analysis_note.tex','w').write(t)
for m,d in derived.items(): print(m, {k:round(v,1) for k,v in d.items()})
