# rebuild_tables.py -- regenerate EVERY release-derived table in the note and supplement as a
# whole tabular, from the release dumps, closure sidecars, unfolder logs, released curves and
# the counting log. Locator: the first \begin{tabular} AFTER the label, which must lie before
# the enclosing \end{table} / \end{center}. (The previous patcher located tables by proximity
# and, for captions longer than 400 characters, rewrote the PRECEDING table.)
import re,sys,csv,math,glob,collections
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'
D='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/systdump'; LG='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
CLF=sys.argv[1]; CNT=sys.argv[2]
CL={(r['family'],r['config'],r['obs']):r for r in csv.DictReader(open(CLF),delimiter='\t')}
def load(f):
    d={}
    for line in open(f):
        if line.startswith('[SYSTDUMP]'): _,k,v=line.split(); d[k]=float(v)
    return d
def det(d): return d['detVar_total'] if 'detVar_total' in d else math.sqrt(sum(v*v for k,v in d.items() if k.startswith('detVar') and not k.endswith('_total')))
CF={'FHC5':'fhc5','RHCFULL':'rhcfull','COMB':'comb'}; CFS=['FHC5','RHCFULL','COMB']
INCL=['pmu','ppi2bin','costhmu','costhpi','thmupi','thetamu']; P1=['Whad','dpt2bin','dalphat2bin','dphit2bin','pn2bin','pmu','ppi2bin','costhmu','costhpi','thmupi']; P1ALL=['Wpipr']+P1
DU={}
for c in CFS:
    for o in INCL: DU[('incl',c,o)]=load(f"{D}/{CF[c]}_{'ppi' if o=='ppi2bin' else o}.dump")
    for o in P1ALL: DU[('1p',c,o)]=load(f"{D}/1p_{CF[c]}_{o}.dump")
LAB={'pmu':r'$p_\mu$','ppi2bin':r'$p_\pi$','costhmu':r'$\cos\theta_\mu$','costhpi':r'$\cos\theta_\pi$','thmupi':r'$\theta_{\mu\pi}$','thetamu':r'$\theta_\mu$',
     'Wpipr':r'$W_{\pi p}$','Whad':r'$W_\mathrm{had}$','dpt2bin':r'$\delta p_T$','dalphat2bin':r'$\delta\alpha_T$','dphit2bin':r'$\delta\phi_T$','pn2bin':r'$p_n$'}
# model chi2 from unfolder logs
MODELS=['MicroBooNE Tune','GENIE','GiBUU','NEUT','NuWro','truth']
def chi2log(fam,c,o):
    f=f"{LG}{'1p_' if fam=='1p' else ''}{CF[c]}_{'ppi' if (fam=='incl' and o=='ppi2bin') else o}.raw"; out={}
    for line in open(f):
        m=re.match(r'^(MicroBooNE Tune|GENIE|GiBUU|NEUT|NuWro|truth): .* = ([0-9.e+-]+)/(\d+) bin',line)
        if m and m.group(1) not in out: out[m.group(1)]=(float(m.group(2)),int(m.group(3)))
    return out
def locate(text,label):
    i=text.index('\\label{%s}'%label); a=text.index('\\begin{tabular}',i)
    ends=[x for x in (text.find('\\end{table}',i),text.find('\\end{center}',i)) if x>0]
    assert a<min(ends), f"{label}: tabular not inside the same environment"
    b=text.index('\\end{tabular}',a)+len('\\end{tabular}'); return a,b
def replace(text,label,body):
    a,b=locate(text,label); return text[:a]+body+text[b:]
def T(cols,header,rows): return '\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n%s\n\\bottomrule\n\\end{tabular}'%(cols,header,'\n'.join(rows))
note=open(R+'analysis_note.tex').read(); supp=open(R+'technical_supplement.tex').read()
# --- tab:fhc / rhc / comb
for lab,c in [('tab:fhc','FHC5'),('tab:rhc','RHCFULL'),('tab:comb','COMB')]:
    rows=[]
    for o in INCL:
        d=DU[('incl',c,o)]; cl=CL[('incl',c,o)]
        rows.append(f"      {LAB[o]:<18} & ${d['sigma_int']:.3f}$ & ${d['flux_total']:.1f}\\%$ & ${det(d):.1f}\\%$ & ${d['PredTotal']:.1f}\\%$ & ${float(cl['chi2']):.2f}/{cl['ndf']}$ ($p={float(cl['pval']):.2f}$) \\\\")
    note=replace(note,lab,T('lccccc','Observable & $\\sigma_\\mathrm{int}$ & Flux & Detector & Total & Closure $\\chi^2/\\mathrm{ndf}$',rows))
# --- tab:sigint_all
rows=[f"{LAB[o]:<17} & "+' & '.join(f"${DU[('incl',c,o)]['sigma_int']:.3f}$" for c in CFS)+' \\\\' for o in INCL[:5]]
rows.append('\\midrule'); rows.append(f"{LAB['thetamu']:<17} & "+' & '.join(f"${DU[('incl',c,'thetamu')]['sigma_int']:.3f}$" for c in CFS)+' \\\\')
note=replace(note,'tab:sigint_all',T('lccc','Observable & FHC & RHC & Combined',rows))
# --- tab:chi2_incl and tab:chi2_theta
def chirow(beam,o,c,first):
    x=chi2log('incl',c,o); cells=[f"${x[m][0]:.2f}/{x[m][1]}$" for m in ['truth','MicroBooNE Tune','GENIE','GiBUU','NEUT','NuWro']]
    return f"    {beam if first else '':<4} & {LAB[o]:<22} & "+' & '.join(cells)+' \\\\'
rows=[]
for beam,c in [('FHC','FHC5'),('RHC','RHCFULL'),('comb','COMB')]:
    for k,o in enumerate(INCL): rows.append(chirow(beam,o,c,k==0))
    if c!='COMB': rows.append('    \\midrule')
note=replace(note,'tab:chi2_incl',T('llcccccc','Beam & Observable & truth & uB tune & GENIE & GiBUU & NEUT & NuWro',rows))
rows=[]
for beam,c in [('FHC','FHC5'),('RHC','RHCFULL'),('Comb','COMB')]:
    for o in ['costhmu','thetamu']:
        x=chi2log('incl',c,o); rows.append(f"{beam} & {LAB[o]:<18} & "+' & '.join(f"${x[m][0]:.2f}/{x[m][1]}$" for m in ['truth','MicroBooNE Tune','GENIE','GiBUU','NEUT','NuWro'])+' \\\\')
    if c!='COMB': rows.append('\\midrule')
note=replace(note,'tab:chi2_theta',T('llcccccc','Beam & Variable & truth & uB tune & GENIE & GiBUU & NEUT & NuWro',rows))
# --- tab:wtki_*
for lab,c in [('tab:wtki_fhc','FHC5'),('tab:wtki_rhc','RHCFULL'),('tab:wtki_comb','COMB')]:
    rows=[]
    for o in P1:
        d=DU[('1p',c,o)]; cl=CL[('1p',c,o)]
        rows.append(f"      {LAB[o]:<18} & ${d['sigma_int']:.3f}$ & ${float(cl['unf_over_truth']):.2f}$ & ${float(cl['chi2']):.2f}/{cl['ndf']}$ \\\\")
        if o=='pn2bin': rows.append('      \\midrule')
    note=replace(note,lab,T('lccc','Observable & $\\sigma_\\mathrm{int}$ & unf./truth & $\\chi^2/\\mathrm{ndf}$',rows))
# --- systbreak tables
ROWS=[('\\textbf{Prediction total}','PredTotal'),('Cross section (GENIE)','xsec_total'),('Flux (PPFX)','flux_total'),('Detector','DET'),('Reinteraction','reint'),('MC stat','MCstats'),('EXT stat','EXTstats'),('Data stat','DataStats'),('POT $+$ targets','POTT'),('\\textbf{Total (incl.\\ data stat)}','total')]
def sb(c):
    rows=[]
    for lab,key in ROWS:
        vals=[]
        for o in INCL:
            d=DU[('incl',c,o)]; v=det(d) if key=='DET' else (math.sqrt(d['POT']**2+d['numTargets']**2) if key=='POTT' else d[key]); vals.append(f"{v:.1f}")
        if lab.startswith('\\textbf'): vals=[f"\\textbf{{{v}}}" for v in vals]
        rows.append(lab+' & '+' & '.join(vals)+' \\\\')
    return T('lcccccc','Source & '+' & '.join(LAB[o] for o in INCL),rows)
note=replace(note,'tab:systbreak',sb('FHC5'))
for lab,c in [('tab:systbreak_fhc','FHC5'),('tab:systbreak_rhc','RHCFULL'),('tab:systbreak_comb','COMB')]: supp=replace(supp,lab,sb(c))
# --- tab:systematics (ranges)
def rng(fam,key,bold=False):
    vs=[(det(DU[(fam,c,o)]) if key=='DET' else DU[(fam,c,o)][key]) for c in CFS for o in (INCL if fam=='incl' else P1ALL)]
    return (f"$\\mathbf{{{min(vs):.1f}}}$--$\\mathbf{{{max(vs):.1f}}}$" if bold else f"${min(vs):.1f}$--${max(vs):.1f}$")
SR=[('Flux (PPFX multisims)','\\texttt{weightsFlux}','flux_total',0),('Detector response','Dedicated samples','DET',0),('Cross-section model','\\texttt{weightsGenie}','xsec_total',0),('Hadron re-interaction','\\texttt{weightsReint}','reint',0),('POT counting','Beam toroids','POT',0),('Target count','FV geometry','numTargets',0),None,('MC statistics','universe spread','MCstats',0),('EXT statistics','beam-off sample','EXTstats',0),('Data statistics','thrown fake data','DataStats',0),None,('\\textbf{Prediction total}','quadrature sum','PredTotal',1),('\\textbf{Total}','incl.\\ data stats','total',1)]
rows=['      \\midrule' if r is None else f"      {r[0]:<24}& {r[1]:<22}& {rng('incl',r[2],r[3])} & {rng('1p',r[2],r[3])} \\\\" for r in SR]
note=replace(note,'tab:systematics',T('L{3.4cm} L{3.4cm} c c','Source & Branch / method & Inclusive (\\%) & Proton-tagged (\\%)',rows))
# --- tab:ppi_partial
def part(c):
    rows=[l.split('\t') for l in open(R+f'data_release/curves_incl_{c}_ppi2bin.tsv') if not l.startswith('#') and not l.startswith('bin')]
    return [(float(r[3])*float(r[4]), float(r[3])*float(r[5])) for r in rows]
pp={c:part(c) for c in CFS}
rows=[r"$\sigma(0.175<p_\pi<0.205\GeVc)$ & "+' & '.join(f"${pp[c][0][0]:.3f}\\pm{pp[c][0][1]:.3f}$" for c in CFS)+' \\\\',
      r"$\sigma(p_\pi>0.205\GeVc)$        & "+' & '.join(f"${pp[c][1][0]:.3f}\\pm{pp[c][1][1]:.3f}$" for c in CFS)+' \\\\']
note=replace(note,'tab:ppi_partial',T('lccc','Region & FHC & RHC & Combined',rows))
# --- tab:cutcount + supplement tab:total_xsec from the counting log
vals=collections.OrderedDict(); cur=None; gens={}
for line in open(CNT):
    m=re.match(r'=== (CC1mu1piXp|CC1mu1pi1p) (FHC|RHC|Combined) ===',line)
    if m: cur=(m.group(1),m.group(2)); vals[cur]={}
    elif cur:
        for k,pat in [('eff',r'eff=([\d.]+)'),('bkg',r'N_bkg=([\d.]+)'),('nsel',r'N_sel=([\d.]+)'),('sig',r'sigma_count = ([\d.]+)'),('stat',r'data-stat ([\d.]+)%')]:
            mm=re.search(pat,line)
            if mm: vals[cur][k]=float(mm.group(1))
        if 'generators' in line: gens[cur]={g:float(v) for g,v in re.findall(r'(genie|gibuu|neut|nuwro)=([\d.]+)',line)}
SY={('CC1mu1piXp','FHC'):DU[('incl','FHC5','pmu')]['PredTotal'],('CC1mu1piXp','RHC'):DU[('incl','RHCFULL','pmu')]['PredTotal'],('CC1mu1piXp','Combined'):DU[('incl','COMB','pmu')]['PredTotal'],
    ('CC1mu1pi1p','FHC'):DU[('1p','FHC5','pmu')]['PredTotal'],('CC1mu1pi1p','RHC'):DU[('1p','RHCFULL','pmu')]['PredTotal'],('CC1mu1pi1p','Combined'):DU[('1p','COMB','pmu')]['PredTotal']}
def fmt(x):
    if x>=1000: return f"{int(round(x)):,}".replace(',', '\\,')
    return f"{x:.1f}" if x<10 else f"{int(round(x))}"
rows=['    \\multicolumn{10}{l}{\\emph{Inclusive $\\mathrm{CC}\\,1\\mu1\\pi Xp$}} \\\\']
TOT={}
for fam,title in [('CC1mu1piXp',None),('CC1mu1pi1p','    \\multicolumn{10}{l}{\\emph{Proton-tagged $\\mathrm{CC}\\,1\\mu1\\pi1p$ ($p_p>0.3$~GeV/$c$)}} \\\\')]:
    if title: rows.append('    \\midrule'); rows.append(title)
    for cfg in ['FHC','RHC','Combined']:
        v=vals[(fam,cfg)]; sy=SY[(fam,cfg)]; tot=math.sqrt(sy*sy+v['stat']**2); TOT[(fam,cfg)]=tot
        rows.append(f"    {cfg:<8} & ${fmt(v['nsel'])}$ & ${fmt(v['bkg'])}$ & ${100*v['eff']:.1f}$ & ${v['sig']:.3f}\\pm{v['sig']*tot/100:.3f}$ & ${v['stat']:.1f}$ & ${sy:.1f}$ & ${tot:.1f}$ & ${v['sig']:.3f}$ & $1.000$ \\\\")
note=replace(note,'tab:cutcount',T('lrrccccccc','Config & $N_\\mathrm{sel}$ & $N_\\mathrm{bkg}$ & $\\varepsilon$ [\\%] & $\\sigma$ & stat [\\%] & syst [\\%] & total [\\%] & truth & $\\sigma$/truth',rows))
rows=[]
for cfg in ['FHC','RHC','Combined']:
    v=vals[('CC1mu1piXp',cfg)]; g=gens[('CC1mu1piXp',cfg)]; err=v['sig']*TOT[('CC1mu1piXp',cfg)]/100
    cells=[f"${g[k]:.2f}\\,({abs(v['sig']-g[k])/err:.1f}\\sigma)$" for k in ['genie','gibuu','neut','nuwro']]
    rows.append(f"    {cfg:<8} & ${v['sig']:.2f} \\pm {err:.2f}$ & "+' & '.join(cells)+' \\\\')
supp=replace(supp,'tab:total_xsec',T('lccccc','Config & Measurement & GENIE & GiBUU & NEUT & NuWro',rows))
# --- supplement W_pipr six-bin closure
rows=[f"{n:<8} & ${DU[('1p',c,'Wpipr')]['sigma_int']:.3f}$ & ${float(CL[('1p',c,'Wpipr')]['unf_over_truth']):.2f}$ & ${float(CL[('1p',c,'Wpipr')]['chi2']):.2f}/6$ \\\\" for n,c in [('FHC','FHC5'),('RHC','RHCFULL'),('Combined','COMB')]]
supp=replace(supp,'tab:wpipr_sixbin_closure',T('lccc','Config & $\\sigma_\\mathrm{int}$ & unf./truth & $\\chi^2/\\mathrm{ndf}$',rows))
open(R+'analysis_note.tex','w').write(note); open(R+'technical_supplement.tex','w').write(supp)
sig=[DU[('incl',c,o)]['sigma_int'] for c in CFS for o in INCL]; print("rebuilt; incl sigma_int range %.3f-%.3f; cut-and-count:"%(min(sig),max(sig)), {k:round(v['sig'],3) for k,v in vals.items()})
