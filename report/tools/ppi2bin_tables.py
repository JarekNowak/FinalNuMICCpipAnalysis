# ppi2bin_tables.py -- regenerate the supplement's p_pi tables: tab:ppi2bin (inclusive, three regions since
# 2026-09-26) and tab:ppi2bin_1p (proton-tagged, two regions)
# from the release: unfolded values with total uncertainty and A_C-smeared truth from the closure
# sidecars, tune and generator curves from curves_*.tsv, sigma_int / chi2 / p from current_results.
import csv, ROOT
ROOT.gROOT.SetBatch(True)
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; P='/data/uboone/processed/'
CF=[('FHC','FHC5','fhc5'),('RHC','RHCFULL','rhcfull'),('Combined','COMB','comb')]
res={(r['family'],r['config'],r['observable']):r for f in ('current_results.tsv','current_results_1p_new.tsv') for r in csv.DictReader(open(R+'results/'+f),delimiter='\t')}
def side(fam,T,key='ppi2bin',nb=2):
    f=ROOT.TFile.Open(P+f"closure_hists_xsec_{'ccpi1p_' if fam=='1p' else ''}{T}_{key}.root"); u=f.Get('h_unfolded_nuwro'); t=f.Get('h_fakedata_truth')
    out=[(u.GetBinContent(b),u.GetBinError(b),t.GetBinContent(b)) for b in range(1,nb+1)]; f.Close(); return out
def curves(fam,T,key='ppi2bin'):
    rows=[l.split('\t') for l in open(R+f'data_release/curves_{fam}_{T}_{key}.tsv') if l[0].isdigit()]
    hdr=[l for l in open(R+f'data_release/curves_{fam}_{T}_{key}.tsv') if l.startswith('bin')][0].split()
    return {c:[float(r[hdr.index(c)]) for r in rows] for c in hdr if c.endswith('_smeared')}
t=open(R+'notes/technical_supplement.tex').read()
def replace_tab(label,body):
    global t; i=t.index('\\label{%s}'%label); a=t.index('\\begin{tabular}',i); b=t.index('\\end{tabular}',a)+len('\\end{tabular}'); t=t[:a]+body+t[b:]
# inclusive: three regions since 2026-09-26 (binnings of the 0.50 migration criterion; files keyed ppi3bin)
NB=3
S={T:side('incl',T,'ppi3bin',NB) for _,T,_ in CF}; C={T:curves('incl',T,'ppi3bin') for _,T,_ in CF}
def cells(T,fn): return ' & '.join(fn(T,b) for b in range(NB))
rows=[]
rows.append('unfolded data & '+' & '.join(cells(T,lambda T,b:f"${S[T][b][0]:.2f}$")+f" & ${float(res[('incl',c,'ppi3bin')]['sigma_int']):.3f}$" for _,T,c in CF)+' \\\\')
rows.append('uncertainty   & '+' & '.join(cells(T,lambda T,b:f"$\\pm{S[T][b][1]:.2f}$")+' & ' for _,T,c in CF)+' \\\\')
rows.append('$A_C$ truth   & '+' & '.join(cells(T,lambda T,b:f"${S[T][b][2]:.2f}$")+' & ' for _,T,c in CF)+' \\\\')
rows.append('uB tune       & '+' & '.join(cells(T,lambda T,b:f"${C[T]['tune_smeared'][b]:.2f}$")+' & ' for _,T,c in CF)+' \\\\\\midrule')
for g in ['GENIE','GiBUU','NEUT','NuWro']:
    rows.append(f'{g:<5} & '+' & '.join(cells(T,lambda T,b:f"${C[T][g+'_smeared'][b]:.2f}$")+' & ' for _,T,c in CF)+' \\\\')
rows.append('\\midrule\n$\\chi^2/\\mathrm{ndf}$ & '+' & '.join(f"\\multicolumn{{{NB+1}}}{{{'c|' if k<2 else 'c'}}}{{${float(res[('incl',c,'ppi3bin')]['chi2_truth']):.2f}/{NB}$ ($p={float(res[('incl',c,'ppi3bin')]['p_truth']):.2f}$)}}" for k,(_,T,c) in enumerate(CF))+' \\\\')
hb=' & '.join([f'bin {b+1}' for b in range(NB)]+['$\\sigma_\\mathrm{int}$'])
body=('\\begin{tabular}{l'+'|'.join(['c'*(NB+1)]*3)+'}\n\\toprule\n & '+' & '.join(f"\\multicolumn{{{NB+1}}}{{{'c|' if k<2 else 'c'}}}{{{lab}}}" for k,(lab,_,_) in enumerate(CF))
      +' \\\\\n & '+' & '.join([hb]*3)+' \\\\\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}')
replace_tab('tab:ppi2bin',body)
# proton-tagged
S1={T:side('1p',T) for _,T,_ in CF}; C1={T:curves('1p',T) for _,T,_ in CF}
rows=[]
rows.append('unfolded data & '+' & '.join(f"${S1[T][0][0]:.3f}$ & ${S1[T][1][0]:.3f}$" for _,T,c in CF)+' \\\\')
rows.append('uncertainty   & '+' & '.join(f"$\\pm{S1[T][0][1]:.3f}$ & $\\pm{S1[T][1][1]:.3f}$" for _,T,c in CF)+' \\\\')
rows.append('$A_C$ truth   & '+' & '.join(f"${S1[T][0][2]:.3f}$ & ${S1[T][1][2]:.3f}$" for _,T,c in CF)+' \\\\')
rows.append('uB tune       & '+' & '.join(f"${C1[T]['tune_smeared'][0]:.3f}$ & ${C1[T]['tune_smeared'][1]:.3f}$" for _,T,c in CF)+' \\\\\\midrule')
rows.append('unfolded/truth & '+' & '.join(f"${S1[T][0][0]/S1[T][0][2]:.2f}$ & ${S1[T][1][0]/S1[T][1][2]:.2f}$" for _,T,c in CF)+' \\\\\\midrule')
rows.append('$\\sigma_\\mathrm{int}$ & '+' & '.join(f"\\multicolumn{{2}}{{{'c|' if k<2 else 'c'}}}{{${float(res[('1p',c,'ppi2bin')]['sigma_int']):.3f}$}}" for k,(_,T,c) in enumerate(CF))+' \\\\')
rows.append('$\\chi^2/\\mathrm{ndf}$ & '+' & '.join(f"\\multicolumn{{2}}{{{'c|' if k<2 else 'c'}}}{{${float(res[('1p',c,'ppi2bin')]['chi2_truth']):.2f}/2$ ($p={float(res[('1p',c,'ppi2bin')]['p_truth']):.2f}$)}}" for k,(_,T,c) in enumerate(CF))+' \\\\')
body='\\begin{tabular}{lcc|cc|cc}\n\\toprule\n & \\multicolumn{2}{c|}{FHC} & \\multicolumn{2}{c|}{RHC} & \\multicolumn{2}{c}{Combined} \\\\\n & bin 1 & bin 2 & bin 1 & bin 2 & bin 1 & bin 2 \\\\\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}'
replace_tab('tab:ppi2bin_1p',body)
open(R+'notes/technical_supplement.tex','w').write(t)
print("incl:",{T:(round(S[T][0][0],2),round(S[T][1][0],2),round(S[T][2][0],2),res[('incl',c,'ppi3bin')]['sigma_int'],res[('incl',c,'ppi3bin')]['chi2_truth']) for _,T,c in CF})
print("1p ratios:",{T:(round(S1[T][0][0]/S1[T][0][2],2),round(S1[T][1][0]/S1[T][1][2],2),res[('1p',c,'ppi2bin')]['chi2_truth']) for _,T,c in CF})
