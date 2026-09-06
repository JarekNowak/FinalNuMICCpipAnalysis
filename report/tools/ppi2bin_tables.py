# ppi2bin_tables.py -- regenerate the supplement's two-bin p_pi tables (tab:ppi2bin, tab:ppi2bin_1p)
# from the release: unfolded values with total uncertainty and A_C-smeared truth from the closure
# sidecars, tune and generator curves from curves_*.tsv, sigma_int / chi2 / p from current_results.
import csv, ROOT
ROOT.gROOT.SetBatch(True)
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; P='/data/uboone/processed/'
CF=[('FHC','FHC5','fhc5'),('RHC','RHCFULL','rhcfull'),('Combined','COMB','comb')]
res={(r['family'],r['config'],r['observable']):r for f in ('current_results.tsv','current_results_1p_new.tsv') for r in csv.DictReader(open(R+f),delimiter='\t')}
def side(fam,T):
    f=ROOT.TFile.Open(P+f"closure_hists_xsec_{'ccpi1p_' if fam=='1p' else ''}{T}_ppi2bin.root"); u=f.Get('h_unfolded_nuwro'); t=f.Get('h_fakedata_truth')
    out=[(u.GetBinContent(b),u.GetBinError(b),t.GetBinContent(b)) for b in (1,2)]; f.Close(); return out
def curves(fam,T):
    rows=[l.split('\t') for l in open(R+f'data_release/curves_{fam}_{T}_ppi2bin.tsv') if l[0].isdigit()]
    hdr=[l for l in open(R+f'data_release/curves_{fam}_{T}_ppi2bin.tsv') if l.startswith('bin')][0].split()
    return {c:[float(r[hdr.index(c)]) for r in rows] for c in hdr if c.endswith('_smeared')}
t=open(R+'technical_supplement.tex').read()
def replace_tab(label,body):
    global t; i=t.index('\\label{%s}'%label); a=t.index('\\begin{tabular}',i); b=t.index('\\end{tabular}',a)+len('\\end{tabular}'); t=t[:a]+body+t[b:]
# inclusive
S={T:side('incl',T) for _,T,_ in CF}; C={T:curves('incl',T) for _,T,_ in CF}
def row(lab,vals): return f"{lab} & "+' & '.join(vals)+' \\\\'
def trip(fn): return [x for _,T,_ in CF for x in (fn(T,0),fn(T,1),'')]
rows=[row('unfolded data',[x for _,T,_ in CF for x in (f"${S[T][0][0]:.2f}$",f"${S[T][1][0]:.2f}$",f"${float(res[('incl',_c,'ppi2bin')]['sigma_int']):.3f}$")]) if False else None]
rows=[]
rows.append('unfolded data & '+' & '.join(f"${S[T][0][0]:.2f}$ & ${S[T][1][0]:.2f}$ & ${float(res[('incl',c,'ppi2bin')]['sigma_int']):.3f}$" for _,T,c in CF)+' \\\\')
rows.append('uncertainty   & '+' & '.join(f"$\\pm{S[T][0][1]:.2f}$ & $\\pm{S[T][1][1]:.2f}$ & " for _,T,c in CF)+' \\\\')
rows.append('$A_C$ truth   & '+' & '.join(f"${S[T][0][2]:.2f}$ & ${S[T][1][2]:.2f}$ & " for _,T,c in CF)+' \\\\')
rows.append('uB tune       & '+' & '.join(f"${C[T]['tune_smeared'][0]:.2f}$ & ${C[T]['tune_smeared'][1]:.2f}$ & " for _,T,c in CF)+' \\\\\\midrule')
for g in ['GENIE','GiBUU','NEUT','NuWro']:
    rows.append(f'{g:<5} & '+' & '.join(f"${C[T][g+'_smeared'][0]:.2f}$ & ${C[T][g+'_smeared'][1]:.2f}$ & " for _,T,c in CF)+' \\\\')
rows.append('\\midrule\n$\\chi^2/\\mathrm{ndf}$ & '+' & '.join(f"\\multicolumn{{3}}{{{'c|' if k<2 else 'c'}}}{{${float(res[('incl',c,'ppi2bin')]['chi2_truth']):.2f}/2$ ($p={float(res[('incl',c,'ppi2bin')]['p_truth']):.2f}$)}}" for k,(_,T,c) in enumerate(CF))+' \\\\')
body='\\begin{tabular}{lccc|ccc|ccc}\n\\toprule\n & \\multicolumn{3}{c|}{FHC} & \\multicolumn{3}{c|}{RHC} & \\multicolumn{3}{c}{Combined} \\\\\n & bin 1 & bin 2 & $\\sigma_\\mathrm{int}$ & bin 1 & bin 2 & $\\sigma_\\mathrm{int}$\n & bin 1 & bin 2 & $\\sigma_\\mathrm{int}$ \\\\\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}'
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
open(R+'technical_supplement.tex','w').write(t)
print("incl:",{T:(round(S[T][0][0],2),round(S[T][1][0],2),res[('incl',c,'ppi2bin')]['sigma_int'],res[('incl',c,'ppi2bin')]['chi2_truth']) for _,T,c in CF})
print("1p ratios:",{T:(round(S1[T][0][0]/S1[T][0][2],2),round(S1[T][1][0]/S1[T][1][2],2),res[('1p',c,'ppi2bin')]['chi2_truth']) for _,T,c in CF})
