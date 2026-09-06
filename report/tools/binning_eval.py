# binning_eval.py -- compare candidate finer binnings (built by slurm_extfix.sbatch into
# processed/rebuild_extfix) with the released extraction: A_C row sums, closure chi2/p, sigma_int,
# data-statistical and prediction-total fractions, per-bin uncertainties.
import re, math, sys, ROOT
ROOT.gROOT.SetBatch(True)
RB='/data/uboone/processed/rebuild_extfix/'; P='/data/uboone/processed/'; LG='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/fdfix/'
CF={'fhc5':'FHC5','rhcfull':'RHCFULL'}
def dump(path):
    d={}; chi=None
    for l in open(path):
        if l.startswith('[SYSTDUMP]'): _,k,v=l.split(); d[k]=float(v)
        m=re.search(r'truth: .* = ([0-9.e+-]+)/(\d+) bins?, p-value = ([0-9.e+-]+)',l)
        if m and chi is None: chi=(float(m.group(1)),int(m.group(2)),float(m.group(3)))
    return d,chi
def side(path):
    f=ROOT.TFile.Open(path); u=f.Get('h_unfolded_nuwro'); a=f.Get('h_A_C'); nb=u.GetNbinsX()
    rows=[sum(a.GetBinContent(j,i) for j in range(1,nb+1)) for i in range(1,nb+1)]  # row i = smeared bin (y), sum over true bins (x)
    rel=[100*u.GetBinError(i)/abs(u.GetBinContent(i)) if u.GetBinContent(i)!=0 else float('nan') for i in range(1,nb+1)]
    f.Close(); return nb,rows,rel
for base,cands in [('costhpi',['costhpi5','costhpi6']),('thmupi',['thmupi6','thmupi7'])]:
    for cfg,T in CF.items():
        d,chi=dump(LG+f'{cfg}_{base}.raw'); nb,rows,rel=side(P+f'closure_hists_xsec_{T}_{base}.root')
        print(f"{T:8s} {base:9s} released {nb}b  sig {d['sigma_int']:.3f} Pred {d['PredTotal']:.1f}% DataStat {d['DataStats']:.1f}%  chi2 {chi[0]:.2f}/{chi[1]} p={chi[2]:.2f}  A_C rows {min(rows):.2f}-{max(rows):.2f}  bin unc {min(rel):.0f}-{max(rel):.0f}%")
        for c in cands:
            try: d,chi=dump(RB+f'unfold_ccpi_{T}_{c}.log'); nb,rows,rel=side(RB+f'closure_hists_xsec_ccpi_{T}_{c}.root')
            except Exception as e: print(f"{T:8s} {c:9s} not available ({e})"); continue
            print(f"{T:8s} {c:9s} cand     {nb}b  sig {d['sigma_int']:.3f} Pred {d['PredTotal']:.1f}% DataStat {d['DataStats']:.1f}%  chi2 {chi[0]:.2f}/{chi[1]} p={chi[2]:.2f}  A_C rows {min(rows):.2f}-{max(rows):.2f}  bin unc {min(rel):.0f}-{max(rel):.0f}%")
