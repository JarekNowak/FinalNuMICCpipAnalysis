# total_export.py -- export the six one-bin total cross-section extractions (inclusive and
# proton-tagged x FHC/RHC/combined) to data_release/total_xsec.tsv: central value, full
# universe-propagated uncertainty and its breakdown, realised fake-data truth, CV tune, A_C.
import ROOT, re, sys
ROOT.gROOT.SetBatch(True)
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/'; P='/data/uboone/processed/'
CF={'fhc5':'FHC5','rhcfull':'RHCFULL','comb':'COMB'}
def dump(f):
    d={}
    for l in open(f):
        if l.startswith('[SYSTDUMP]'): _,k,v=l.split(); d[k]=float(v)
    return d
def chi2(f):
    for l in open(f):
        m=re.search(r'truth: .* = ([0-9.e+-]+)/1 bin',l)
        if m: return float(m.group(1))
    return float('nan')
out=open(R+'report/data_release/total_xsec.tsv','w')
out.write('# One-bin total flux-averaged cross section (10^-38 cm^2/Ar): one true bin (all signal, cos(theta_mu) in [-1,1],\n'
          '# no overflow) and one reco bin, Wiener filter off (plain inversion, A_C = 1). Percentages are fractional\n'
          '# uncertainties of the prediction; sigma is the extraction on the Poisson-thrown fake data, truth its realised truth.\n')
cols=['family','config','sigma','err_total','stat_pct','syst_pct','total_pct','flux_pct','detVar_pct','xsec_pct','reint_pct','POT_pct','numTargets_pct','MCstats_pct','EXTstats_pct','truth_realised','tune_cv','sigma_over_truth','chi2_truth','A_C']
out.write('\t'.join(cols)+'\n')
for fam,pre,fp in [('incl','',''),('1p','1p_','ccpi1p_')]:
    for c in ['fhc5','rhcfull','comb']:
        d=dump(R+f'logs/systdump/{pre}{c}_total.dump'); f=ROOT.TFile.Open(P+f'closure_hists_xsec_{fp}{CF[c]}_total.root')
        u=f.Get('h_unfolded_nuwro'); t=f.Get('h_fakedata_truth'); g=f.Get('h_genie_tune'); a=f.Get('h_A_C')
        w=u.GetBinWidth(1); sig=u.GetBinContent(1)*w; err=u.GetBinError(1)*w; tr=t.GetBinContent(1)*w; tune=g.GetBinContent(1)*w
        assert abs(sig-d['sigma_int'])<1e-3*max(1,sig), (sig,d['sigma_int'])
        tot=(d['PredTotal']**2+d['DataStats']**2)**0.5
        row=[fam,c,f'{sig:.4f}',f'{err:.4f}',f"{d['DataStats']:.2f}",f"{d['PredTotal']:.2f}",f'{tot:.2f}',f"{d['flux_total']:.2f}",f"{d['detVar_total']:.2f}",f"{d['xsec_total']:.2f}",f"{d['reint']:.2f}",f"{d['POT']:.2f}",f"{d['numTargets']:.2f}",f"{d['MCstats']:.2f}",f"{d['EXTstats']:.2f}",f'{tr:.4f}',f'{tune:.4f}',f'{sig/tr:.4f}',f'{chi2(R+f"logs/fdfix/{pre}{c}_total.raw"):.4f}',f'{a.GetBinContent(1,1):.6f}']
        out.write('\t'.join(row)+'\n'); f.Close()
out.close(); print('wrote total_xsec.tsv')
