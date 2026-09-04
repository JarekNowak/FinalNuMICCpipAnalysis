import csv,math,sys
D='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/logs/systdump'; R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'
cl={(r['family'],r['config'],r['obs']):r for r in csv.DictReader(open(sys.argv[1]),delimiter='\t')}
def load(f):
    d={}
    for line in open(f):
        if line.startswith('[SYSTDUMP]'): _,k,v=line.split(); d[k]=float(v)
    return d
def det(d): return d['detVar_total'] if 'detVar_total' in d else math.sqrt(sum(v*v for k,v in d.items() if k.startswith('detVar') and not k.endswith('_total')))
CF={'FHC5':'fhc5','RHCFULL':'rhcfull','COMB':'comb'}
hdr='family\tobservable\tconfig\tsigma_int\tPredTotal_pct\tdetVar_pct\tflux_pct\txsec_pct\tchi2_truth\tp_truth\n'
out=open(R+'current_results.tsv','w'); out1=open(R+'current_results_1p_new.tsv','w'); out.write(hdr); out1.write(hdr)
for fam,obs,pref in [('incl',['pmu','ppi2bin','costhmu','costhpi','thmupi','thetamu'],''),('1p',['Whad','Wpipr','costhmu','costhpi','dalphat2bin','dphit2bin','dpt2bin','pmu','pn2bin','ppi2bin','thmupi'],'1p_')]:
    for cfg in ['FHC5','RHCFULL','COMB']:
        for o in obs:
            lo='ppi' if (fam=='incl' and o=='ppi2bin') else o
            d=load(f'{D}/{pref}{CF[cfg]}_{lo}.dump'); c=cl[(fam,cfg,o)]
            row=f"{fam}\t{o}\t{CF[cfg]}\t{d['sigma_int']:.4e}\t{d['PredTotal']:.1f}\t{det(d):.1f}\t{d['flux_total']:.1f}\t{d['xsec_total']:.1f}\t{float(c['chi2']):.3f}\t{float(c['pval']):.3f}\n"
            out.write(row)
            if fam=='1p': out1.write(row)
out.close(); out1.close(); print("tsv ok")
