# index_extractions.py -- regenerate data_release/index_extractions.tsv from the current release
# (current_results.tsv + current_results_1p_new.tsv, curves_*.tsv, cov/ directories). One row per
# released differential extraction; status from the note's scope matrix.
import csv, os, glob
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'; D=R+'data_release/'
CF={'fhc5':'FHC5','rhcfull':'RHCFULL','comb':'COMB'}
STATUS={'Wpipr':'withdrawn (differential shape; integral only)','ppi2bin':'primary (two-bin, coarse)'}
rows=[]
for fn,fam in [('current_results.tsv','incl'),('current_results_1p_new.tsv','1p')]:
    for r in csv.DictReader(open(R+fn),delimiter='\t'):
        if r['family']!=fam: continue
        tag=f"{fam}_{CF[r['config']]}_{r['observable']}"
        cur=D+f'curves_{tag}.tsv'; ac=D+f'A_C_{tag}.tsv'; cov=D+f'cov/{tag}'
        assert os.path.exists(cur) and os.path.exists(ac) and os.path.isdir(cov), tag
        nb=sum(1 for l in open(cur) if l[0].isdigit()); nc=len(glob.glob(cov+'/*.txt'))
        st=STATUS.get(r['observable'],'primary' if fam=='incl' else 'secondary')
        rows.append([tag,fam,CF[r['config']],r['observable'],str(nb),f"{float(r['sigma_int']):.4e}",r['PredTotal_pct'],r['detVar_pct'],f"{float(r['chi2_truth']):.3f}",f"{float(r['p_truth']):.3f}",os.path.basename(ac),os.path.basename(cur),str(nc),st])
rows.sort()
with open(D+'index_extractions.tsv','w') as o:
    o.write('# Definitive index of every released differential extraction (release 2026-09-06).\n# One row per extraction; A_C / curves / cov files are all present for each. total_pct is the\n# prediction-total fractional uncertainty (bin-averaged, excluding data statistics); chi2/p are the\n# closure against the realised fake-data truth with the complete configured covariance. The six one-bin totals are\n# indexed separately in total_xsec.tsv.\n')
    o.write('tag\tfamily\tconfig\tobservable\tbins\tsigma_int\ttotal_pct\tdetVar_pct\tchi2\tp\tA_C\tcurves\tcov_components\tstatus\n')
    for r in rows: o.write('\t'.join(r)+'\n')
print(len(rows),'rows')
