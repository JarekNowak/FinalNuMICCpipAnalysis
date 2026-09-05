# ensemble_stat_pulls.py OBS -- pull mean/width and 68/95% coverage of the fake-data ensemble against
# (a) the full covariance (sidecar errors) and (b) the data-statistical covariance alone
# (mat_table_cov_DataStats.txt from slurm/ens_statcov.sh), plus the integral offset of the ensemble
# mean from the fixed central-value reference (h_genie_tune, smeared by each member's own A_C).
import sys, math, ROOT
ROOT.gROOT.SetBatch(True)
obs=sys.argv[1]; P='/data/uboone/processed/ens/'
def cov(path):
    n=None; C={}
    for l in open(path):
        p=l.split()
        if p[0]=='numXbins': n=int(p[1])
        elif p[0] not in ('numYbins','xbin') and len(p)==3: C[(int(p[0]),int(p[1]))]=float(p[2])
    return n,C
pull_full=[]; pull_stat=[]; ints=[]; refs=[]; nb=None
for t in range(1,33):
    f=ROOT.TFile.Open(P+f'closure_hists_xsec_{obs}_t{t}.root')
    if not f or f.IsZombie(): continue
    u=f.Get('h_unfolded_nuwro'); r=f.Get('h_genie_tune')
    n,Cs=cov(P+f'statcov_{obs}_t{t}/unfold_output/mat_table_cov_DataStats.txt'); _,Ct=cov(P+f'statcov_{obs}_t{t}/unfold_output/mat_table_cov_total.txt')
    if nb is None: nb=u.GetNbinsX(); pull_full=[[] for _ in range(nb)]; pull_stat=[[] for _ in range(nb)]
    iu=ir=0.
    for b in range(1,nb+1):
        w=u.GetBinWidth(b); e=u.GetBinError(b); d=u.GetBinContent(b)-r.GetBinContent(b)
        iu+=u.GetBinContent(b)*w; ir+=r.GetBinContent(b)*w
        k=e/math.sqrt(Ct[(b-1,b-1)])          # unit conversion sidecar-error / table-error
        es=math.sqrt(Cs[(b-1,b-1)])*k
        pull_full[b-1].append(d/e); pull_stat[b-1].append(d/es)
    ints.append(iu); refs.append(ir); f.Close()
def summ(P):
    allp=[x for b in P for x in b]; n=len(allp); m=sum(allp)/n; w=math.sqrt(sum((x-m)**2 for x in allp)/(n-1))
    c68=sum(1 for x in allp if abs(x)<1)/n; c95=sum(1 for x in allp if abs(x)<1.96)/n; return m,w,c68,c95
m=sum(ints)/len(ints); ref=sum(refs)/len(refs); sd=math.sqrt(sum((x-m)**2 for x in ints)/(len(ints)-1))
print(f"{obs}: members {len(ints)}, bins {nb}")
for lab,P in [('statistical only',pull_stat),('full',pull_full)]:
    mm,w,c68,c95=summ(P); print(f"  {lab:<18} pull mean {mm:+.2f} width {w:.2f} cov68 {100*c68:.1f}% cov95 {100*c95:.1f}%")
    print("    per-bin means:", ' '.join(f"{sum(b)/len(b):+.2f}" for b in P))
print(f"  integral: ensemble mean {m:.4f} vs reference {ref:.4f} -> {100*(m/ref-1):+.2f}% ({(m-ref)/(sd/math.sqrt(len(ints))):+.1f} sigma on the mean; throw-to-throw sd {100*sd/ref:.1f}%)")
