"""mcs_eval.py -- the MCS momentum-scale term from the DATA-SIDE variation, for any released inclusive
observable (and the one-bin totals).

The fake data are re-thrown with every MCS-measured muon momentum (reco != range, 63% of the selected
sample) scaled by 1 +- 0.05 while the simulation, the response and the selection flag stay nominal
(macros/mcs_scale_fakedata.C; extractions in rebuild_alt/ccpi_<CFG>_<obs>_mcsdata{up,dn}05). The earlier
bin-configuration scaling moved data and simulation together and therefore bounded only the
regulariser's response to a remapped variable, not the systematic; it is superseded by this.

For observables other than p_mu the scale enters only through the p_mu acceptance (the selection
cuts on the muon candidate momentum), so the effect is a small, nearly flat normalisation shift.

Per bin: half-difference D_i = (x_up - x_dn)/2 of the unfolded result, quoted as a fraction of the
nominal and in units of the released total uncertainty; the covariance C_ij = D_i D_j (fully
correlated unisim) is written in the release units (bin integrals, like cov_*.txt) as
data_release/cov/incl_<CFG>_<obs>/cov_MCSscale.txt, with cov_total_plusMCS.txt = cov_total + C.
For the one-bin totals (obs=total) there is no cov directory: the half-difference is reported against
err_total of data_release/total_xsec.tsv and nothing is written except the summary row.

    python3 report/tools/mcs_eval.py [FHC5|RHCFULL|COMB] [obs]      (default FHC5 pmu)
    python3 report/tools/mcs_eval.py all                            (every CFG x observable -> summary tsv)

RHC and combined (2026-09-24): the RHC fake data are scaled the same way
(mcs_scale_fakedata.C(s,tag,"rhc")); the combined variation scales FHC and RHC together, i.e. one
common MCS scale, which is the correlated assumption (same detector, same estimator).
"""
import numpy as np, uproot, os, sys
RB='/data/uboone/processed/rebuild_alt/'; LIVE='/data/uboone/processed/'
R='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/'
OBS_ALL=['pmu','costhmu','ppi3bin','costhpi','thmupi','total']   # 2026-09-26: theta_mu dropped, p_pi in three regions
def cov(path):
    n=None; C=None
    for l in open(path):
        p=l.split()
        if p[0]=='numXbins': n=int(p[1]); C=np.zeros((n,n))
        elif p[0] not in ('numYbins','xbin') and len(p)==3: C[int(p[0]),int(p[1])]=float(p[2])
    return C
def write(path,C,n):
    with open(path,'w') as o:
        o.write(f'numXbins {n}\nnumYbins {n}\nxbin  ybin  z\n')
        for i in range(n):
            for j in range(n): o.write(f'{i}  {j}  {C[i,j]:.17e}\n')
def total_err(CFG):
    for l in open(R+'data_release/total_xsec.tsv'):
        p=l.split('\t')
        if p[0]=='incl' and p[1]=={'FHC5':'fhc5','RHCFULL':'rhcfull','COMB':'comb'}[CFG]: return float(p[2]),float(p[3])
    raise KeyError(CFG)
CONF='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/configs/'
BINCFG={'ppi3bin':'ccpi_ppi_bin_config_3bin.txt'}
def uses_mcs_branch(obs):
    """True when the live reco bin definitions read the scaled branch. The variation rescales only
    CC1mu1piXp_candidate_muon_mom_reco and keeps the selection flag nominal, so an observable whose reco
    bins do not read it sees bit-identical fake data at any binning: its term is zero by construction."""
    t=open(CONF+BINCFG.get(obs,f'ccpi_{obs}_bin_config_opt.txt')).read()
    return 'muon_mom' in t
def evaluate(CFG,obs,quiet=False):
    up_f=RB+f'closure_hists_xsec_ccpi_{CFG}_{obs}_mcsdataup05.root'; dn_f=RB+f'closure_hists_xsec_ccpi_{CFG}_{obs}_mcsdatadn05.root'
    nom=uproot.open(LIVE+f'closure_hists_xsec_{CFG}_{obs}.root')['h_unfolded_nuwro']
    v=nom.values(); e=nom.errors(); w=np.diff(nom.axis().edges()); n=len(v)
    if obs!='total' and not uses_mcs_branch(obs):
        up=dn=v.copy()                          # zero by construction (2026-09-26, binnings of the 0.50 criterion)
        if not quiet: print(f'== {CFG} {obs}: reco binning does not read the muon momentum -> term identically zero')
    else:
        if not (os.path.exists(up_f) and os.path.exists(dn_f)): return None
        up=uproot.open(up_f)['h_unfolded_nuwro'].values(); dn=uproot.open(dn_f)['h_unfolded_nuwro'].values()
        assert len(up)==n and len(dn)==n, f'{CFG} {obs}: MCS extraction has {len(up)} bins, release {n}'
    D=(up-dn)/2
    if obs=='total':
        sig,err=total_err(CFG); e=np.array([err/w[0]])     # sidecar is per unit width; tsv is the integral
        assert abs(v[0]*w[0]-sig)<2e-4*max(1,sig), (v[0]*w[0],sig)
    I=(v*w).sum(); Iu=(up*w).sum()/I; Id=(dn*w).sum()/I; H=((up-dn)*w).sum()/2/I
    res=dict(cfg=CFG,obs=obs,nbins=n,up_int=Iu,dn_int=Id,half_int=H,max_frac=np.max(abs(D)/v),max_sig=np.max(abs(D)/e),
             frac=D/v,sig=D/e)
    if quiet: return res
    print(f'== {CFG} {obs}')
    print('bin   nominal   up/nom   dn/nom   half-diff/nom   half-diff/sigma_tot')
    for i in range(n): print(f'{i+1:>3}  {v[i]:8.4f}  {up[i]/v[i]:7.3f}  {dn[i]/v[i]:7.3f}  {D[i]/v[i]:+13.3f}  {D[i]/e[i]:+12.2f}')
    print(f'integral: nominal {I:.4f}  up {Iu:.3f}  dn {Id:.3f}  half-diff {H:+.3%}')
    print(f'mean |half-diff|/nom {np.mean(abs(D)/v):.3%}, max |half-diff|/sigma {np.max(abs(D)/e):.2f}')
    if obs!='total':
        REL=R+f'data_release/cov/incl_{CFG}_{obs}/'
        Db=D*w                      # release units: bin integrals
        C=np.outer(Db,Db); Ct=cov(REL+'cov_total.txt')
        assert Ct.shape==(n,n)
        write(REL+'cov_MCSscale.txt',C,n); write(REL+'cov_total_plusMCS.txt',Ct+C,n)
        print('total uncertainty per bin, before -> after adding the MCS term (% of bin):')
        print(' '.join(f'{100*np.sqrt(Ct[i,i])/(v[i]*w[i]):.1f}->{100*np.sqrt(Ct[i,i]+C[i,i])/(v[i]*w[i]):.1f}' for i in range(n)))
        print('wrote', REL+'cov_MCSscale.txt', 'and cov_total_plusMCS.txt')
    return res
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='all':
        rows=[]
        for CFG in ['FHC5','RHCFULL','COMB']:
            for obs in OBS_ALL:
                r=evaluate(CFG,obs)
                if r is None: print(f'== {CFG} {obs}: not extracted yet'); continue
                rows.append(r)
        out=R+'data_release/mcs_scale_2026-09-26.tsv'
        with open(out,'w') as o:
            o.write('# Data-side MCS momentum-scale (+-5%) term per released inclusive extraction: integral shift for the up/down\n'
                    '# variations, the half-difference of the integral, and the largest per-bin half-difference as a fraction of\n'
                    '# the nominal and in units of the released total uncertainty (mcs_eval.py; cov_MCSscale.txt per directory).\n')
            o.write('config\tobservable\tbins\tup_int\tdn_int\thalf_int_pct\tmax_bin_pct\tmax_bin_sigma\tper_bin_pct\n')
            for r in rows: o.write(f"{r['cfg']}\t{r['obs']}\t{r['nbins']}\t{r['up_int']:.4f}\t{r['dn_int']:.4f}\t{100*r['half_int']:+.2f}\t{100*r['max_frac']:.2f}\t{r['max_sig']:.2f}\t{' '.join(f'{100*x:+.1f}' for x in r['frac'])}\n")
        print('wrote',out)
    else:
        CFG=sys.argv[1] if len(sys.argv)>1 else 'FHC5'; obs=sys.argv[2] if len(sys.argv)>2 else 'pmu'
        evaluate(CFG,obs)
