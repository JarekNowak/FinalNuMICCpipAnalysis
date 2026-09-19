"""mcs_eval.py -- the MCS momentum-scale term for d sigma/dp_mu (FHC), from the DATA-SIDE variation.

The fake data are re-thrown with every MCS-measured muon momentum (reco != range, 63% of the selected
sample) scaled by 1 +- 0.05 while the simulation, the response and the selection flag stay nominal
(macros/mcs_scale_fakedata.C; extractions in rebuild_alt/ccpi_FHC5_pmu_mcsdata{up,dn}05). The earlier
bin-configuration scaling moved data and simulation together and therefore bounded only the
regulariser's response to a remapped variable, not the systematic; it is superseded by this.

Per bin: half-difference D_i = (x_up - x_dn)/2 of the unfolded result, quoted as a fraction of the
nominal and in units of the released total uncertainty; the covariance C_ij = D_i D_j (fully
correlated unisim) is written in the release units (bin integrals, like cov_*.txt) as
data_release/cov/incl_FHC5_pmu/cov_MCSscale.txt, with cov_total_plusMCS.txt = cov_total + C.

    python3 report/tools/mcs_eval.py
"""
import numpy as np, uproot, os
RB='/data/uboone/processed/rebuild_alt/'; LIVE='/data/uboone/processed/'
REL='/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/report/data_release/cov/incl_FHC5_pmu/'
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
nom=uproot.open(LIVE+'closure_hists_xsec_FHC5_pmu.root')['h_unfolded_nuwro']
v=nom.values(); e=nom.errors(); w=np.diff(nom.axis().edges()); n=len(v)
up=uproot.open(RB+'closure_hists_xsec_ccpi_FHC5_pmu_mcsdataup05.root')['h_unfolded_nuwro'].values()
dn=uproot.open(RB+'closure_hists_xsec_ccpi_FHC5_pmu_mcsdatadn05.root')['h_unfolded_nuwro'].values()
D=(up-dn)/2
print('bin   nominal   up/nom   dn/nom   half-diff/nom   half-diff/sigma_tot')
for i in range(n): print(f'{i+1:>3}  {v[i]:8.4f}  {up[i]/v[i]:7.3f}  {dn[i]/v[i]:7.3f}  {D[i]/v[i]:+13.3f}  {D[i]/e[i]:+12.2f}')
I=(v*w).sum(); print(f'integral: nominal {I:.4f}  up {(up*w).sum()/I:.3f}  dn {(dn*w).sum()/I:.3f}  half-diff {((up-dn)*w).sum()/2/I:+.3%}')
print(f'mean |half-diff|/nom {np.mean(abs(D)/v):.3%}, max |half-diff|/sigma {np.max(abs(D)/e):.2f}')
Db=D*w                      # release units: bin integrals
C=np.outer(Db,Db); Ct=cov(REL+'cov_total.txt')
assert Ct.shape==(n,n)
write(REL+'cov_MCSscale.txt',C,n); write(REL+'cov_total_plusMCS.txt',Ct+C,n)
print('total uncertainty per bin, before -> after adding the MCS term (% of bin):')
print(' '.join(f'{100*np.sqrt(Ct[i,i])/(v[i]*w[i]):.1f}->{100*np.sqrt(Ct[i,i]+C[i,i])/(v[i]*w[i]):.1f}' for i in range(n)))
print('wrote', REL+'cov_MCSscale.txt', 'and cov_total_plusMCS.txt')
