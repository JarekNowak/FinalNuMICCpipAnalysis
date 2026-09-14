# val_pulls_fig.py -- figure of the fake-data ensemble pulls (statistical covariance) for the three
# ensemble extractions, against a unit Gaussian. Same pull construction as ensemble_stat_pulls.py.
#     python3 report/tools/val_pulls_fig.py OUT.eps
import sys, math, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
P='/data/uboone/processed/ens/'
def cov(path):
    C={}
    for l in open(path):
        p=l.split()
        if p[0] not in ('numXbins','numYbins','xbin') and len(p)==3: C[(int(p[0]),int(p[1]))]=float(p[2])
    return C
def pulls(obs):
    out=[]
    for t in range(1,33):
        f=ROOT.TFile.Open(P+f'closure_hists_xsec_{obs}_t{t}.root')
        if not f or f.IsZombie(): continue
        u=f.Get('h_unfolded_nuwro'); r=f.Get('h_genie_tune')
        Cs=cov(P+f'statcov_{obs}_t{t}/unfold_output/mat_table_cov_DataStats.txt'); Ct=cov(P+f'statcov_{obs}_t{t}/unfold_output/mat_table_cov_total.txt')
        for b in range(1,u.GetNbinsX()+1):
            e=u.GetBinError(b); k=e/math.sqrt(Ct[(b-1,b-1)]); es=math.sqrt(Cs[(b-1,b-1)])*k
            out.append((u.GetBinContent(b)-r.GetBinContent(b))/es)
        f.Close()
    return out
specs=[('pmu','d#sigma/dp_{#mu}'),('costhmu','d#sigma/dcos#theta_{#mu}'),('total','one-bin total')]
c=ROOT.TCanvas('c','',1500,520); c.Divide(3,1,0.004,0.004); keep=[]
for i,(obs,ttl) in enumerate(specs):
    p=pulls(obs); n=len(p); m=sum(p)/n; w=math.sqrt(sum((x-m)**2 for x in p)/(n-1))
    c68=sum(1 for x in p if abs(x)<1)/n; c95=sum(1 for x in p if abs(x)<1.96)/n
    print(f"{obs}: n={n} mean {m:+.2f} width {w:.2f} cov68 {100*c68:.1f}% cov95 {100*c95:.1f}%")
    c.cd(i+1); ROOT.gPad.SetLeftMargin(0.15); ROOT.gPad.SetBottomMargin(0.15); ROOT.gPad.SetTopMargin(0.1)
    nb=12 if n>100 else 8
    h=ROOT.TH1D(f'h_{obs}','',nb,-4,4); [h.Fill(max(-3.999,min(3.999,x))) for x in p]; keep.append(h)
    h.SetLineColor(ROOT.TColor.GetColor('#0072B2')); h.SetFillColor(ROOT.TColor.GetColor('#9ECAE9')); h.SetLineWidth(2)
    g=ROOT.TF1(f'g_{obs}','[0]*TMath::Gaus(x,0,1,1)',-4,4); g.SetParameter(0,n*h.GetBinWidth(1)); g.SetLineColor(ROOT.kBlack); g.SetLineStyle(2); g.SetLineWidth(2); keep.append(g)
    h.SetMaximum(1.45*max(h.GetMaximum(),g.GetMaximum(-4,4)))
    h.GetXaxis().SetTitle('pull (statistical covariance)'); h.GetYaxis().SetTitle('bins #times members')
    for ax in (h.GetXaxis(),h.GetYaxis()): ax.SetTitleSize(0.055); ax.SetLabelSize(0.048)
    h.Draw('hist'); g.Draw('same')
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.058); tx.DrawLatex(0.18,0.92,ttl)
    tx.SetTextSize(0.048); tx.DrawLatex(0.18,0.82,f'width {w:.2f},  mean {m:+.2f}'); tx.DrawLatex(0.18,0.75,f'68%: {100*c68:.0f}%   95%: {100*c95:.0f}%')
    keep.append(tx)
c.SaveAs(sys.argv[1])
