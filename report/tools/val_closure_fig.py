# val_closure_fig.py -- figure of the fake-data closure (unfolded / realised-truth integral) for the
# 18 inclusive extractions, from the TSV written by closure_tables.C.
#     python3 report/tools/val_closure_fig.py closure_summary.tsv OUT.eps
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.startswith('incl\t')]
obs=['pmu','ppi2bin','costhmu','costhpi','thmupi','thetamu']
lab={'pmu':'p_{#mu}','ppi2bin':'p_{#pi}','costhmu':'cos#theta_{#mu}','costhpi':'cos#theta_{#pi}','thmupi':'#theta_{#mu#pi}','thetamu':'#theta_{#mu}'}
cfg=[('FHC5','FHC','#0072B2',20),('RHCFULL','RHC','#D55E00',21),('COMB','combined','#009E73',22)]
R={(r[1],r[2]):float(r[5]) for r in rows}
vals=[R[(c,o)] for c,_,_,_ in cfg for o in obs if (c,o) in R]
mean=sum(vals)/len(vals); print(f"{len(vals)} extractions, mean {mean:.4f}, range {min(vals):.3f}-{max(vals):.3f}")
c=ROOT.TCanvas('c','',1100,600); c.SetLeftMargin(0.13); c.SetBottomMargin(0.13); c.SetRightMargin(0.03); c.SetTopMargin(0.08)
lo=min(0.9,min(vals)-0.02); hi=max(1.1,max(vals)+0.03)
fr=ROOT.TH1D('fr','',len(obs),0,len(obs)); fr.SetMinimum(lo); fr.SetMaximum(hi)
for i,o in enumerate(obs): fr.GetXaxis().SetBinLabel(i+1,lab[o])
fr.GetXaxis().SetLabelSize(0.085); fr.GetYaxis().SetTitle('unfolded / truth'); fr.GetYaxis().SetTitleSize(0.07); fr.GetYaxis().SetLabelSize(0.06); fr.GetYaxis().SetTitleOffset(0.8)
fr.Draw('axis')
one=ROOT.TLine(0,1,len(obs),1); one.SetLineWidth(2); one.Draw()
ml=ROOT.TLine(0,mean,len(obs),mean); ml.SetLineStyle(2); ml.SetLineColor(ROOT.kGray+2); ml.SetLineWidth(2); ml.Draw()
lg=ROOT.TLegend(0.15,0.76,0.8,0.91); lg.SetNColumns(4); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(0.06); keep=[one,ml,lg]
for j,(c_,name,col,mk) in enumerate(cfg):
    g=ROOT.TGraph(); k=0
    for i,o in enumerate(obs):
        if (c_,o) in R: g.SetPoint(k,i+0.3+0.2*j,R[(c_,o)]); k+=1
    g.SetMarkerStyle(mk); g.SetMarkerSize(2.4); g.SetMarkerColor(ROOT.TColor.GetColor(col)); g.Draw('P'); lg.AddEntry(g,name,'p'); keep.append(g)
lg.AddEntry(ml,f'mean {mean:.3f}','l'); lg.Draw()
tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.06); tx.DrawLatex(0.13,0.935,'MicroBooNE NuMI fake-data closure')
c.SaveAs(sys.argv[2])
