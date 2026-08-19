// ---------------------------------------------------------------------------
// The NuMI beam geometry as MicroBooNE sees it, drawn to settle one distinction:
// the 7.7 degree OFF-AXIS angle (where the detector sits relative to the NuMI
// beamline) and the 28.6 degree angle between the arriving neutrinos and the
// DETECTOR z axis are angles between different pairs of directions. Only the
// second enters the analysis.
//
// Panels (a) and (b) are drawn to scale - the user-units-per-pixel are equalised
// in x and y - so every marked angle is the true angle. Panel (b) is the x-z
// projection; the quoted 28.6 degrees is the full 3D angle.
// ---------------------------------------------------------------------------

// Equalise x and y scales so arcs stay circular and drawn angles are honest.
static double sq_range( TPad* p, double x1, double x2, double yc ) {
  p->Update();
  double w = p->GetWw()*p->GetAbsWNDC(), h = p->GetWh()*p->GetAbsHNDC();
  double yr = (x2-x1)*h/w;
  p->Range( x1, yc-0.5*yr, x2, yc+0.5*yr );
  return 0.5*yr;
}
static void dr_arrow( double x1,double y1,double x2,double y2,int col,double w,double sz=0.020 ) {
  TArrow* a=new TArrow(x1,y1,x2,y2,sz,"|>");
  a->SetLineColor(col); a->SetFillColor(col); a->SetLineWidth(w); a->Draw();
}
static void dr_tex( double x,double y,const char* s,int col,double sz,int align=11,double ang=0 ) {
  TLatex* t=new TLatex(x,y,s);
  t->SetTextColor(col); t->SetTextSize(sz); t->SetTextAlign(align);
  t->SetTextAngle(ang); t->SetTextFont(42); t->Draw();
}
static void dr_ndc( double x,double y,const char* s,int col,double sz,int align=11 ) {
  TLatex* t=new TLatex(x,y,s); t->SetNDC();
  t->SetTextColor(col); t->SetTextSize(sz); t->SetTextAlign(align); t->SetTextFont(42); t->Draw();
}
static void dr_arc( double x,double y,double r,double a1,double a2,int col,int style=1 ) {
  TArc* c=new TArc(x,y,r,a1,a2);
  c->SetFillStyle(0); c->SetLineColor(col); c->SetLineStyle(style);
  c->SetLineWidth(2); c->SetNoEdges(); c->Draw("only");
}
static void dr_line( double x1,double y1,double x2,double y2,int col,int style,int w ) {
  TLine* l=new TLine(x1,y1,x2,y2);
  l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(w); l->Draw();
}

void make_numi_geometry_figure() {
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineStyle(1);

  const int kBeam=kAzure+2, kNu=kRed+1, kDet=kGray+3, kMu=kGreen+3;

  TCanvas* c=new TCanvas("c","",1740,640);
  TPad* pa=new TPad("pa","",0.000,0.0,0.350,1.0);
  TPad* pb=new TPad("pb","",0.350,0.0,0.680,1.0);
  TPad* pc=new TPad("pc","",0.680,0.0,1.000,1.0);
  pa->Draw(); pb->Draw(); pc->Draw();

  // ===================== (a) beamline geometry, to scale =====================
  pa->cd();
  sq_range( pa, -40., 850., 40. );
  const double LZ=672.7, LT=91.1;                       // metres

  dr_arrow( 0,0, 800,0, kBeam, 3 );
  dr_tex( 800, -42, "NuMI beamline axis", kBeam, 0.036, 31 );

  dr_arrow( 0,0, LZ,LT, kNu, 3 );                       // the neutrinos' path
  dr_arc( 0,0, 300., 0., 7.71, kNu );
  dr_tex( 316, 24, "7.7#circ", kNu, 0.044, 12 );
  dr_tex( 316, -6, "off-axis angle", kGray+3, 0.030, 12 );

  dr_line( LZ,0, LZ,LT, kGray+2, 2, 1 );
  dr_tex( LZ+16, 0.5*LT, "91 m", kGray+3, 0.032, 12 );

  TMarker* tg=new TMarker(0,0,20); tg->SetMarkerSize(1.5); tg->Draw();
  dr_tex( 16, -42, "NuMI target", kBlack, 0.036, 11 );

  dr_tex( 300, 84, "678.8 m", kNu, 0.032, 21, 7 );

  TBox* db=new TBox(LZ-40,LT-26,LZ+40,LT+26);
  db->SetFillColor(kOrange-2); db->Draw();
  TBox* dbo=new TBox(LZ-40,LT-26,LZ+40,LT+26);
  dbo->SetFillStyle(0); dbo->SetLineColor(kBlack); dbo->SetLineWidth(2); dbo->Draw("l");
  dr_tex( LZ, LT+42, "MicroBooNE", kBlack, 0.036, 21 );

  dr_ndc( 0.045, 0.93, "(a)  beamline geometry", kBlack, 0.050 );
  dr_ndc( 0.045, 0.875, "where the detector sits relative", kGray+3, 0.033 );
  dr_ndc( 0.045, 0.835, "to the NuMI beam", kGray+3, 0.033 );
  dr_ndc( 0.045, 0.115, "This angle says nothing about", kGray+3, 0.031 );
  dr_ndc( 0.045, 0.075, "how the TPC is oriented.", kGray+3, 0.031 );

  // ================= (b) the same neutrinos in the TPC, to scale =============
  pb->cd();
  sq_range( pb, -400., 1150., 128. );
  const double VX=430., VY=128.;
  const double THN=28.62*TMath::DegToRad();             // nu, above z
  const double THM=-12.0*TMath::DegToRad();             // a muon, below z

  TBox* tpc=new TBox(0,0,1036,256); tpc->SetFillColor(kOrange-9); tpc->Draw();
  TBox* tpo=new TBox(0,0,1036,256);
  tpo->SetFillStyle(0); tpo->SetLineColor(kBlack); tpo->SetLineWidth(2); tpo->Draw("l");
  dr_tex( 20, 222, "TPC", kBlack, 0.032, 11 );

  // detector z axis through the vertex
  dr_line( VX-380, VY, VX+600, VY, kDet, 2, 3 );
  dr_arrow( VX+500, VY, VX+660, VY, kDet, 3 );
  dr_tex( VX+430, VY+14, "detector z", kDet, 0.037, 31 );

  // the neutrino: arriving, and its onward ray
  double L1=620.;
  dr_arrow( VX-L1*TMath::Cos(THN), VY-L1*TMath::Sin(THN), VX, VY, kNu, 4 );
  dr_line( VX, VY, VX+470*TMath::Cos(THN), VY+470*TMath::Sin(THN), kNu, 2, 2 );
  dr_tex( VX-L1*TMath::Cos(THN)-20, VY-L1*TMath::Sin(THN)-10, "#nu_{#mu}", kNu, 0.050, 31 );

  // the nu-to-z angle: the number that matters
  dr_arc( VX,VY, 260., 0., 28.62, kNu );
  dr_tex( VX+300, VY+74, "28.6#circ", kNu, 0.046, 12 );

  // a muon, and the two angles it would be assigned
  double L2=560.;
  dr_arrow( VX,VY, VX+L2*TMath::Cos(THM), VY+L2*TMath::Sin(THM), kMu, 3 );
  dr_tex( VX+L2*TMath::Cos(THM)+18, VY+L2*TMath::Sin(THM)-8, "#mu", kMu, 0.050, 11 );
  dr_arc( VX,VY, 225., -12., 0., kDet, 1 );             // theta from z
  dr_arc( VX,VY, 415., -12., 28.62, kNu, 2 );           // theta from the neutrino

  dr_ndc( 0.055, 0.93, "(b)  inside the detector", kBlack, 0.050 );
  dr_ndc( 0.055, 0.875, "the same #nu, in TPC coordinates", kGray+3, 0.033 );
  dr_ndc( 0.055, 0.20, "the same #mu is assigned", kBlack, 0.033 );
  dr_ndc( 0.055, 0.150, "#theta_{#mu} = 12#circ  from detector z", kDet, 0.036 );
  dr_ndc( 0.055, 0.100, "#theta_{#mu} = 40.6#circ  from the #nu", kNu, 0.036 );
  dr_ndc( 0.055, 0.045, "the generators mean the second", kGray+3, 0.031 );

  // ===================== (c) what it does to a distribution ==================
  pc->cd();
  gPad->SetLeftMargin(0.185); gPad->SetBottomMargin(0.155);
  gPad->SetTopMargin(0.215);  gPad->SetRightMargin(0.045);
  TFile* fh=TFile::Open("/home/t2k/nowak/.claude/jobs/a2cc9f4f/tmp/costh_hists.root");
  TH1D* hz=(TH1D*)fh->Get("cmu_z");
  TH1D* hb=(TH1D*)fh->Get("cmu_b");
  hz->Scale(1e-3); hb->Scale(1e-3);
  hz->SetLineColor(kDet); hz->SetLineWidth(3); hz->SetLineStyle(2);
  hb->SetLineColor(kNu);  hb->SetLineWidth(3);
  hz->SetTitle("");
  hz->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hz->GetYaxis()->SetTitle("true signal events (10^{3})");
  hz->GetXaxis()->SetTitleSize(0.055); hz->GetYaxis()->SetTitleSize(0.052);
  hz->GetXaxis()->SetLabelSize(0.046); hz->GetYaxis()->SetLabelSize(0.046);
  hz->GetXaxis()->SetTitleOffset(1.20); hz->GetYaxis()->SetTitleOffset(1.55);
  hz->SetMaximum( 1.20*TMath::Max(hz->GetMaximum(),hb->GetMaximum()) );
  hz->Draw("hist"); hb->Draw("hist same");

  double edges[4]={0.45,0.65,0.80,0.90};
  for ( int i=0;i<4;++i )
    dr_line( edges[i],0, edges[i],hz->GetMaximum()*0.72, kGray+1, 3, 1 );

  TLegend* lg=new TLegend(0.22,0.60,0.66,0.755);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.043);
  lg->AddEntry(hz,"about detector z","l");
  lg->AddEntry(hb,"about the #nu direction","l");
  lg->Draw();

  dr_ndc( 0.045, 0.93, "(c)  the consequence", kBlack, 0.050 );
  dr_ndc( 0.045, 0.875, "true CC1#mu1#pi signal, FHC", kGray+3, 0.033 );
  dr_ndc( 0.045, 0.825, "64.6% of events change bin", kNu, 0.036 );
  dr_ndc( 0.60, 0.50, "dotted: analysis", kGray+2, 0.030 );
  dr_ndc( 0.60, 0.455, "bin edges", kGray+2, 0.030 );

  c->SaveAs("figures/numi_geometry.pdf");
  c->SaveAs("/home/t2k/nowak/.claude/jobs/a2cc9f4f/tmp/numi_geometry.png");
  printf("  wrote figures/numi_geometry.pdf\n");
}
