// dsigma_style.h -- presentation shared by the differential cross-section result figures
// (dsigma_current.C: inclusive; dsigma_ccpi1p.C: proton-tagged). Written 2026-09-27 for the
// review points on the result figures:
//   * panel titles and axis labels in physics notation, with units, looked up from the file key
//     (d sigma / d x per unit of x on the vertical axis);
//   * text sized in pixels (font precision 3) so that it keeps its printed size whatever the pad
//     layout, and one legend per figure (a free grid cell or a strip above the panels) instead of a
//     crowded legend in every panel;
//   * the open-ended top bin drawn hatched and labelled with its physical range: its plotted upper
//     edge is a conventional boundary and the drawn height is the bin integral / plotted width;
//   * review terminology (CV pseudo-data), and a panel without generator curves says so, since the
//     figure legend lists the generators whenever any of its panels draws them.
// Presentation only. The plotted numbers are read from the closure sidecars exactly as before, and
// every drawn curve is printed as a "[plotted]" line so that a regeneration can be diffed.
#pragma once
#include "TH1D.h"
#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "THLimitsFinder.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

struct DsObs { const char* key; const char* sym; const char* unit; bool open_top; };

// Notation of every observable key of the result files. open_top marks an upper truth bin without an
// upper edge (the starred edges of report/tables/binning.tex and binning_1p.tex); the edge values
// themselves are always read from the histograms.
inline const DsObs* ds_obs(const std::string& key) {
  static const DsObs T[] = {
    {"pmu",         "p_{#mu}",          "GeV/c",     true },
    {"ppi",         "p_{#pi}",          "GeV/c",     true },
    {"ppi2bin",     "p_{#pi}",          "GeV/c",     true },
    {"ppi3bin",     "p_{#pi}",          "GeV/c",     true },
    {"costhmu",     "cos#theta_{#mu}",  "",          false},
    {"costhpi",     "cos#theta_{#pi}",  "",          false},
    {"thmupi",      "#theta_{#mu#pi}",  "rad",       false},
    {"Wpipr",       "W_{#pi p}",        "GeV/c^{2}", true },
    {"Whad",        "W_{had}",          "GeV/c^{2}", true },
    {"dpt2bin",     "#deltap_{T}",      "GeV/c",     true },
    {"dalphat2bin", "#delta#alpha_{T}", "deg",       false},
    {"dphit3bin",   "#delta#phi_{T}",   "deg",       false},
    {"pn2bin",      "p_{n}",            "GeV/c",     true },
    {"pp",          "p_{p}",            "GeV/c",     true },
  };
  for ( const auto& o : T ) if ( key == o.key ) return &o;
  static const DsObs unknown = {"", "x", "", false};
  printf("  [warn] dsigma_style: no notation for key '%s'\n", key.c_str());
  return &unknown;
}

inline std::string ds_xtitle(const DsObs* o) {
  std::string s = o->sym;
  if ( *o->unit ) s += std::string(" [") + o->unit + "]";
  return s;
}

// d sigma / d x [10^-38 cm^2 / (unit of x) / Ar]
inline std::string ds_ytitle(const DsObs* o) {
  std::string u = o->unit;
  if ( u.find('/') != std::string::npos ) u = "(" + u + ")";
  return std::string("d#sigma/d") + o->sym + " [10^{-38} cm^{2}" + ( u.empty() ? "" : "/" + u ) + "/Ar]";
}

inline std::string ds_open_label(const DsObs* o, double lo) {
  char b[128]; snprintf(b, sizeof b, "%s > %g%s%s", o->sym, lo, *o->unit ? " " : "", o->unit);
  return b;
}

// Text sizes in pixels (font 43), pad margins as fractions of the pad.
struct DsStyle {
  int title_px, label_px, binlabel_px, legend_px, ptitle_px, note_px;
  double lm, rm, tm, bm, xoff, yoff;
};

// One panel: the curves as drawn (already styled), the observable key and the open-bin edge.
struct DsPanel {
  std::string key;
  TH1D* data = nullptr; TH1D* truth = nullptr; TH1D* tune = nullptr;
  std::vector<TH1D*> gen; std::vector<int> gidx;
  double open_lo = -1.;      // physical lower edge of the open top bin (< 0: no open bin)
  double open_hi = -1.;      // its conventional (plotting) upper edge; stated in the panel when the
                             // axis does not show it (index_axis)
  bool index_axis = false;   // bins drawn at equal width, labelled by their ranges
  std::string gen_note;      // drawn in a panel without generator curves (set by the caller)
};

static const char* DS_GLAB[4] = {"GENIE", "GiBUU", "NEUT", "NuWro"};

inline void ds_dump(const char* cfg, const char* fig, const char* obs, const char* curve, TH1* h) {
  if ( !h ) return;
  printf("[plotted] %s %s %s %s n=%d", cfg, fig, obs, curve, h->GetNbinsX());
  for ( int b = 1; b <= h->GetNbinsX(); ++b )
    printf(" | %.9g %.9g %.12g %.12g", h->GetXaxis()->GetBinLowEdge(b), h->GetXaxis()->GetBinUpEdge(b),
           h->GetBinContent(b), h->GetBinError(b));
  printf("\n"); fflush(stdout);
}

// The hatched fill that marks an open-ended bin, in the panels and in the legend.
inline void ds_open_fill(TAttFill* a) { a->SetFillStyle(3554); a->SetFillColor(kGray + 1); }

// Draw one panel into the current pad. The vertical range (set by the caller) is unchanged.
inline void ds_draw_panel(DsPanel& P, const DsStyle& S, const char* cfg, const char* fig,
                          std::vector<TObject*>& keep) {
  const DsObs* O = ds_obs(P.key);
  TH1D* h = P.data;
  gPad->SetLeftMargin(S.lm); gPad->SetRightMargin(S.rm);
  gPad->SetTopMargin(S.tm);  gPad->SetBottomMargin(S.bm);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(ds_xtitle(O).c_str());
  h->GetYaxis()->SetTitle(ds_ytitle(O).c_str());
  for ( TAxis* a : { h->GetXaxis(), h->GetYaxis() } ) {
    a->SetTitleFont(43); a->SetLabelFont(43);
    a->SetTitleSize(S.title_px); a->SetLabelSize(S.label_px);
    a->SetNdivisions(505);
  }
  if ( P.index_axis ) h->GetXaxis()->SetLabelSize(S.binlabel_px);
  h->GetXaxis()->SetTitleOffset(S.xoff); h->GetYaxis()->SetTitleOffset(S.yoff);
  h->GetYaxis()->SetMaxDigits(3);
  h->Draw("E1");
  const int n = h->GetNbinsX();
  const double fx1 = S.lm, fx2 = 1. - S.rm, fy1 = S.bm, fy2 = 1. - S.tm;
  const double padh = gPad->GetWh() * gPad->GetAbsHNDC();          // pad height in pixels
  const double note_y = fy2 - 0.02;                                 // top of the note stack (NDC)
  // Notes at the top right of the frame, above every curve (the 1.35 headroom of the vertical range
  // keeps the curves below 74% of the frame height): the physical range of an open-ended top bin;
  // on an equal-width axis, where the conventional width of that bin is not visible, the width;
  // and, in a panel without generator curves, why they are absent.
  const int small_px = int( 0.85 * S.note_px );
  std::vector<std::pair<std::string, int>> notes;
  if ( P.open_lo >= 0. ) notes.push_back({ ds_open_label(O, P.open_lo), S.note_px });
  if ( P.open_lo >= 0. && P.index_axis && P.open_hi > P.open_lo )
    notes.push_back({ Form("height = integral / %g %s", P.open_hi - P.open_lo, O->unit), small_px });
  if ( P.gen.empty() && !P.gen_note.empty() ) notes.push_back({ P.gen_note, small_px });
  double stack_px = 0.; for ( auto& nt : notes ) stack_px += 1.3 * nt.second;
  if ( P.open_lo >= 0. ) {
    // hatched column over the open bin, stopping just below the notes
    const double top_ndc = note_y - stack_px / padh - 0.01;
    const double ytop = h->GetMinimum() + ( top_ndc - fy1 ) / ( fy2 - fy1 ) * ( h->GetMaximum() - h->GetMinimum() );
    TBox* box = new TBox(h->GetXaxis()->GetBinLowEdge(n), h->GetMinimum(), h->GetXaxis()->GetXmax(), ytop);
    ds_open_fill(box); box->SetLineWidth(0); box->Draw(); keep.push_back(box);
  }
  if ( P.truth ) P.truth->Draw("hist same");
  if ( P.tune )  P.tune->Draw("hist same");
  for ( auto g : P.gen ) g->Draw("hist same");
  h->Draw("E1 same");
  h->Draw("sameaxis");   // ticks and frame back on top of the hatched fill
  TLatex t; t.SetNDC(); t.SetTextFont(43);
  // panel title: the observable
  t.SetTextSize(S.ptitle_px); t.SetTextAlign(21);
  t.DrawLatex(0.5 * (fx1 + fx2), fy2 + 0.022, O->sym);
  // the note stack, right-aligned at the top of the frame
  t.SetTextAlign(33); t.SetTextColor(kGray + 3);
  double ny = note_y;
  for ( auto& nt : notes ) { t.SetTextSize(nt.second); t.DrawLatex(fx2 - 0.02, ny, nt.first.c_str()); ny -= 1.3 * nt.second / padh; }
  t.SetTextColor(kBlack);
  // ROOT labels only "round" ticks, so an axis that does not start at zero would read as starting
  // at the first labelled tick: write the lower limit under the left edge of the frame.
  // It goes in the row of the tick labels, or one row lower if the first labelled tick is too close.
  const double xmin = h->GetXaxis()->GetXmin(), xmax = h->GetXaxis()->GetXmax();
  if ( !P.index_axis && xmin > 0. ) {
    double bl, bh, bw; int nb;
    THLimitsFinder::Optimize(xmin, xmax, 5, bl, bh, nb, bw, "");
    double tick = bl; while ( tick < xmin - 1e-9 * ( xmax - xmin ) ) tick += bw;
    const double padw = gPad->GetWw() * gPad->GetAbsWNDC();
    const double gap_px = ( tick - xmin ) / ( xmax - xmin ) * ( fx2 - fx1 ) * padw;
    const double need_px = 0.55 * S.label_px * ( strlen(Form("%g", xmin)) + strlen(Form("%g", tick)) ) / 2. + 8.;
    const double y = S.bm - 0.012 - ( gap_px < need_px ? 1.15 * S.label_px / padh : 0. );
    t.SetTextSize(S.label_px); t.SetTextAlign(23);
    t.DrawLatex(fx1, y, Form("%g", xmin));
  }
  ds_dump(cfg, fig, P.key.c_str(), "data", h); ds_dump(cfg, fig, P.key.c_str(), "truth", P.truth);
  ds_dump(cfg, fig, P.key.c_str(), "tune", P.tune);
  for ( size_t k = 0; k < P.gen.size(); ++k ) ds_dump(cfg, fig, P.key.c_str(), DS_GLAB[P.gidx[k]], P.gen[k]);
}

// The legend of a figure: every curve drawn in any of its panels, and the open-bin fill.
// strip = true: entries in rows across a strip above the panels; false: one column in a free cell.
inline void ds_draw_legend(const std::vector<DsPanel*>& panels, const DsStyle& S, bool strip,
                           std::vector<TObject*>& keep) {
  TH1* data = nullptr; TH1* truth = nullptr; TH1* tune = nullptr; TH1* gen[4] = {};
  bool open = false;
  for ( auto p : panels ) {
    if ( !data ) data = p->data;
    if ( !truth ) truth = p->truth;
    if ( !tune ) tune = p->tune;
    for ( size_t k = 0; k < p->gen.size(); ++k ) if ( !gen[p->gidx[k]] ) gen[p->gidx[k]] = p->gen[k];
    if ( p->open_lo >= 0. ) open = true;
  }
  TBox* fill = new TBox(0, 0, 1, 1); ds_open_fill(fill); fill->SetLineWidth(0); keep.push_back(fill);
  auto fmt = [&](TLegend* l) {
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextFont(43); l->SetTextSize(S.legend_px);
    keep.push_back(l);
  };
  if ( strip ) {
    // three rows: the measurement and the two references; the generators (when any panel draws
    // them); the open-bin fill
    TLegend* a = new TLegend(0.02, 0.68, 0.98, 0.99); fmt(a); a->SetNColumns(3); a->SetMargin(0.14);
    if ( data )  a->AddEntry(data,  "Unfolded CV pseudo-data", "lep");
    if ( truth ) a->AddEntry(truth, "Pseudo-data truth (A_{C}-smeared)", "l");
    if ( tune )  a->AddEntry(tune,  "MicroBooNE tune", "l");
    a->Draw();
    int ng = 0; for ( int g = 0; g < 4; ++g ) if ( gen[g] ) ++ng;
    if ( ng ) {
      TLegend* g2 = new TLegend(0.02, 0.36, 0.98, 0.67); fmt(g2); g2->SetNColumns(4); g2->SetMargin(0.2);
      for ( int g = 0; g < 4; ++g ) if ( gen[g] ) g2->AddEntry(gen[g], DS_GLAB[g], "l");
      g2->Draw();
    }
    if ( open ) {
      TLegend* b = new TLegend(0.02, 0.03, 0.98, 0.34); fmt(b); b->SetMargin(0.055);
      b->AddEntry(fill, "Open-ended top bin, drawn over a conventional width: height = bin integral / plotted width", "f");
      b->Draw();
    }
  } else {
    int ne = ( data ? 1 : 0 ) + ( truth ? 2 : 0 ) + ( tune ? 1 : 0 ) + ( open ? 2 : 0 );
    for ( int g = 0; g < 4; ++g ) if ( gen[g] ) ++ne;
    TLegend* a = new TLegend(0.06, std::max(0.04, 0.96 - 0.1 * ne), 0.98, 0.96); fmt(a); a->SetMargin(0.20);
    if ( data )  a->AddEntry(data,  "Unfolded CV pseudo-data", "lep");
    if ( truth ) a->AddEntry(truth, "#splitline{Pseudo-data truth}{(A_{C}-smeared)}", "l");
    if ( tune )  a->AddEntry(tune,  "MicroBooNE tune", "l");
    for ( int g = 0; g < 4; ++g ) if ( gen[g] ) a->AddEntry(gen[g], DS_GLAB[g], "l");
    if ( open ) a->AddEntry(fill, "#splitline{Open-ended top bin:}{height = integral / plotted width}", "f");
    a->Draw();
  }
}

// One figure: ncol x nrow panels. With a free grid cell the legend goes there; otherwise into a
// strip of strip_px pixels across the top. Saves to `out`.
inline void ds_figure(const char* out, std::vector<DsPanel*> panels, int ncol, int nrow, int W, int Hpanel,
                      int strip_px, const DsStyle& S, const char* cfg, std::vector<TObject*>& keep) {
  const int np = (int)panels.size();
  const bool strip = ( np >= ncol * nrow );
  const int H = nrow * Hpanel + ( strip ? strip_px : 0 );
  TCanvas c(Form("ds_%s_%d", cfg, (int)keep.size()), "", W, H);
  const double ytop = strip ? 1. - double(strip_px) / H : 1.;
  TString fig = gSystem->BaseName(out); fig.ReplaceAll(".pdf", "");
  std::vector<TPad*> pads;
  for ( int r = 0; r < nrow; ++r ) for ( int k = 0; k < ncol; ++k ) {
    const double y2 = ytop - r * ytop / nrow, y1 = ytop - ( r + 1 ) * ytop / nrow;
    TPad* p = new TPad(Form("p_%s_%d_%d", fig.Data(), r, k), "", double(k) / ncol, y1, double(k + 1) / ncol, y2);
    p->SetFillStyle(4000); c.cd(); p->Draw(); pads.push_back(p);
  }
  for ( int i = 0; i < np; ++i ) { pads[i]->cd(); ds_draw_panel(*panels[i], S, cfg, fig.Data(), keep); }
  if ( strip ) {
    c.cd(); TPad* lp = new TPad(Form("leg_%s", fig.Data()), "", 0., ytop, 1., 1.);
    lp->SetFillStyle(4000); lp->Draw(); lp->cd(); ds_draw_legend(panels, S, true, keep);
  } else {
    pads[np]->cd(); ds_draw_legend(panels, S, false, keep);
  }
  c.SaveAs(out);
  printf("wrote %s\n", out);
}
