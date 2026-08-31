// export_matrices.C -- dump the additional-smearing matrices A_C as plain text.
//
// A Wiener-SVD result estimates A_C x_true, not x_true, so a prediction can only be
// compared against it after being multiplied by the same A_C. That makes A_C part of
// the measurement: without it the released central values are not usable. This writes
// one TSV per extraction plus an index, so the matrices can be applied without ROOT.
//
//   root -l -b -q 'macros/export_matrices.C("../report/data_release")'
#include <sys/stat.h>

// The central-value double-weighting fix (commit 51af326, 2026-08-30 13:26) moved every
// cross section by 3.8-12.9%, and the Wiener filter -- hence A_C -- depends on the data
// covariance, so a sidecar written before it is not the released measurement. Exporting
// one silently is how a stale number reaches a data release, so the cutoff is enforced
// here rather than trusted to whoever runs this.
static const time_t FAKEDATA_FIX_EPOCH = 1788092802; // git show -s --format=%ct 51af326

static void dump_one(const char* path, const char* tag, const char* outdir, FILE* index){
  struct stat st;
  if ( stat(path, &st) != 0 ) { printf("  MISSING   %s\n", path); return; }
  if ( st.st_mtime < FAKEDATA_FIX_EPOCH ) {
    printf("  STALE     %-30s (written %.19s, predates the fake-data fix) -- SKIPPED\n",
           tag, ctime(&st.st_mtime));
    return;
  }
  TFile* f = TFile::Open(path);
  if (!f || f->IsZombie()) { printf("  MISSING  %s\n", path); return; }
  TH2D* ac = (TH2D*)f->Get("h_A_C");
  if (!ac) { printf("  no h_A_C  %s\n", tag); f->Close(); return; }

  int nx = ac->GetNbinsX(), ny = ac->GetNbinsY();
  TString out = TString::Format("%s/A_C_%s.tsv", outdir, tag);
  FILE* fp = fopen(out.Data(), "w");
  fprintf(fp, "# Additional smearing matrix A_C for %s\n", tag);
  fprintf(fp, "# Compare a prediction p to the published result as (A_C . p), summing over\n");
  fprintf(fp, "# the TRUE index. Row = smeared bin i, column = true bin j, value = A_C[i][j].\n");
  fprintf(fp, "# %d smeared bins x %d true bins.\n", ny, nx);
  fprintf(fp, "smeared_bin");
  for (int j = 1; j <= nx; ++j) fprintf(fp, "\ttrue_%d", j);
  fprintf(fp, "\n");
  double rowsum_min = 1e30, rowsum_max = -1e30;
  for (int i = 1; i <= ny; ++i) {
    fprintf(fp, "%d", i);
    double rs = 0.;
    for (int j = 1; j <= nx; ++j) { double v = ac->GetBinContent(j, i); fprintf(fp, "\t%.10g", v); rs += v; }
    fprintf(fp, "\n");
    rowsum_min = std::min(rowsum_min, rs); rowsum_max = std::max(rowsum_max, rs);
  }
  fclose(fp);

  // A_C is not norm-preserving; the row-sum range is the honest measure of by how much,
  // and it is why summing a Wiener-SVD differential result is not the physical total.
  fprintf(index, "%s\t%d\t%d\t%.6g\t%.6g\t%s\n", tag, ny, nx, rowsum_min, rowsum_max,
          gSystem->BaseName(path));
  printf("  wrote A_C_%-28s  %2dx%-2d  row sums %.3f-%.3f\n", tag, ny, nx, rowsum_min, rowsum_max);
  f->Close();
}

void export_matrices(const char* outdir = "../report/data_release"){
  mkdir(outdir, 0755);
  TString idx = TString::Format("%s/index_A_C.tsv", outdir);
  FILE* index = fopen(idx.Data(), "w");
  fprintf(index, "# Index of released additional-smearing matrices.\n");
  fprintf(index, "# rowsum_min/max show how far A_C departs from norm preservation.\n");
  fprintf(index, "tag\tn_smeared_bins\tn_true_bins\trowsum_min\trowsum_max\tsource_file\n");

  const char* cfgs[3]  = {"FHC5","RHCFULL","COMB"};
  const char* incl[6]  = {"pmu","ppi2bin","costhmu","costhpi","thmupi","thetamu"};
  const char* p1p[10]  = {"pmu","ppi2bin","costhmu","costhpi","thmupi","Wpipr","Whad",
                          "dpt2bin","dphit2bin","dalphat2bin"};

  for (auto c : cfgs) {
    for (auto o : incl) {
      TString p = TString::Format("/data/uboone/processed/closure_hists_xsec_%s_%s.root", c, o);
      if (gSystem->AccessPathName(p)) continue;
      dump_one(p, TString::Format("incl_%s_%s", c, o).Data(), outdir, index);
    }
  }
  for (auto c : cfgs) {
    for (auto o : p1p) {
      TString p = TString::Format("/data/uboone/processed/closure_hists_xsec_ccpi1p_%s_%s.root", c, o);
      if (gSystem->AccessPathName(p)) continue;
      dump_one(p, TString::Format("1p_%s_%s", c, o).Data(), outdir, index);
    }
  }
  fclose(index);
  printf("\n  index: %s\n", idx.Data());
}
