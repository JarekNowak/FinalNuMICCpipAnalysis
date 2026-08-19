// Integrated (single-bin) W_pipr generator predictions, to accompany
// configs/ccpi1p_Wpipr_bin_config_1bin.txt. W_pipr does not support a shape
// measurement at this exposure (see that config's README), so it is quoted as an
// integrated cross section and the predictions must be integrated to match.
//
// Follows the same convention as rebin_fte.C: the FTE bin content is the per-bin
// cross section (dsigma/dx * width), so collapsing to one bin means SUMMING the
// contents, with errors added in quadrature - not integrating a density. The output
// histogram spans [0,1] with a single unit-width bin, exactly as remap() produces.
void make_wpipr_1bin(){
  const char* gens[4]={"genie","gibuu","neut","nuwro"};
  for ( auto g : gens ) {
    TString in = Form("%s_wtki_fte.root", g);
    TFile* f = TFile::Open(in);
    if ( !f || f->IsZombie() ) { printf("  %-6s : cannot open %s\n", g, in.Data()); continue; }
    TH1D* h = (TH1D*)f->Get("Wpipr_fte");
    if ( !h ) { printf("  %-6s : no Wpipr_fte\n", g); f->Close(); continue; }
    double c = 0., e2 = 0.;
    for ( int i = 1; i <= h->GetNbinsX(); ++i ) {
      c  += h->GetBinContent(i);
      e2 += h->GetBinError(i) * h->GetBinError(i);
    }
    TString out = Form("%s_wtki_1bin_fte.root", g);
    TFile o( out, "RECREATE" );
    TH1D h1( "Wpipr_fte", "", 1, 0., 1. );
    h1.SetBinContent( 1, c );
    h1.SetBinError  ( 1, TMath::Sqrt(e2) );
    h1.Write();
    o.Close();
    printf("  %-6s : %d bins, sum %.4f -> 1 bin, content %.4f  (%s)\n",
      g, h->GetNbinsX(), c, h1.GetBinContent(1), out.Data());
    f->Close();
  }
}
