// combine_ext.C -- generic version of newg4/combine_newg4.C + make_fte.C (FHC), combine_rhc.C +
// make_fte.C (RHC) and combine_comb_fte.C (COMB), for an arbitrary comma-separated observable
// list. The arithmetic is copied expression-for-expression so that, on the released
// observables, the output is bit-identical to the released FTE files (checked by verify_ext.C).
//   root -l -b -q 'gen2d/combine_ext.C("numu.root","numubar.root","out_fte.root","fhc","pmu,costhpi_costhmu")'
//   root -l -b -q 'gen2d/combine_ext.C+' -e 'combine_comb_ext("fhc_fte.root","rhc_fte.root","comb_fte.root","pmu")'
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cmath>

void combine_ext( const char* numu_f, const char* numubar_f, const char* out_f, const char* mode, const char* obs_csv ) {
  const bool rhc = TString( mode ) == "rhc";
  const double PN = rhc ? 2.92348e-10 : 4.43515e-10, PNB = rhc ? 3.52298e-10 : 2.37644e-10, PT = PN + PNB, A = 40.0;
  TFile fn( numu_f ), fnb( numubar_f ), fo( out_f, "recreate" );
  TObjArray* names = TString( obs_csv ).Tokenize( "," );
  for ( int i = 0; i < names->GetEntries(); ++i ) {
    TString n = ( (TObjString*)names->At( i ) )->GetString();
    TH1D* a = (TH1D*)fn.Get( n ); TH1D* b = (TH1D*)fnb.Get( n );
    if ( !a || !b ) { printf( "  missing %s\n", n.Data() ); continue; }
    int nb = a->GetNbinsX();
    TH1D* fte = new TH1D( n + "_fte", "", nb, 0, nb );
    for ( int k = 1; k <= nb; ++k ) {
      double v = A*(a->GetBinContent(k)*PN+b->GetBinContent(k)*PNB)/PT/1e-38;
      double e = A*sqrt(pow(a->GetBinError(k)*PN,2)+pow(b->GetBinError(k)*PNB,2))/PT/1e-38;
      fte->SetBinContent( k, v * a->GetBinWidth(k) );
      fte->SetBinError(   k, e * a->GetBinWidth(k) );
    }
    fte->SetDirectory( &fo ); fte->Write();
    printf( "  %-18s %s sigma_int=%.4f [1e-38 cm2/Ar]\n", n.Data(), mode, fte->Integral() );
  }
  fo.Close();
}

void combine_comb_ext( const char* fhc_fte, const char* rhc_fte, const char* out_fte, const char* obs_csv ) {
  const double PhiF = 4.43515e-10 + 2.37644e-10, PhiR = 2.92348e-10 + 3.52298e-10;
  const double POTF = 8.857e20, POTR = 1.1082e21;
  const double wF = PhiF*POTF, wR = PhiR*POTR, WT = wF+wR;
  TFile fF( fhc_fte ), fR( rhc_fte ), fo( out_fte, "recreate" );
  TObjArray* names = TString( obs_csv ).Tokenize( "," );
  for ( int i = 0; i < names->GetEntries(); ++i ) {
    TString nm = ( (TObjString*)names->At( i ) )->GetString() + "_fte";
    TH1D* a = (TH1D*)fF.Get( nm ); TH1D* b = (TH1D*)fR.Get( nm );
    if ( !a || !b ) { printf( "  missing %s\n", nm.Data() ); continue; }
    TH1D* h = (TH1D*)a->Clone( nm ); h->Reset(); h->SetDirectory( &fo );
    for ( int k = 1; k <= a->GetNbinsX(); ++k ) {
      h->SetBinContent(k,(wF*a->GetBinContent(k)+wR*b->GetBinContent(k))/WT);
      h->SetBinError(k, sqrt(pow(wF*a->GetBinError(k),2)+pow(wR*b->GetBinError(k),2))/WT);
    }
    h->Write();
  }
  fo.Close(); printf( "wrote %s\n", out_fte );
}
