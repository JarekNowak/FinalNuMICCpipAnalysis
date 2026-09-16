// verify_ext.C -- max relative difference, bin by bin, between histogram NAME in two files.
//   root -l -b -q 'gen2d/verify_ext.C("a.root","b.root","pmu,ppi")'   -> [VERIFY] lines
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cmath>
void verify_ext( const char* fa, const char* fb, const char* names_csv ) {
  TFile A( fa ), B( fb );
  TObjArray* names = TString( names_csv ).Tokenize( "," );
  for ( int i = 0; i < names->GetEntries(); ++i ) {
    TString n = ( (TObjString*)names->At( i ) )->GetString();
    TH1* a = (TH1*)A.Get( n ), *b = (TH1*)B.Get( n );
    if ( !a || !b ) { printf( "[VERIFY] %-16s MISSING (%s %s)\n", n.Data(), a ? "" : fa, b ? "" : fb ); continue; }
    if ( a->GetNbinsX() != b->GetNbinsX() ) { printf( "[VERIFY] %-16s NBINS %d vs %d\n", n.Data(), a->GetNbinsX(), b->GetNbinsX() ); continue; }
    double md = 0.;
    for ( int k = 1; k <= a->GetNbinsX(); ++k ) {
      double x = a->GetBinContent( k ), y = b->GetBinContent( k ), s = std::max( std::abs( x ), std::abs( y ) );
      if ( s > 0 ) md = std::max( md, std::abs( x - y ) / s );
    }
    printf( "[VERIFY] %-16s max rel diff %.2e %s\n", n.Data(), md, md < 1e-12 ? "SAME" : "DIFFERENT" );
  }
}
