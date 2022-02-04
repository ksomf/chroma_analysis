//  g++ $( /home/eddie/jzanotti/opt/install/chroma-scalar-double/bin/chroma-config --cxxflags) num_to_mom.cpp -o num_to_mom
//  g++ $( /home/accounts/jzanotti/opt/install-chroma-double-openmpi/chroma-FH/bin/chroma-config --cxxflags) num_to_mom.cpp -o num_to_mom

#include <stdio.h>
#include "multi.h"
#include "sftmom.h"
#include "sftmom.cc"

using namespace qa;
int main( int argc, char **argv ){
  multi1d<int> mom_offset(4);
  mom_offset[0] = 0;
  mom_offset[1] = 0;
  mom_offset[2] = 0;
  mom_offset[3] = 0;
 
  FILE *file = fopen( "../chroma_momentum_numbers.h", "w+" );

  unsigned long max_mom2 = 17;
  for( int i = 0; i <= max_mom2; ++i ){
    SftMom test(i,mom_offset,false,3);
    fprintf( file, "i64 momentum_number_to_momenta_%d[][3] = {", i );
    for( int j = 0; j < test.numMom(); ++j ){
      multi1d<int> moms = test.numToMom( j );
      !j ? fprintf( file, "  " ) : fprintf( file, ", " );
      //IT TURNS OUT POSITIVE MOMENTUM MEANS NEGATIVE MOMENTA, CHROMA WHY
      fprintf( file, "{%2d,%2d,%2d}", -moms[0], -moms[1], -moms[2] );
    }
    fprintf( file, " };\n" );
  }
  fprintf( file, "i64 (*momentum_number_to_momenta[])[3] = {", max_mom2+1 );
  for( unsigned long i = 0; i <= max_mom2; ++i ){
    !i ? fprintf( file, "  " ) : fprintf( file, ", " );
    fprintf( file, "momentum_number_to_momenta_%lu", i );
  }
  fprintf( file, "};\n" );

  fprintf( file, "u64 mom2_max_to_momentum_number_count[%lu] = {", max_mom2+1 );
  for( unsigned long i = 0; i <= max_mom2; ++i ){
    !i ? fprintf( file, "  " ) : fprintf( file, ", " );
    fprintf( file, "hlARRAYCOUNT(momentum_number_to_momenta_%lu)", i );
  }
  fprintf( file, "};\n" );
  fclose( file );
  //unsigned long mom2_max = 13;
  //SftMom test(mom2_max,mom_offset,false,3);
  //for( int i = 0; i < test.numMom(); ++i ){
  //  multi1d<int> moms = test.numToMom( i );
  //  printf( "{%2d,%2d,%2d}\n", -moms[0], -moms[1], -moms[2] );
  //}
  //

  return( 0 );
}
