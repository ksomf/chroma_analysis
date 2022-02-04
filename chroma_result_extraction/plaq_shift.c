//-------------------------------------------------------------------------------------------//
//-- Kim Somfleth <kim.somfleth@internode.on.net>                                          --//
//--                        LIME FILE EXTRACTER FOR FEYNMAN HELLMANN                       --//
//-------------------------------------------------------------------------------------------//

#include "../hl.h"
#include "chroma_lime_reader.h"
#include <stdio.h>

#define NumC 3.0

typedef struct Trajectory_Shift{
  c8  *trajectory;
  u64  dim[4];
#if defined( CLOVER_GMUNU )
  r64 *e_cross_b;
  r64 *b2_minus_e2;
#else
  r64  action;
  r64 *action_shifts;
  r64 *plaq_slices;
  r64  mean_plaquette;
  r64  mean_s_plaquette;
  r64  mean_t_plaquette;
  r64  mean_difference;
#endif
} Trajectory_Shift;

static void PrintMatrix3( r64 *a ){
  for( u64 i = 0; i < 3; ++i ){
    printf( "(%f,%f,%f)\n", a[i*3], a[i*3+1], a[i*3+2] );
  }
  printf( "\n" );
}

static void PrintMatrix3i( r64 *a_r, r64 *a_i ){
  for( u64 i = 0; i < 3; ++i ){
    printf( "(%f+%fi,%f+%fi,%f+%fi)\n", a_r[i*3], a_i[i*3], a_r[i*3+1], a_i[i*3+1], a_r[i*3+2], a_i[i*3+2] );
  }
  printf( "\n" );
}

//-- double U[NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]; <http://plone.jldg.org/wiki/images/f/f3/Ildg-file-format.pdf> --//
#define SWIZ_FLD 2
static m3c64 SwizzleMemory( hlMem_Pool *mem, r64 *data, u64 nx, u64 ny, u64 nz, u64 nt ){
  m3c64 result = hlM3C64Generate( mem, 4 * (nx+2*SWIZ_FLD)*(ny+2*SWIZ_FLD)*(nz+2*SWIZ_FLD)*(nt+2*SWIZ_FLD) );
  u64 tstride  = nz*ny*nx*4*3*3*2;
  u64 zstride  =    ny*nx*4*3*3*2;
  u64 ystride  =       nx*4*3*3*2;
  u64 xstride  =          4*3*3*2;
  u64 mustride =            3*3*2;
  u64 mu_s_stride = (nt+2*SWIZ_FLD)*(nz+2*SWIZ_FLD)*(ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
  u64 t_s_stride  =                 (nz+2*SWIZ_FLD)*(ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
  u64 z_s_stride  =                                 (ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
  u64 y_s_stride  =                                                 (nx+2*SWIZ_FLD);
  u64 x_s_stride  =                                                               1;
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( result.real.vals ); ++matrix_index ){
    for( u64 mu = 0; mu < 4; ++mu ){
      u64 mupos  = mu*mustride;
      u64 mulpos = mu*mu_s_stride;
      for( i64 it = -SWIZ_FLD; it < (i64)nt + SWIZ_FLD; ++it ){
        u64 tpos  = ((it+nt)%nt)*tstride  + mupos;
        u64 tlpos = (it+SWIZ_FLD)*t_s_stride + mulpos;
        for( i64 iz = -SWIZ_FLD; iz < (i64)nz + SWIZ_FLD; ++iz ){
          u64 zpos  = ((iz+nz)%nz)*zstride  + tpos;
          u64 zlpos = (iz+SWIZ_FLD)*z_s_stride + tlpos;
          for( i64 iy = -SWIZ_FLD; iy < (i64)ny + SWIZ_FLD; ++iy ){
            u64 ypos  = ((iy+ny)%ny)*ystride  + zpos;
            u64 ylpos = (iy+SWIZ_FLD)*y_s_stride + zlpos;
            for( i64 ix = -SWIZ_FLD; ix < (i64)nx + SWIZ_FLD; ++ix ){
              u64 base_index = ((ix+nx)%nx)*xstride  + ypos;
              u64 link_index = (ix+SWIZ_FLD)*x_s_stride + ylpos;
              result.real.vals[ matrix_index ][ link_index ] = data[ base_index + matrix_index * 2     ];
              result.imag.vals[ matrix_index ][ link_index ] = data[ base_index + matrix_index * 2 + 1 ];
            }
          }
        }
      }
    }
  }
 
  /*
//-- CHECK THE SWIZZLE WAS DONE CORRECTLY --//
#define GETFROMDATA( it, iz, iy, ix, mu ) (data + it*tstride + iz*zstride + iy*ystride + ix*xstride + mu*mustride)
#define GETFROMSWIZ( it, iz, iy, ix, mu ) ((it+SWIZ_FLD)*t_s_stride + (iz+SWIZ_FLD)*z_s_stride + (iy+SWIZ_FLD)*y_s_stride + (ix+SWIZ_FLD)*x_s_stride + mu*mu_s_stride)
  for( u64 mu = 0; mu < 4; ++mu ){
    for( u64 it = 0; it < nt; ++it ){
      for( u64 iz = 0; iz < nz; ++iz ){
        for( u64 iy = 0; iy < ny; ++iy ){
          for( u64 ix = 0; ix < nx; ++ix ){
            r64 *point = GETFROMDATA( it, iz, iy, ix, mu );
            u64 offset = GETFROMSWIZ( it, iz, iy, ix, mu );
            //printf( "%lu, %f, %f\n", offset, point[ 0 ], result.real.vals[ 0 ][ offset ] );
            hlASSERT( point[  0 ] == result.real.vals[ 0 ][ offset ] );
            hlASSERT( point[  1 ] == result.imag.vals[ 0 ][ offset ] );
            hlASSERT( point[  2 ] == result.real.vals[ 1 ][ offset ] );
            hlASSERT( point[  3 ] == result.imag.vals[ 1 ][ offset ] );
            hlASSERT( point[  4 ] == result.real.vals[ 2 ][ offset ] );
            hlASSERT( point[  5 ] == result.imag.vals[ 2 ][ offset ] );
            hlASSERT( point[  6 ] == result.real.vals[ 3 ][ offset ] );
            hlASSERT( point[  7 ] == result.imag.vals[ 3 ][ offset ] );
            hlASSERT( point[  8 ] == result.real.vals[ 4 ][ offset ] );
            hlASSERT( point[  9 ] == result.imag.vals[ 4 ][ offset ] );
            hlASSERT( point[ 10 ] == result.real.vals[ 5 ][ offset ] );
            hlASSERT( point[ 11 ] == result.imag.vals[ 5 ][ offset ] );
            hlASSERT( point[ 12 ] == result.real.vals[ 6 ][ offset ] );
            hlASSERT( point[ 13 ] == result.imag.vals[ 6 ][ offset ] );
            hlASSERT( point[ 14 ] == result.real.vals[ 7 ][ offset ] );
            hlASSERT( point[ 15 ] == result.imag.vals[ 7 ][ offset ] );
            hlASSERT( point[ 16 ] == result.real.vals[ 8 ][ offset ] );
            hlASSERT( point[ 17 ] == result.imag.vals[ 8 ][ offset ] );
          }
        }
      }
    }
  }
  */
  return( result );
}

i32 main( i32 argc, c8 **argv ){
  hlMem_Pool mem;
  hlMem_Pool tmp;
  hlInitializeProgramMemory( &mem, GIGABYTES( 4 ), &tmp, GIGABYTES( 16 ) );
  r64 beta = 6.0; //atof( argv[1] );
  printf( "USING BETA=%f\n", beta );

  u64 trajectory_count = argc-2;

  Trajectory_Shift *shifts = hlPUSHARRAY( &mem, trajectory_count, Trajectory_Shift );
  for( u64 file_index = argc-trajectory_count; file_index < argc; ++file_index ){
    hlResetMemoryPool( &tmp );
    c8 *path   = argv[ file_index ];
    hlFile file = hlFileRead( &tmp, path );
    u64 trajectory_index = file_index - (argc-trajectory_count);

    c8 *ident_str;
    if( hlStrEndsWithC( ".lime", path ) ){ //-- NEW STYLE qcdsf_b6p00_dp00_24x48.%d.%05d.lime --//
      c8 *filename = hlStrAfterLastC( '/', path );
      c8 *end_name = hlStrAfterLastC( 'l', filename ) - 2;
      ident_str = hlStrCopy( &mem, (u64)(end_name-filename), filename );
    }else{ //-- OLD STYLE GENERATED LIME FILE NAME WITH qcdsf_b6p00_dp00_24x48.lime%d --//
      c8 *filename = hlStrAfterLastC( '/', path );
      c8 *before_suffix = hlStrAfterLastC( '.', filename ) - 1;
      c8 *start_num     = hlStrAfterLastC( 'e', before_suffix );
      u64 config_num    = strtoul( start_num, 0, 10 );
          ident_str     = hlStrCopy( &mem, (u64)(before_suffix-filename) + 10, filename );
      snprintf( ident_str + (u64)(before_suffix-filename), 10, ".1.%05lu", config_num );
    }
    printf( "PARSING: %s", ident_str );
    Lime_Header *ildg_format_header = ReadLimeHeader( &file.contents );
    while( !hlStrStartsWithC( "ildg", ildg_format_header->type ) ){
      ildg_format_header = ReadLimeHeader( &file.contents );
    }
    hlXML format_xml = hlParseXML( &tmp, &mem, (c8 *)ildg_format_header->data, ildg_format_header->data_size-1 );
    u64 nx = hlXMLQueryNum( &tmp, &tmp, &format_xml, "lx" );
    u64 ny = hlXMLQueryNum( &tmp, &tmp, &format_xml, "ly" );
    u64 nz = hlXMLQueryNum( &tmp, &tmp, &format_xml, "lz" );
    u64 nt = hlXMLQueryNum( &tmp, &tmp, &format_xml, "lt" );
    u64 four_volume  = nx*ny*nz*nt;
    u64 three_volume = nx*ny*nz;

    Lime_Header *ildg_bin_header = ReadLimeHeader( &file.contents );
    printf( "\nREADING LIME FILE OF SIZE: %lu, THAT IS: %lu r64 OR PER SITE: %lu\n", ildg_bin_header->data_size, ildg_bin_header->data_size / sizeof(r64), ildg_bin_header->data_size / sizeof(r64) / nx / ny / nz / nt );
    hlASSERT( ildg_bin_header->data_size / sizeof(r64) == nx*ny*nz*nt*4*3*3*2 );
    BSwap64Array( (r64 *)ildg_bin_header->data, ildg_bin_header->data_size / sizeof(r64) );
    hlDEBUG_STR( "MEMORY READ\n" );
    m3c64 gauge_links = SwizzleMemory( &tmp, (r64 *)ildg_bin_header->data, nx, ny, nz, nt );
    hlDEBUG_STR( "MEMORY SWIZZLED\n" );
    u64 fmunu_count = 6*nt*nz*ny*nx;
#if defined( CLOVER_GMUNU )
    r64 *action_real   = hlPUSHARRAY( &tmp, fmunu_count, r64 );
    m3c64 clover_gmunu = hlM3C64Generate( &tmp, fmunu_count );
#else
    r64 *wilson_one_by_one = hlPUSHARRAY( &tmp, fmunu_count, r64 );
#endif

    m3c64 scratch1  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch2  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch3  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch4  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch5  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch6  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch7  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch8  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch9  = hlM3C64Generate( &tmp, nx );
    m3c64 scratch10 = hlM3C64Generate( &tmp, nx );
    m3c64 workspace = hlM3C64Generate( &tmp, nx );
    m3c64 adj_tmp   = hlM3C64Generate( &tmp, nx );

    u64 plane_stride  = nt*nz*ny*nx;
    u64 t_stride      =    nz*ny*nx;
    u64 z_stride      =       ny*nx;
    u64 y_stride      =          nx;
    u64 x_stride      =           1;
    u64 mu_s_stride = (nt+2*SWIZ_FLD)*(nz+2*SWIZ_FLD)*(ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
    u64 t_s_stride  =                 (nz+2*SWIZ_FLD)*(ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
    u64 z_s_stride  =                                 (ny+2*SWIZ_FLD)*(nx+2*SWIZ_FLD);
    u64 y_s_stride  =                                                 (nx+2*SWIZ_FLD);
    u64 x_s_stride  =                                                               1;
    u64 shift_factor[] = { x_s_stride, y_s_stride, z_s_stride, t_s_stride };         
    hlDEBUG_STR( "RUNNING PLAQUETTE CALCULATION\n" );
    for( u64 mu = 1, iplane = 0; mu < 4; ++mu ){
      u64 mu_shift = shift_factor[ mu ];
      for( u64 nu = 0; nu < mu; ++nu, ++iplane ){
        u64 nu_shift = shift_factor[ nu ];
        for( u64 it = 0; it < nt; ++it ){
          for( u64 iz = 0; iz < nz; ++iz ){
            for( u64 iy = 0; iy < ny; ++iy ){
              u64 st_swiz_stride = (it+SWIZ_FLD)*t_s_stride + (iz+SWIZ_FLD)*z_s_stride + (iy+SWIZ_FLD)*y_s_stride + SWIZ_FLD;
              u64 dest_stride    = iplane*plane_stride + it*t_stride + iz*z_stride + iy*y_stride;

              u64 in_mu = st_swiz_stride+mu*mu_s_stride;
              u64 in_nu = st_swiz_stride+nu*mu_s_stride;

#if defined( CLOVER_GMUNU )
              m3c64 dest_clover_plaq = hlM3C64CopyWithOffset( dest_stride, nx, &clover_gmunu );
              m3c64 u_mu      = hlM3C64CopyWithOffset( in_mu           , nx, &gauge_links );
              m3c64 u_nu_pmu  = hlM3C64CopyWithOffset( in_nu + mu_shift, nx, &gauge_links );
              m3c64 u_mu_pnu  = hlM3C64CopyWithOffset( in_mu + nu_shift, nx, &gauge_links );
              m3c64 u_nu      = hlM3C64CopyWithOffset( in_nu           , nx, &gauge_links );

              m3c64 u_mu_mmu_pnu = hlM3C64CopyWithOffset( in_mu - mu_shift + nu_shift, nx, &gauge_links );
              m3c64 u_nu_mmu     = hlM3C64CopyWithOffset( in_nu - mu_shift           , nx, &gauge_links );
              m3c64 u_mu_mmu     = hlM3C64CopyWithOffset( in_mu - mu_shift           , nx, &gauge_links );

              m3c64 u_nu_mmu_mnu = hlM3C64CopyWithOffset( in_nu - mu_shift - nu_shift, nx, &gauge_links );
              m3c64 u_mu_mmu_mnu = hlM3C64CopyWithOffset( in_mu - mu_shift - nu_shift, nx, &gauge_links );
              m3c64 u_nu_mnu     = hlM3C64CopyWithOffset( in_nu - nu_shift           , nx, &gauge_links );

              m3c64 u_mu_mnu     = hlM3C64CopyWithOffset( in_mu - nu_shift           , nx, &gauge_links );
              m3c64 u_nu_pmu_mnu = hlM3C64CopyWithOffset( in_nu + mu_shift - nu_shift, nx, &gauge_links );
              
              hlM3C64Mul(    &workspace, &scratch1, &u_mu    , &u_nu_pmu );
              hlM3C64Mul(    &workspace, &scratch2, &u_nu    , &u_mu_pnu );
              hlM3C64MulAdj( &workspace, &scratch3, &scratch1, &scratch2 );

              hlM3C64MulAdj( &workspace, &scratch1, &u_nu    , &u_mu_mmu_pnu );
              hlM3C64MulAdj( &workspace, &scratch2, &scratch1, &u_nu_mmu     );
              hlM3C64Mul(    &workspace, &scratch4, &scratch2, &u_mu_mmu     );
              //hlM3C64Mul(    &workspace, &scratch1, &u_mu_mmu, &u_nu         );
              //hlM3C64Mul(    &workspace, &scratch2, &u_nu_mmu, &u_mu_mmu_pnu );
              //hlM3C64MulAdj( &workspace, &scratch4, &scratch1, &scratch2     );

              hlM3C64Adj( &adj_tmp, &u_mu_mmu );
              hlM3C64MulAdj( &workspace, &scratch1, &adj_tmp , &u_nu_mmu_mnu );
              hlM3C64Mul(    &workspace, &scratch2, &scratch1, &u_mu_mmu_mnu );
              hlM3C64Mul(    &workspace, &scratch5, &scratch2, &u_nu_mnu     );
              //hlM3C64Mul(    &workspace, &scratch1, &u_mu_mmu_mnu, &u_nu_mnu );
              //hlM3C64Mul(    &workspace, &scratch2, &u_nu_mmu_mnu, &u_mu_mmu );
              //hlM3C64MulAdj( &workspace, &scratch5, &scratch1    , &scratch2 );

              hlM3C64Adj( &adj_tmp, &u_nu_mnu );
              hlM3C64Mul(    &workspace, &scratch1, &adj_tmp , &u_mu_mnu     );
              hlM3C64Mul(    &workspace, &scratch2, &scratch1, &u_nu_pmu_mnu );
              hlM3C64MulAdj( &workspace, &scratch6, &scratch2, &u_mu         );
              //hlM3C64Mul(    &workspace, &scratch1, &u_mu_mnu, &u_nu_pmu_mnu );
              //hlM3C64Mul(    &workspace, &scratch2, &u_nu_mnu, &u_mu         );
              //hlM3C64MulAdj( &workspace, &scratch6, &scratch1, &scratch2     );

              hlM3C64Add( &scratch7, &scratch3, &scratch4 );
              hlM3C64Add( &scratch8, &scratch5, &scratch6 );
              hlM3C64Add( &scratch9, &scratch7, &scratch8 );
              hlM3C64Adj( &scratch10, &scratch9 );
              hlM3C64Sub( &dest_clover_plaq, &scratch9, &scratch10 );
#else
              m3c64 u_mu     = hlM3C64CopyWithOffset( in_mu           , nx, &gauge_links );
              m3c64 u_nu_pmu = hlM3C64CopyWithOffset( in_nu + mu_shift, nx, &gauge_links );
              m3c64 u_nu     = hlM3C64CopyWithOffset( in_nu           , nx, &gauge_links );
              m3c64 u_mu_pnu = hlM3C64CopyWithOffset( in_mu + nu_shift, nx, &gauge_links );

              hlM3C64Mul( &workspace, &scratch1, &u_mu, &u_nu_pmu );
              hlM3C64Mul( &workspace, &scratch2, &u_nu, &u_mu_pnu );
              hlM3C64RealTraceMulAdj( wilson_one_by_one+dest_stride, &scratch1, &scratch2 );
#endif
            }
          }
        }
      }
    }
    hlDEBUG_STR( "PLAQUETTES CALCULATED\n" );


#if defined( CLOVER_GMUNU )
    hlDEBUG_STR( "NORMALISING AND SUBBING TRACE\n" );
    hlM3C64ScalerMulInplace( &clover_gmunu, 1.0 / 8.0 );
    hlM3C64SubTraceInplace( &clover_gmunu );
    hlDEBUG_STR( "DONE\n" );

    r64 *tmp_real1   = hlPUSHARRAY( &tmp, four_volume, r64 );
    r64 *tmp_real2   = hlPUSHARRAY( &tmp, four_volume, r64 );
    r64 *e_cross_b   = hlPUSHARRAY( &mem, four_volume, r64 );
    r64 *b2_minus_e2 = hlPUSHARRAY( &mem, four_volume, r64 );

    //printf( "\n" );
    //hlM3C64Print( &clover_gmunu, 0 );
    //hlM3C64Print( &clover_gmunu, 1 );
    //hlM3C64Print( &clover_gmunu, 2 );
    //hlM3C64Print( &clover_gmunu, 3 );
    //hlM3C64Print( &clover_gmunu, 4 );
    //hlM3C64Print( &clover_gmunu, 5 );
    //hlM3C64Print( &clover_gmunu, 6 );
    //hlM3C64Print( &clover_gmunu, 7 );

    //-- Iplane = (mu,nu), (2,1), (3,1), (3,2), (4,1), (4,2), (4,3) --//
    //-- ExB_z = -ReTr( iF41*iF31 + iF42*iF32 )                     --//
    //-- ExB_x =  ReTr( iF42*iF21 + iF43*iF31 )                     --//
    u64 plane_stride21 = 0*plane_stride;
    u64 plane_stride31 = 1*plane_stride;
    u64 plane_stride32 = 2*plane_stride;
    u64 plane_stride41 = 3*plane_stride;
    u64 plane_stride42 = 4*plane_stride;
    u64 plane_stride43 = 5*plane_stride;

    hlDEBUG_STR( "COPYING TENSORS\n" );
    m3c64 f21 = hlM3C64CopyWithOffset( plane_stride21, four_volume, &clover_gmunu );
    m3c64 f31 = hlM3C64CopyWithOffset( plane_stride31, four_volume, &clover_gmunu );
    m3c64 f32 = hlM3C64CopyWithOffset( plane_stride32, four_volume, &clover_gmunu );
    m3c64 f41 = hlM3C64CopyWithOffset( plane_stride41, four_volume, &clover_gmunu );
    m3c64 f42 = hlM3C64CopyWithOffset( plane_stride42, four_volume, &clover_gmunu );
    m3c64 f43 = hlM3C64CopyWithOffset( plane_stride43, four_volume, &clover_gmunu );
    
    hlDEBUG_STR( "CALCULATING OP2\n" );
    hlM3C64RealTraceMul( tmp_real1, &f42, &f21 );
    hlM3C64RealTraceMul( tmp_real2, &f43, &f31 );
    v64 minus_one = vsetpd( -1.0 );
    hlDEBUG_STR( "SUMMING CALCULATIONS TOGETHOR\n" );
    for( u64 i = 0; i < four_volume / vstep64; ++i ){
      v64 a = vmovapd( tmp_real1 + i*vstep64 );
      v64 b = vmovapd( tmp_real2 + i*vstep64 );
      //v64 res = vmulpd( vaddpd( a, b ), minus_one );
      v64 res = vaddpd( a, b );
      vmovapd( e_cross_b + i*vstep64, res );
    }
    hlDEBUG_STR( "OP2 DONE\n" );
    //printf( "Sample values: %f, %f, %f, %f %f\n", e_cross_b[0], e_cross_b[1], e_cross_b[2], e_cross_b[3], e_cross_b[4] );
    //-- B^2 - E^2 == ReTr( -(iF)^2temporal + 1/3*(iF)^2spatial ) --//
    hlDEBUG_STR( "CALCULATING ACTION\n" );
    hlM3C64RealTraceMul( action_real, &clover_gmunu, &clover_gmunu );
    r64 total      = 0.0;
    r64 check_total = 0.0;
    hlDEBUG_STR( "CALCULATION OP 2 BY SUMMATION OF PLAQUETTES\n" );
    for( u64 i = 0; i < four_volume; ++i ){
      b2_minus_e2[i] = 2.0 / 3.0 * ( action_real[ 0*plane_stride + i ]   // NOTE THIS HAS A MINUS SIGN FROM THE USUAL CONVETION SO BEWARE OF THAT 
                                   + action_real[ 1*plane_stride + i ]   // THIS IS BECAUSE WE USE THE ANTI HERMITIAN FIELD STRENGTH TENSOR 
                                   + action_real[ 2*plane_stride + i ] 
                                   - action_real[ 3*plane_stride + i ] 
                                   - action_real[ 4*plane_stride + i ] 
                                   - action_real[ 5*plane_stride + i ] );
      total += b2_minus_e2[i];
    }
    hlDEBUG_STR( "LOOKING AT ACTION DEBUG INFO\n" );
    for( u64 i = 0; i < four_volume * 6; ++i ){
      check_total += action_real[i];
    }
    printf( "\n" );
    for( u64 i = 0; i < 10; ++i ){
      printf( "%+f, ", b2_minus_e2[i] );
    }
    printf( "\n" );
    printf( "Volume Average: %f\n", total );
    printf( "Fmn Average: %f\n", check_total / four_volume / 6.0 );
    //action *= beta / 12.0 / four_volume / NumC;
    //should_be_zero *= beta / 12.0 / four_volume / NumC;
    Trajectory_Shift *shift = shifts + trajectory_index;
                      shift->dim[ 0 ]    = nx;
                      shift->dim[ 1 ]    = ny;
                      shift->dim[ 2 ]    = nz;
                      shift->dim[ 3 ]    = nt;
                      shift->trajectory  = ident_str;
                      shift->e_cross_b   = e_cross_b;
                      shift->b2_minus_e2 = b2_minus_e2;
#else
    r64 *action_shift = hlPUSHARRAY( &mem, four_volume, r64 );
    r64  t_plaq  = 0.0;
    r64  s_plaq  = 0.0;
    r64  total   = 0.0;
    for( u64 i = 0; i < four_volume; ++i ){
      action_shift[ i ] = 0.0;
      for( u64 iplane = 0; iplane < 6; ++iplane ){
        r64 tmp_plaq = wilson_one_by_one[ iplane*plane_stride + i ] / NumC;
        if( iplane < 3 ){
          action_shift[ i ] += beta * tmp_plaq;
          s_plaq += tmp_plaq;
        }else{
          action_shift[ i ] -= beta * tmp_plaq;
          t_plaq += tmp_plaq;
        }
      }
      total += action_shift[i];
    }
    r64 *w_slices = hlPUSHARRAY( &mem, nt, r64 );
    for( u64 it = 0; it < nt; ++it ){
      w_slices[ it ] = 0.0;
      for( u64 i = 0; i < three_volume; ++i ){
        u64 j = it*three_volume + i;
        for( u64 iplane = 0; iplane < 6; ++iplane ){
          w_slices[ it ] += wilson_one_by_one[ iplane*plane_stride + j ] / NumC;
        }
      }
      w_slices[ it ] /= three_volume * 6;
    }
    r64 action = beta * ( s_plaq - t_plaq );
    r64 mean_t_plaq = t_plaq / four_volume / 3.0;
    r64 mean_s_plaq = s_plaq / four_volume / 3.0;
    r64 mean_plaq = (mean_t_plaq+mean_s_plaq) / 2.0;
    Trajectory_Shift *shift = shifts + trajectory_index;
                      shift->dim[ 0 ] = nx;
                      shift->dim[ 1 ] = ny;
                      shift->dim[ 2 ] = nz;
                      shift->dim[ 3 ] = nt;
                      shift->mean_plaquette   = mean_plaq;
                      shift->mean_s_plaquette = mean_s_plaq;
                      shift->mean_t_plaquette = mean_t_plaq;
                      shift->action           = action;
                      shift->mean_difference  = mean_s_plaq - mean_t_plaq;
                      shift->trajectory       = ident_str;
                      shift->action_shifts    = action_shift;
                      shift->plaq_slices      = w_slices;
    r64 average_gluon_action = beta / 3.0 * ( 1.0 - mean_plaq  );
    printf( " %f,(%f,%f)\n", mean_plaq, mean_t_plaq, mean_s_plaq );
    printf( " %f\n", average_gluon_action );
    printf( "SAMPLPE WILSON PLAQ VALUES:" );
    for( u64 i = 0; i < 10; ++i ){
      printf( "%+f, ", shift->action_shifts[i] );
    }
    printf( "\n" );
    printf( "Total: %f\n", total );
   
/*
#if defined( IMP_WILSON_ACTION )
           tmp_plaq += (5.0/3.0) *( 1.0 - wilson_one_by_one[ stride + i ]/NumC;
                     - (1.0/12.0)*( 1.0 - wilson_one_by_two[ stride + i ]/NumC;
                     - (1.0/12.0)*( 1.0 - wilson_two_by_one[ stride + i ]/NumC;
           tmp_one_by_one += 1.0-wilson_one_by_one[ stride + i ]/NumC;
           tmp_one_by_two += 1.0-wilson_one_by_two[ stride + i ]/NumC;
           tmp_two_by_one += 1.0-wilson_two_by_one[ stride + i ]/NumC;
           printf( "->>%f,%f,%f\n", 1.0 - wilson_one_by_one[ stride + i ]/3  
           //                       , 1.0 - wilson_one_by_two[ stride + i ]/3
           //                       , 1.0 - wilson_two_by_one[ stride + i ]/3 );
#endif
*/
#endif
  }   

  printf( "CONSTRUCTING JSON\n" );
  hlJSON json = {};
  hlJSONGenerateObject( &mem, &json, 10 );
#if defined( CLOVER_GMUNU )
  hlJSONPushStringField( &json, "dat", "O[nconf][cloverB2mE2,cloverExB_x][nt][nz][ny][nx]" );
  hlJSON *action_ident = hlJSONPushArrayField( &mem, &json, "actions", 2 );
  hlJSONPushStringValue( action_ident, "cloverB2mE2" );
  hlJSONPushStringValue( action_ident, "cloverExB_x" );
#else
  r64 mean_plaq        = 0.0;
  r64 mean_t_plaq      = 0.0;
  r64 mean_s_plaq      = 0.0;
  r64 total_action     = 0.0;
  r64 mean_difference  = 0.0;
  hlJSONPushStringField( &json, "op", "wilsonB2mE2" );
  hlJSONPushStringField( &json, "dat", "action_shift[nconf][nt][nz][ny][nx]" );
#endif
  hlJSON *keys = hlJSONPushArrayField(  &mem, &json, "keys", trajectory_count );
  for( u64 trajectory_index = 0; trajectory_index < trajectory_count; ++trajectory_index ){
    Trajectory_Shift *shift = shifts + trajectory_index;
#if defined( CLOVER_GMUNU )
#else
    mean_plaq        += shift->mean_plaquette;
    mean_t_plaq      += shift->mean_t_plaquette;
    mean_s_plaq      += shift->mean_s_plaquette;
    mean_difference  += shift->mean_difference;
    total_action     += shift->action;
#endif
    hlJSONPushStringValue( keys, shift->trajectory );
  }
  hlJSON *dim = hlJSONPushArrayField( &mem, &json, "dim", hlARRAYCOUNT( shifts->dim ) );
  for( u64 dim_index = 0; dim_index < hlARRAYCOUNT( shifts->dim ); ++dim_index ){
    hlJSONPushNumberValue( dim, shifts->dim[ dim_index ] );
  }
#if defined( CLOVER_GMUNU )
#else
  total_action     /= trajectory_count;
  mean_plaq        /= trajectory_count;
  mean_t_plaq      /= trajectory_count;
  mean_s_plaq      /= trajectory_count;
  mean_difference  /= trajectory_count;
  
  hlJSONPushDoubleField( &json, "total_action" , total_action );
  hlJSONPushDoubleField( &json, "mean_plaq"       , mean_plaq        );
  hlJSONPushDoubleField( &json, "mean_t_plaq"     , mean_t_plaq      );
  hlJSONPushDoubleField( &json, "mean_s_plaq"     , mean_s_plaq      );
  hlJSONPushDoubleField( &json, "mean_difference" , mean_difference  );
  hlJSON *first_plaq_slices = hlJSONPushArrayField( &mem, &json, "first_slice", shifts->dim[3] );
  for( u64 it = 0; it < shifts->dim[3]; ++it ){
    hlJSONPushDoubleValue( first_plaq_slices, shifts->plaq_slices[it] );
  }
#endif

  c8 *filename = "action_shift.yml";
  c8 *data_filename = "action_shift.dat";
  printf( "Writing Files: %s,%s\n", filename, data_filename );
  FILE *out_file = fopen( filename, "w+" );
  hlJSONPrint( out_file, &json );
  fclose( out_file );

  FILE *data_file = fopen( data_filename, "wb" );
  for( u64 trajectory_index = 0; trajectory_index < trajectory_count; ++trajectory_index ){
    Trajectory_Shift *shift = shifts + trajectory_index;
#if defined( CLOVER_GMUNU )
    u64 four_volume = shift->dim[0]*shift->dim[1]*shift->dim[2]*shift->dim[3];
    fwrite( shift->b2_minus_e2, sizeof(r64), four_volume, data_file );
    fwrite( shift->e_cross_b  , sizeof(r64), four_volume, data_file );
#else
    fwrite( shift->action_shifts, sizeof(r64), shift->dim[0]*shift->dim[1]*shift->dim[2]*shift->dim[3], data_file );
#endif
  }
  fclose( data_file );
  hlPrintProgramMemory( &tmp, &mem );
  return( 0 );
}
