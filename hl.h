//-------------------------------------------------------------------------------------------//
//-- Kim Somfleth <kim.somfleth@internode.on.net>                                          --//
//--                      SINGLE HEADER LIBARY                                             --//
//-------------------------------------------------------------------------------------------//

//-- HL META WISHLIST --//
//-- PRINT STRUCT     --//
//-- ENUM2MAP         --//
//-- ->next structs   --//
#if !defined( _HARMLESS_H )
#include <stdlib.h>

#if defined( _MSC_VER )
#include <windows.h>
#elif defined( __GNUC__ )
#include <dlfcn.h> //-- CHECKED TO BE CRT-LESS SAFE --//
#endif

#if !(defined( hlAVX512 ) || defined( hlAVX ) || defined( hlSSE ))
  #if defined( __AVX512__ )
    #define hlAVX512
    #define hlSIMD
  #elif defined( __AVX__ )
    #define hlAVX
    #define hlSIMD
  #endif
#endif
#if defined( __cplusplus )
  #define hlEXTERN extern "C" 
#else
  #define hlEXTERN extern
#endif

#define _STR( str ) #str
#define STR( str ) _STR( str )

#define hlUSEWITHCAUTION(_1,_2,name,...) name

#if defined( _MSC_VER )
hlEXTERN int _fltused = 0; //-- WTF C --//
#elif defined( __GNUC__ )
hlEXTERN int _fltused;     //-- WTF C --//
#endif


//-- THESE ARE CRT-LESS SAFE HEADER FILES --//
#include <stddef.h> //-- size_t                                 --//
#include <stdint.h> //-- int8_t -- int64_t, uint8_t -- u64_t --//
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef char     c8;
typedef u32      b32;
typedef float    r32;
typedef double   r64;
#include <stdarg.h> //-- va_arg, va_start, va_end, va_arg       --//
#if defined( _MSC_VER )
  #include <Intrin.h>      //-- ASSEMBLY INTRINSICS --//
#elif defined( __GNUC__ )
  r64 strtod( const c8 *str, c8 **end );
  long atol( const c8 *str );
  hlEXTERN r64 atof( const c8 *str );
  //#include <x86intrin.h>   //-- ASSEMBLY INTRINSICS --//
#endif

#define hlARRAYCOUNT( array ) ( sizeof( array ) / sizeof( (array)[ 0 ] ) )

//-------------------------------------------------------------------------------------------//
//--                                     ASSERTS                                           --//
//-------------------------------------------------------------------------------------------//

#if defined( hlDEBUG )
  #define __hlSTATIC_ASSERT3( expr, msg ) typedef char static_assertion_##msg[(!!( expr ) ) * 2 - 1]
  #define __hlSTATIC_ASSERT2( expr, ln ) __hlSTATIC_ASSERT3( expr, cc_assert_line_##ln )
  #define __hlSTATIC_ASSERT1( expr, ln ) __hlSTATIC_ASSERT2( expr, ln )
  #define hlSTATIC_ASSERT( expression )  __hlSTATIC_ASSERT1( expression, __LINE__ )
#if defined( _MSC_VER )
  #define hlWRITE_DEBUG_STR( str ) ( OutputDebugStringA( str ) )
  #define hlASSERT( expression ) if( !(expression ) ){ *(i32 *)0 = 0; }
#elif defined( __GNUC__ )
#if defined( __cplusplus ) 
#include <stdio.h> //--TODO--//
  #define hlWRITE_DEBUG_STR( str )( { char write_str[] = str; printf( "%s", write_str ); } )
  #define hlDEBUG_STR( str ) hlWRITE_DEBUG_STR( "[DEBUG] " str )
  #define hlERROR_STR( str ) hlWRITE_DEBUG_STR( "[ERROR] " str )
#else
  u64 write( i32 dest, void *buffer, u64 buffer_size ); //-- INCLUDES UNISTR --//
  #define hlWRITE_DEBUG_STR( str ) ( { char write_str[] = str; write( 1, write_str, sizeof( write_str ) ); } )
  #define hlDEBUG_STR( str ) hlWRITE_DEBUG_STR( "[DEBUG] " str )
  #define hlERROR_STR( str ) hlWRITE_DEBUG_STR( "[ERROR] " str )
#endif
  #define _hlASSERTION_STRING2( file, line ) "ASSERTION FAILED! " file ":" #line "\n"
  #define _hlASSERTION_STRING1( file, line ) _hlASSERTION_STRING2( file, line ) 
  #define ASSERTION_STRING _hlASSERTION_STRING1( __FILE__, __LINE__ )

  #define hlASSERT( expression ) if( !(expression) ){ hlERROR_STR( ASSERTION_STRING ); *(volatile i32 *)0 = 0; }
#endif

#else
  #define hlDEBUG_STR( str )
  #define hlERROR_STR( str )
  #define hlSTATIC_ASSERT( expression )
  #define hlASSERT( expression )
#endif

#define hlCOUNTER_DEBUG (hlDEBUG_STR( "COUNTER" STR( __COUNTER__ ) "\n" ))

hlSTATIC_ASSERT( sizeof( u32    ) == 4 );
hlSTATIC_ASSERT( sizeof( u64    ) == 8 );
hlSTATIC_ASSERT( sizeof( i32    ) == 4 );
hlSTATIC_ASSERT( sizeof( i64    ) == 8 );
hlSTATIC_ASSERT( sizeof( size_t ) == 8 );
hlSTATIC_ASSERT( sizeof( r32    ) == 4 );
hlSTATIC_ASSERT( sizeof( r64    ) == 8 );
hlSTATIC_ASSERT( sizeof( void * ) == 8 );
hlSTATIC_ASSERT( sizeof( c8     ) == 1 );

#if !(defined( __x86_64__ ) || defined( _M_AMD64 ))
#error REQUIRE 64 BIT OS
#endif

#define R64_DIGIT_PRECISION 15                      //-- NUMBER OF DIGITS IN PRECISION --//
#define R64_DIGIT_MANTISSA  53                      //-- NUMBER OF DIGITS IN MANTISSA  --//
#define R64_EPSILON         2.2204460492503131E-16  //-- 1.0 + REAL64_EPSILON != 1.0   --//
#define R64_MIN             2.2250738585072014E-308 //-- MINIMUM REAL64                --//
#define R64_MAX             1.7976931348623157E+308 //-- MAXIMUM REAL64                --//

#define R32_LARGEST_ODD     16777215                //-- POINT AT WHICH (u)int32 MORE PRECISE THAN real32 --//
#define R32_DIGIT_PRECISION 6                       //-- NUMBER OF DIGITS IN PRECISION                    --//
#define R32_DIGIT_MANTISSA  24                      //-- NUMBER OF DIGITS IN MANTISSA                     --//
#define R32_EPSILON         1.19209290E-07f         //-- 1.0 + REAL32_EPSILON != 1.0                      --//
#define R32_MIN             1.17549435E-38f         //-- MINIMUM REAL32                                   --//
#define R32_MAX             3.40282347E+38f         //-- MAXIMUM REAL32                                   --//

#define KILOBYTES( a ) (( 1024LL *            a   ))
#define MEGABYTES( a ) (( 1024LL * KILOBYTES( a ) ))
#define GIGABYTES( a ) (( 1024LL * MEGABYTES( a ) ))

//-------------------------------------------------------------------------------------------//
//--                                      INTRINSICS                                       --//
//-------------------------------------------------------------------------------------------//
#if defined( __clang__ )
#define bswap32 __builtin_bswap32
#define bswap64 __builtin_bswap64
#else
#define bswap32 __builtin_bswap32 //WHATTODO
#define bswap64 __builtin_bswap64
#endif
static u16 bswap16( u16 data ){ return ( ( ( data & 0x00ff ) << 8 ) | ( ( data & 0xff00 ) >> 8  ) ); }
static void BSwap64Array( void *arr, u64 count ){
  i64 *array = (i64 *)arr;
  for( u64 index = 0; index < count; ++index ){
    array[ index ] = bswap64( array[ index ] );
  }
}

#define max_vstep32 16
#define max_vstep64 8
#if defined( hlAVX512 )
  #include <immintrin.h>
  #define vstep32 16
  #define vstep64 8
  #define v32          __m512
  #define v64          __m512d
  #define v32i         __m512i
  #define vmulps       _mm512_mul_ps
  #define vmulpd       _mm512_mul_pd
  #define vaddps       _mm512_add_ps
  #define vaddpd       _mm512_add_pd
  #define vsubps       _mm512_sub_ps
  #define vsubpd       _mm512_sub_pd
  #define vmovapd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm512_store_pd , _mm512_load_pd  )( __VA_ARGS__ )
  #define vmovupd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm512_storeu_pd, _mm512_loadu_pd )( __VA_ARGS__ )
  #define vxorpd       _mm512_xor_pd
  #define vxorps       _mm512_xor_ps
  #define vzeropd      _mm512_setzero_pd
  #define vzerops      _mm512_setzero_ps
  #define vsetpd       _mm512_set1_pd
  #define vsetps       _mm512_set1_ps
  #define vblendpd     _mm512_blend_pd
  #define vblendps     _mm512_blend_ps
#elif defined( hlAVX )
  #define vstep32 8
  #define vstep64 4
  #define v32          __m256
  #define v64          __m256d
  #define v32i         __m256i
  #define vmulps       _mm256_mul_ps
  #define vmulpd       _mm256_mul_pd
  #define vaddps       _mm256_add_ps
  #define vaddpd       _mm256_add_pd
  #define vsubps       _mm256_sub_ps
  #define vsubpd       _mm256_sub_pd
  #define vmovapd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm256_store_pd , _mm256_load_pd  )( __VA_ARGS__ )
  #define vmovupd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm256_storeu_pd, _mm256_loadu_pd )( __VA_ARGS__ )
  #define vxorpd       _mm256_xor_pd
  #define vxorps       _mm256_xor_ps
  #define vzeropd      _mm256_setzero_pd
  #define vzerops      _mm256_setzero_ps
  #define vsetpd       _mm256_set1_pd
  #define vsetps       _mm256_set1_ps
  #define vblendpd     _mm256_blend_pd
  #define vblendps     _mm256_blend_ps
#elif defined( hlSSE )
  #include <xmmintrin.h> //-- SSE  --//
  //#include <emmintrin.h> //-- SSE2 --//
  //#include <pmmintrin.h> //-- SSE3 --//
  //#include <tmmintrin.h> //-- SSE4 --//
#if defined( _MSC_VER )
  //#include <nmmintrin.h> //-- SSE4.2 FOR WINDOWS --//
#endif
  
  #define vstep32 4
  #define vstep64 2
  #define v32          __m128
  #define v64          __m128d
  #define v32i         __m128i
  #define vmulps       _mm_mul_ps
  #define vmulpd       _mm_mul_pd
  #define vaddps       _mm_add_ps
  #define vaddpd       _mm_add_pd
  #define vsubps       _mm_sub_ps
  #define vsubpd       _mm_sub_pd
  #define vmovapd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm_store_pd , _mm_load_pd  )( __VA_ARGS__ )
  #define vmovupd(...) hlUSEWITHCAUTION( __VA_ARGS__, _mm_storeu_pd, _mm_loadu_pd )( __VA_ARGS__ )
  #define vxorpd       _mm_xor_pd
  #define vxorps       _mm_xor_ps
  #define vzeropd      _mm_setzero_pd
  #define vzerops      _mm_setzero_ps
  #define vsetpd       _mm_set1_pd
  #define vsetps       _mm_set1_ps
  #define vblendpd     _mm_blend_pd
  #define vblendps     _mm_blend_ps
#endif


//-------------------------------------------------------------------------------------------//
//--                                     MEM POOLS                                         --//
//-------------------------------------------------------------------------------------------//

typedef struct {
  u8   *base;
  u64   size;
  u64   used;
#if defined( hlDEBUG )
  i32   count;
#endif
} hlMem_Pool;

static void hlInitializeMemoryPool( hlMem_Pool *pool, u8 *base, u64 size ){
  pool->size  = size;
  pool->base  = base;
  pool->used  = 0;
#if defined( hlDEBUG )
  pool->count = 0;
#endif
}

static void hlResetMemoryPool( hlMem_Pool *pool ){
  pool->used = 0;
}

static void hlZeroMemory( void *memory, u64 size ){
#if 0
  u64 ywords_to_zero = size / sizeof( __m256 );
  u64 bytes_to_zero  = size % sizeof( __m256 );

  u8 *mem = (u8 *)memory;
  v64 zero = vzerops();
  
  for( u64 yword_index = 0; yword_index < ywords_to_zero; ++yword_index ){
    vmovups( (r32 *)mem, zero );
    mem = mem + sizeof( __m256 );
  }

  for( u64 byte_index = 0; byte_index < bytes_to_zero; ++byte_index ){
    *mem++ = 0;
  }
#else
  u8 *value = (u8 *)memory;
  for( u64 i = 0; i < size; ++i ){
    *value++ = 0;
  }
#endif
}

static b32 hlMemCmp( void *mem1, void *mem2, u64 size ){
  b32 result = 1;
  u8 *val1 = (u8 *)mem1;
  u8 *val2 = (u8 *)mem2;
  for( u64 i = 0; result && i < size; ++i ){
    result &= *val1++ == *val2++;
  }
  return(result);
}

#define ALIGN_U8_PTR_TO_8_BYTES( ptr ) ( (u8 *)( (u64)(ptr + 0x7) & ~0x7 ) )
#define ALIGN_VOID_PTR_TO_32_BYTES( ptr ) ( (void *)( (u64)(ptr + 0x1f) & ~0x1f ) )
#define hlPUSHARRAY( pool, count, type ) ((type *)_hlPushSize( pool, (sizeof( type ) * count ) ))
#define hlPUSHTYPE(  pool,        type ) ((type *)_hlPushSize( pool, sizeof( type ) ))
#define hlPUSHSIZE(  pool,        size ) (        _hlPushSize( pool, size ))
static void *_hlPushSize( hlMem_Pool *pool, u64 size ){
  size = (size + 0x1f) & ~0x1f;
  hlASSERT( ALIGN_VOID_PTR_TO_32_BYTES( pool->base + pool->used ) == pool->base + pool->used );
  //printf( "POOL DEBUG (%lu) %lu / %lu\n", size, pool->used, pool->size );
  hlASSERT( ( pool->used + size ) <= pool->size );
  void *result = pool->base + pool->used; 
  hlZeroMemory( result, size );
  
  pool->used += size;

  return( result );
}

void *memcpy(void *dest, const void *src, unsigned long bytes ); //--TODO--//
static void hlMemCopy( void *src, void *dest, u64 bytes ){
  memcpy( dest, src, bytes );
}

//-------------------------------------------------------------------------------------------//
//--                               CONVENIENCE LIBRARY                                     --//
//-------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------//
//--                                    MATH LIBRARY                                       --//
//-------------------------------------------------------------------------------------------//

b32 hlIsnan( r64 val ){
  return( val != val );
}

typedef struct hlDbl_Bits {
  union{
    r64 dbl;
    i64 bits;
  };
} hlDbl_Bits;


static i64 hlAbsi64( i64 val ){
  i64 result = val < 0 ? -val : val;
  return result;
}

static r32 hlAbs32( r32 val ){
  r32 result = val < 0.0f ? -val : val;
  return( result );
}

static r64 hlAbs64( r64 val ){
  r64 result = val < 0.0f ? -val : val;
  return( result );
}

typedef struct {
  union {
    r64 *vals[9];
    struct {
      r64 *row[3][3];
    };
    struct {
      r64 *row0[3];
      r64 *row1[3];
      r64 *row2[3];
    };
    struct {
      r64 *r11;
      r64 *r12;
      r64 *r13;
      r64 *r21;
      r64 *r22;
      r64 *r23;
      r64 *r31;
      r64 *r32;
      r64 *r33;
    };
  };
  u64 count;
} m3r64;

typedef struct {
  m3r64 real;
  m3r64 imag;
} m3c64;

#include <stdio.h>
static void hlM3C64Print( m3c64 *a, u64 count ){
  for( u64 row_index = 0; row_index < hlARRAYCOUNT( a->real.row ); ++row_index ){
    printf( "(%+f %+fi,%+f %+fi,%+f %+fi)\n", a->real.row[row_index][ 0 ][ count ], a->imag.row[row_index][ 0 ][ count ]
                                            , a->real.row[row_index][ 1 ][ count ], a->imag.row[row_index][ 1 ][ count ]
                                            , a->real.row[row_index][ 2 ][ count ], a->imag.row[row_index][ 2 ][ count ] );
  }
  printf( "\n" );
}

static m3r64 hlM3R64Generate( hlMem_Pool *pool, u64 count ){
  m3r64 result;
        result.count = count;
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( result.vals ); ++matrix_index ){
    result.vals[ matrix_index ] = hlPUSHARRAY( pool, count, r64 );
  }
  return( result );
}

static m3c64 hlM3C64Generate( hlMem_Pool *pool, u64 count ){
  m3c64 result;
        result.real = hlM3R64Generate( pool, count );
        result.imag = hlM3R64Generate( pool, count );
  return( result );
}

static m3r64 hlM3R64CopyWithOffset( u64 offset, u64 count, m3r64 *base ){
  m3r64 result;
        result.count = count;
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( result.vals ); ++matrix_index ){
    result.vals[ matrix_index ] = base->vals[ matrix_index ] + offset;
  }
  return( result );
}

static m3c64 hlM3C64CopyWithOffset( u64 offset, u64 count, m3c64 *base ){
  m3c64 result;
        result.real = hlM3R64CopyWithOffset( offset, count, &base->real );
        result.imag = hlM3R64CopyWithOffset( offset, count, &base->imag );
  return( result );
}

static void hlM3R64Add( m3r64 *dest, m3r64 *a, m3r64 *b ){
  hlASSERT( dest->count % max_vstep64 == 0 );
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->vals ); ++matrix_index ){
#if defined( hlSIMD )
    for( u64 i = 0; i < dest->count / vstep64; ++i ){
       v64 ar = vmovapd( a->vals[ matrix_index ] + i*vstep64 );
       v64 br = vmovapd( b->vals[ matrix_index ] + i*vstep64 ); 
       v64 res = vaddpd( ar, br );
       vmovapd( dest->vals[ matrix_index ] + i*vstep64, res );
    }
#else
    for( u64 i = 0; i < dest->count; ++i ){
      dest->vals[ matrix_index ][ i ] = a->vals[ matrix_index ][ i ] + b->vals[ matrix_index ][ i ];
    }
#endif
  }
}

static void hlM3C64Add( m3c64 *dest, m3c64 *a, m3c64 *b ){
  hlM3R64Add( &dest->real, &a->real, &b->real );
  hlM3R64Add( &dest->imag, &a->imag, &b->imag );
}

static void hlM3R64Sub( m3r64 *dest, m3r64 *a, m3r64 *b ){
  hlASSERT( dest->count % max_vstep64 == 0 );
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->vals ); ++matrix_index ){
#if defined( hlSIMD )
    for( u64 i = 0; i < dest->count / vstep64; ++i ){
       v64 ar = vmovapd( a->vals[ matrix_index ] + i*vstep64 );
       v64 br = vmovapd( b->vals[ matrix_index ] + i*vstep64 ); 
       v64 res = vsubpd( ar, br );
       vmovapd( dest->vals[ matrix_index ] + i*vstep64, res );
    }
#else
    for( u64 i = 0; i < dest->count; ++i ){
      dest->vals[ matrix_index ][ i ] = a->vals[ matrix_index ][ i ] - b->vals[ matrix_index ][ i ];
    }
#endif
  }
}

static void hlM3C64Sub( m3c64 *dest, m3c64 *a, m3c64 *b ){
  hlM3R64Sub( &dest->real, &a->real, &b->real );
  hlM3R64Sub( &dest->imag, &a->imag, &b->imag );
}

static void hlM3R64ScalarMulInplace( m3r64 *a, r64 b ){
  v64 factor = vsetpd( b );
  for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->vals ); ++matrix_index ){
    for( u64 i = 0; i < a->count / vstep64; ++i ){
      v64 res = vmovapd( a->vals[ matrix_index ] + i*vstep64 );
      vmovapd( a->vals[ matrix_index ] + i*vstep64, vmulpd( res, factor ) );
    }
  }
}

static void hlM3C64ScalerMulInplace( m3c64 *a, r64 b ){
  hlM3R64ScalarMulInplace( &a->real, b );
  hlM3R64ScalarMulInplace( &a->imag, b );
}

static void hlM3R64Adj( m3r64 *dest, m3r64 *a ){
#if defined( hlSIMD )
  for( u64 row_index = 0; row_index < hlARRAYCOUNT( dest->row ); ++row_index ){
    for( u64 col_index = 0; col_index < hlARRAYCOUNT( dest->row0 ); ++col_index ){
       for( u64 i = 0; i < dest->count / vstep64; ++i ){
         v64 res = vmovupd( a->row[ row_index ][ col_index ] + i*vstep64 );
         vmovapd( dest->row[ col_index ][ row_index ] + i*vstep64, res );
       }
    }
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED" );
#endif
}

static void hlM3C64Adj( m3c64 *dest, m3c64 *a ){
  hlM3R64Adj( &dest->real, &a->real );
  hlM3R64Adj( &dest->imag, &a->imag );
  hlM3R64ScalarMulInplace( &dest->imag, -1.0 );
}

static void hlM3R64AdjInplace( m3r64 *a ){
  r64 *r10 = a->row[ 1 ][ 0 ];
  r64 *r01 = a->row[ 0 ][ 1 ];
  r64 *r20 = a->row[ 2 ][ 0 ];
  r64 *r02 = a->row[ 0 ][ 2 ];
  r64 *r12 = a->row[ 1 ][ 2 ];
  r64 *r21 = a->row[ 2 ][ 1 ];

  a->row[ 0 ][ 1 ] = r10;
  a->row[ 1 ][ 0 ] = r01;
  a->row[ 0 ][ 2 ] = r20;
  a->row[ 2 ][ 0 ] = r02;
  a->row[ 2 ][ 1 ] = r12;
  a->row[ 1 ][ 2 ] = r21;
}

static void hlM3R64Mul( m3r64 *dest, m3r64 *a, m3r64 *b ){
  hlASSERT( dest->count % max_vstep64 == 0 );
  for( u64 row_index = 0; row_index < hlARRAYCOUNT( dest->row ); ++row_index ){
    for( u64 col_index = 0; col_index < hlARRAYCOUNT( dest->row0 ); ++col_index ){
#if defined( hlSIMD )
      for( u64 i = 0; i < dest->count / vstep64; ++i ){
        v64 res = vzeropd();
        for( u64 inte_index = 0; inte_index < hlARRAYCOUNT( dest->row0 ); ++inte_index ){
          v64 ar = vmovupd( a->row[ row_index  ][ inte_index ] + i*vstep64 );
          v64 br = vmovupd( b->row[ inte_index ][ col_index  ] + i*vstep64 );
          res = vaddpd( res, vmulpd( ar, br ) );
        }
        vmovapd( dest->row[ row_index ][ col_index ] + i*vstep64, res );
      }
#else
      for( u64 i = 0; i < dest->count; ++i ){
        dest->row[ row_index ][ col_index ][ i ] = 0.0;
        for( u64 inte_index = 0; inte_index < hlARRAYCOUNT( dest->row0 ); ++inte_index ){
          dest->row[ row_index ][ col_index ][ i ] += a->row[ row_index  ][ inte_index ][ i ] 
                                                    * b->row[ inte_index ][ col_index  ][ i ];
        }
      }
#endif
    }
  }
}

static void hlM3R64MulAdj( m3r64 *dest, m3r64 *a, m3r64 *b ){
  hlASSERT( dest->count % max_vstep64 == 0 );
  for( u64 row_index = 0; row_index < hlARRAYCOUNT( dest->row ); ++row_index ){
    for( u64 col_index = 0; col_index < hlARRAYCOUNT( dest->row0 ); ++col_index ){
#if defined( hlSIMD )
      for( u64 i = 0; i < dest->count / vstep64; ++i ){
        v64 res = vzeropd();
        for( u64 inte_index = 0; inte_index < hlARRAYCOUNT( dest->row0 ); ++inte_index ){
          v64 ar = vmovupd( a->row[ row_index ][ inte_index ] + i*vstep64 );
          v64 br = vmovupd( b->row[ col_index ][ inte_index ] + i*vstep64 );
          res = vaddpd( res, vmulpd( ar, br ) );
        }
        vmovapd( dest->row[ row_index ][ col_index ] + i*vstep64, res );
      }
#else
      for( u64 i = 0; i < dest->count; ++i ){
        dest->row[ row_index ][ col_index ][ i ] = 0.0;
        for( u64 inte_index = 0; inte_index < hlARRAYCOUNT( dest->row0 ); ++inte_index ){
          dest->row[ row_index ][ col_index ][ i ] += a->row[ row_index ][ inte_index ][ i ] 
                                                    * b->row[ col_index ][ inte_index ][ i ];
        }
      }
#endif
    }
  }
}
static void hlM3R64Trace( r64 *dest, m3r64 *a ){
  hlASSERT( a->count % max_vstep64 == 0 );
#if defined( hlSIMD )
  for( u64 i = 0; i < a->count / vstep64; ++i ){
    v64 r0c0 = vmovapd( a->row0[0] + i*vstep64 );
    v64 r1c1 = vmovapd( a->row1[1] + i*vstep64 );
    v64 r2c2 = vmovapd( a->row2[2] + i*vstep64 );
    v64 res  = vaddpd( r0c0, vaddpd( r1c1, r2c2 ) );
    vmovapd( dest + i*vstep64, res );
  }
#else
  for( u64 i = 0; i < a->count; ++i ){
    dest[ i ] = a->row0[0][ i ] + a->row1[1][ i ] + a->row2[2][ i ];
  }
#endif
}

static void hlM3C64Mul( m3c64 *scratch, m3c64 *dest, m3c64 *a, m3c64 *b ){
  hlM3R64Mul( &dest->real   , &a->real   , &b->real       );
  hlM3R64Mul( &scratch->real, &a->imag   , &b->imag       );
  hlM3R64Sub( &dest->real   , &dest->real, &scratch->real );

  hlM3R64Mul( &dest->imag   , &a->imag   , &b->real       ); 
  hlM3R64Mul( &scratch->imag, &a->real   , &b->imag       ); 
  hlM3R64Add( &dest->imag   , &dest->imag, &scratch->imag );
}

static void hlM3C64MulAdj( m3c64 *scratch, m3c64 *dest, m3c64 *a, m3c64 *b ){
  hlM3R64MulAdj( &dest->real   , &a->real   , &b->real       );
  hlM3R64MulAdj( &scratch->real, &a->imag   , &b->imag       );
  hlM3R64Add(    &dest->real   , &dest->real, &scratch->real );

  hlM3R64MulAdj( &dest->imag   , &a->imag   , &b->real       ); 
  hlM3R64MulAdj( &scratch->imag, &a->real   , &b->imag       ); 
  hlM3R64Sub(    &dest->imag   , &dest->imag, &scratch->imag );
}

static void hlM3R64SubTraceInplace( m3r64 *a ){
  hlASSERT( a->count % max_vstep64 == 0 );
  v64 third = vsetpd( 1.0 / 3.0 );
  for( u64 i = 0; i < a->count / vstep64; ++i ){
    v64 r0c0 = vmovapd( a->row0[0] + i*vstep64 );
    v64 r1c1 = vmovapd( a->row1[1] + i*vstep64 );
    v64 r2c2 = vmovapd( a->row2[2] + i*vstep64 );
    v64 trace  = vmulpd( third, vaddpd( r0c0, vaddpd( r1c1, r2c2 ) ) );
    vmovapd( a->row0[0] + i*vstep64, vsubpd( r0c0, trace ) );
    vmovapd( a->row1[1] + i*vstep64, vsubpd( r1c1, trace ) );
    vmovapd( a->row2[2] + i*vstep64, vsubpd( r2c2, trace ) );
  }
}

static void hlM3C64SubTraceInplace( m3c64 *a ){
  hlM3R64SubTraceInplace( &a->real );
  hlM3R64SubTraceInplace( &a->imag );
}

static void hlM3C64RealTraceMul( r64 *destr, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 resr = vzeropd();
    for( u64 row_index = 0; row_index < hlARRAYCOUNT( a->real.row ); ++row_index ){
      for( u64 col_index = 0; col_index < hlARRAYCOUNT( a->real.row0 ); ++col_index ){
        v64 ar = vmovapd( a->real.row[ row_index ][ col_index ] + i*vstep64 );
        v64 br = vmovapd( b->real.row[ col_index ][ row_index ] + i*vstep64 );
        v64 ai = vmovapd( a->imag.row[ row_index ][ col_index ] + i*vstep64 );
        v64 bi = vmovapd( b->imag.row[ col_index ][ row_index ] + i*vstep64 );
        resr =  vaddpd( resr, vsubpd( vmulpd( ar, br ), vmulpd( ai, bi ) ) );
      }
    }        
    vmovapd( destr + i*vstep64, resr );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}

static void hlM3C64TraceMul( r64 *destr, r64 *desti, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 resr = vzeropd();
    v64 resi = vzeropd();
    for( u64 row_index = 0; row_index < hlARRAYCOUNT( a->real.row ); ++row_index ){
      for( u64 col_index = 0; col_index < hlARRAYCOUNT( a->real.row0 ); ++col_index ){
        v64 ar = vmovapd( a->real.row[ row_index ][ col_index ] + i*vstep64 );
        v64 br = vmovapd( b->real.row[ col_index ][ row_index ] + i*vstep64 );
        v64 ai = vmovapd( a->imag.row[ row_index ][ col_index ] + i*vstep64 );
        v64 bi = vmovapd( b->imag.row[ col_index ][ row_index ] + i*vstep64 );
        resr =  vaddpd( resr, vsubpd( vmulpd( ar, br ), vmulpd( ai, bi ) ) );
        resi =  vaddpd( resi, vaddpd( vmulpd( ar, bi ), vmulpd( ai, br ) ) );
      }
    }        
    vmovapd( destr + i*vstep64, resr );
    vmovapd( desti + i*vstep64, resi );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}

static void hlM3C64TraceMulAdj( r64 *destr, r64 *desti, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 resr = vzeropd();
    v64 resi = vzeropd();
    for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->real.vals ); ++matrix_index ){
      v64 ar = vmovapd( a->real.vals[ matrix_index ] + i*vstep64 );
      v64 br = vmovapd( b->real.vals[ matrix_index ] + i*vstep64 );
      v64 ai = vmovapd( a->imag.vals[ matrix_index ] + i*vstep64 );
      v64 bi = vmovapd( b->imag.vals[ matrix_index ] + i*vstep64 );
      resr =  vaddpd( resr, vaddpd( vmulpd( ar, br ), vmulpd( ai, bi ) ) );
      resi =  vsubpd( resi, vaddpd( vmulpd( ar, bi ), vmulpd( ai, br ) ) );
    }        
    vmovapd( destr + i*vstep64, resr );
    vmovapd( desti + i*vstep64, resi );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}

static void hlM3C64RealTraceMulAdj( r64 *dest, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 res = vzeropd();
    for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->real.vals ); ++matrix_index ){
      v64 ar = vmovapd( a->real.vals[ matrix_index ] + i*vstep64 );
      v64 br = vmovapd( b->real.vals[ matrix_index ] + i*vstep64 );
      v64 ai = vmovapd( a->imag.vals[ matrix_index ] + i*vstep64 );
      v64 bi = vmovapd( b->imag.vals[ matrix_index ] + i*vstep64 );
      res =  vaddpd( res, vaddpd( vmulpd( ar, br ), vmulpd( ai, bi ) ) );
    }        
    vmovapd( dest + i*vstep64, res );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}

static void hlM3C64RealTraceMulTranspose( r64 *dest, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 res = vzeropd();
    for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->real.vals ); ++matrix_index ){
      v64 ar = vmovapd( a->real.vals[ matrix_index ] + i*vstep64 );
      v64 br = vmovapd( b->real.vals[ matrix_index ] + i*vstep64 );
      v64 ai = vmovapd( a->imag.vals[ matrix_index ] + i*vstep64 );
      v64 bi = vmovapd( b->imag.vals[ matrix_index ] + i*vstep64 );
      res = vaddpd( res, vsubpd( vmulpd( ar, br ), vmulpd( ai, bi ) ) );
    }        
    vmovapd( dest + i*vstep64, res );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}

static void hlM3C64ImaginaryTraceMulAdj( r64 *dest, m3c64 *a, m3c64 *b ){
#if defined( hlSIMD )
  for( u64 i = 0; i < a->real.count / vstep64; ++i ){
    v64 res = vzeropd();
    for( u64 matrix_index = 0; matrix_index < hlARRAYCOUNT( a->real.vals ); ++matrix_index ){
      v64 ar = vmovapd( a->real.vals[ matrix_index ] + i*vstep64 );
      v64 br = vmovapd( b->real.vals[ matrix_index ] + i*vstep64 );
      v64 ai = vmovapd( a->imag.vals[ matrix_index ] + i*vstep64 );
      v64 bi = vmovapd( b->imag.vals[ matrix_index ] + i*vstep64 );
      res =  vsubpd( res, vaddpd( vmulpd( ar, bi ), vmulpd( ai, br ) ) );
    }        
    vmovapd( dest + i*vstep64, res );
  }
#else
  hlASSERT( !"NOT YET IMPLEMENTED FOR NON SIMD\n" );
#endif       
}


//-------------------------------------------------------------------------------------------//
//--                                   STRING LIBRARY                                      --//
//-------------------------------------------------------------------------------------------//

static u64 hlStrLengthC( c8 *str ){
  u64 result = 0;
  while( *str++ ){
    ++result; 
  }
  return( result );
}

static b32 hlIsWhitespace( c8 c ){
  b32 result = ( c == ' ' || c == '\t' || c == '\n' );
  return( result );
}

static b32 hlStrStartsWithC( c8 *prefix, c8 *str ){
  b32 result = 1;
  while( result && *prefix ){
    result = *prefix++ == *str++;
  }
  return( result );
}

static b32 hlStrContainsC( c8 *prefix, c8 *str ){
  b32 result = 0;
  u64 prefix_len = hlStrLengthC( prefix );
  while( *str && !result ){
   for( u64 i = 0; i < prefix_len && str[i]; ++i ){
     result = str[i] == prefix[i];
     if( !result ){
       i = prefix_len;
     }
   }
    ++str;
  }

  return( result );
}

static b32 hlStrEndsWithC( c8 *prefix, c8 *str ){
  b32 result = 0;
  u64 prefix_len = hlStrLengthC( prefix );
  u64 str_len    = hlStrLengthC( str    );
  if( str_len >= prefix_len ){
    result = 1;
    for( u64 str_index = 0; str_index < prefix_len; ++str_index ){
      result &= prefix[str_index] == str[str_len-prefix_len+str_index];
    }
  }

  return( result );
}

static b32 hlStrCmpC( c8 *str_a, c8 *str_b ){
  b32 result = 1;
  while( result && (*str_a || *str_b) ){
    result = *str_a++ == *str_b++;
  }
  return( result );
}

static u64 hlStrCountTokkenC( c8 *str, c8 tokken ){
  u64 result = 0;
  while( *str ){
    if( *str++ == tokken ){
      ++result;
    }
  }
  return( result );
}


static u64 hlCountCharacterC( c8 *str, c8 tokken ){
  u64 result = 0;
  for( c8 *s = str; *s; ++s ){
    if( *s == tokken ){
      ++result;
    }
  }
  return( result );
}

static b32 hlIsNumericC( c8 *str ){
  hlASSERT( str );
  b32 result = 1;
  for( c8 *s = str; *s; ++s ){
    char c = *s;
    if( !( ( c >= '0' && c <= '9' ) || c == 'e' || c == '-' || c == '.' || c == ' ' || c == 'n' || c == 'a' || c == '+' ) ){
      result = 0;
      break;
    }
  }
  return( result );
}

static u64 hlCountWhitespaceC( c8 *str ){
  u64 result = 0;
  for( c8 *s = str; *s; ++s ){
    if( hlIsWhitespace( *s ) ){
      ++result;
    }
  }
  return( result );
}

static c8 *hlGobbleWhitespaceC( c8 *str ){
  c8 *result = str;
  while( hlIsWhitespace( *result ) ){
    ++result;
  }
  return( result );
}

static c8 *hlStrGuardedNextSpaceC( c8 *str, c8 guard ){
  c8 *result = str;
  result = hlGobbleWhitespaceC( result );
  while( !hlIsWhitespace( *result ) && *result != guard ){
    ++result;
  }
  return( result );
}

static c8 *hlStrNextSpaceC( c8 *str ){
  c8 *result = str;
  result = hlGobbleWhitespaceC( result );
  while( !hlIsWhitespace( *result ) ){
    ++result;
  }
  return( result );
}

static c8 *hlStrGobbleTilC( c8 *til, c8 *str ){
  c8 *result = str;
  while( !hlStrStartsWithC( til, result ) ){
    ++result;
  }
  result += hlStrLengthC( til );
  return( result );
}

static c8 *hlStrGrabTilC( hlMem_Pool *pool, c8 *str, c8 *til, c8 **output, u64 *output_length ){
  c8 *result = hlStrGobbleTilC( til, str );
  u64 grab_length = (u64)(result - str) - 1;
  *output = (c8 *)hlPUSHSIZE( pool, grab_length + 1 );
  *output_length = grab_length;
  hlMemCopy( str, *output, grab_length );
  (*output)[ grab_length ] = '\0';
  return( result );
}

static u64 hlStrFindFirstOccurenceC( c8 c, c8 *str ){
  u64 result = 0;
  for( c8 *ptr = str; *ptr; ++ptr ){
    if( *ptr == c ){
      result = (u64)(ptr - str);
      break;
    }
  }
  return( result );
}

static u64 hlStrFindNthOccurenceC( c8 c, u64 n, c8 *str ){
  u64 result = 0;
  for( c8 *ptr = str; n && *ptr; ++ptr ){
    if( *ptr == c ){
      result = (u64)(ptr - str);
      --n;
    }
  }
  return( result );
}

static u64 hlStrFindNthToLastOccurenceC( c8 c, u64 n, c8 *str ){
  u64 counter = 0;
  for( c8 *ptr = str; *ptr; ++ptr ){
    if( *ptr == c ){
      ++counter;
    }
  }
  u64 result = hlStrFindNthOccurenceC( c, counter - n + 1, str );
  return( result );
}

static u64 hlStrFindLastOccurenceC( c8 c, c8 *str ){
  u64 result = 0;
  for( c8 *ptr = str; *ptr; ++ptr ){
    if( *ptr == c ){
      result = (u64)(ptr - str);
    }
  }
  return( result );
}

static c8 *hlStrAfterLastC( c8 c, c8 *str ){
  c8 *result = str;
  while( *str ){
    if( *str == c ){
      result = str+1;
    }
    ++str;
  }
  return( result );
}
 
static c8 *hlStrAfterC( c8 c, c8 *str ){
  c8 *result = str;
  while( *str ){
    if( *str == c ){
      result = str+1;
      break;
    }
    ++str;
  }
  return( result );
}

static c8 *hlStrCopy( hlMem_Pool *pool, u64 len, c8 *str ){
  c8 *result;
  result = hlPUSHARRAY( pool, len+1, c8 );
  result[len] = '\0';
  hlMemCopy( str, result, len );
  return( result );
}

static c8 *hlStrCopyC( hlMem_Pool *pool, c8 *str ){
  c8 *result;
  u64 str_length = hlStrLengthC( str ) + 1;
  result = hlPUSHARRAY( pool, str_length, c8 );
  hlMemCopy( str, result, str_length );
  return( result );
}

static c8 *hlStrCatC( hlMem_Pool *pool, c8 *prefix, c8 *postfix ){
  u64 prefix_len  = hlStrLengthC( prefix  );
  u64 postfix_len = hlStrLengthC( postfix );
  c8 *result = hlPUSHARRAY( pool, prefix_len + postfix_len + 1, c8 );
  hlMemCopy( prefix , result           , prefix_len  );
  hlMemCopy( postfix, result+prefix_len, postfix_len );
  result[ prefix_len + postfix_len ] = '\0';
  return( result );
}

typedef struct {
  u64 len;
  c8  val[0];
} hlStr;

static c8 **hlStrSplitC( hlMem_Pool *pool, c8 tok, u64 *res_len, c8 *str ){
  u64 tokken_count = hlStrCountTokkenC( str, tok );
  c8 **result = hlPUSHARRAY( pool, tokken_count, c8 * );
  c8 grab_str[2] = { tok, '\0' };
  u64 dummy_output = 0;
  for( u64 tokken_index = 0; tokken_index < tokken_count; ++tokken_index ){
//static c8 *hlStrGrabTilC( hlMem_Pool *pool, c8 *str, c8 *til, c8 **output, u64 *output_length ){
    str = hlStrGrabTilC( pool, str, grab_str, result + tokken_index, &dummy_output );
  }
  *res_len = tokken_count;
  return( result );
}

//-------------------------------------------------------------------------------------------//
//--                            PORTABLE PLATFORM LIBRARY                                  --//
//--                            BASED OF HANDMADE HERO                                     --//
//-------------------------------------------------------------------------------------------//


#if defined( __GNUC__ )
#include <sys/mman.h>
#endif

u8 *hlVirtualAlloc( u64 base, u64 size ){
#if defined( _MSC_VER )
  u8 *result = (u8 *)VirtualAlloc( (void *)base, size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE );
  hlASSERT( result );
  return( result );
#elif defined( __GNUC__ )
  u8 *result = (u8 *)mmap( (void *)base, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_FIXED | MAP_ANON, -1, 0 );
  hlASSERT( result != (void *)-1 );
  return( result );
#endif
}

b32 hlVirtualFree( u8 *base, u64 size ){
#if defined( _MSC_VER )
  b32 result = VirtualFree( (void *)base, 0, MEM_RELEASE );
  hlASSERT( result );
  return( result );
#elif defined( __GNUC__ )
  b32 result = munmap( (void *)base, size );
  hlASSERT( result >= 0 );
  return( result );
#endif
}

static void hlInitializeProgramMemory( hlMem_Pool *pool1, u64 pool1_size, hlMem_Pool *pool2, u64 pool2_size ){
  hlDEBUG_STR( "VIRTUALLY ALLOCATING MEMORY\n" );
  u8 *memory = hlVirtualAlloc( GIGABYTES( 5 ), pool1_size + pool2_size );
  hlASSERT( memory );
  hlInitializeMemoryPool( pool1, memory             , pool1_size );
  hlInitializeMemoryPool( pool2, memory + pool1_size, pool2_size );
}

static c8 *sizes[] = {
  ""
, "KB"
, "MB"
, "GB"
, "PB"
};

static r64 hlHumanReadableSize( u64 size, c8 **str ){
  u64 power = 0;
  r64 res = (r64)size;
  while( res > 1024.0 ){
    res = res / 1024;
    ++power;
  }
  *str = sizes[power];
  return res;
}

static void hlPrintProgramMemory( hlMem_Pool *tmp, hlMem_Pool *mem ){
  hlDEBUG_STR( "Program Memory Diagnostics\n" );
  c8 *used_str;
  c8 *size_str;
  r64 used_size = hlHumanReadableSize( tmp->used, &used_str );
  r64 max_size  = hlHumanReadableSize( tmp->size, &size_str );
  printf( "TMP MEM: %f%s / %f%s\n", used_size, used_str, max_size, size_str );
  used_size = hlHumanReadableSize( mem->used, &used_str );
  max_size  = hlHumanReadableSize( mem->size, &size_str );
  printf( "MAIN MEM: %f%s / %f%s\n", used_size, used_str, max_size, size_str );
}

typedef struct {
  u64 size;
  u8 *contents;
} hlFile;

#include <stdio.h>//--TODO--//
static hlFile hlFileRead( hlMem_Pool *mem, c8 *path ){
  hlFile result;
#if defined( _MSC_VER )
#elif defined( __GNUC__ )
  FILE *file = fopen( path, "r" );
  fseek( file, 0, SEEK_END );
  result.size = ftell( file );
  fseek( file, 0, SEEK_SET );
  result.contents = hlPUSHSIZE( mem, result.size );
  fread( result.contents, result.size, 1, file );
  fclose( file );
#endif
  return( result );
}

#include <stdio.h>//--TODO--//
static void hlFileWrite( void *data, u64 size, c8 *path ){
#if defined( _MSC_VER )

#elif defined( __GNUC__ )
  FILE *file = fopen( path, "w" );
  fwrite( data, size, 1, file );
  fclose( file );
#endif
}

//static HARMLESS_PLATFORM_COPY_ENTIRE_FILE( Clang64CopyEntireFile ){
//  u64   start_used = pool->used;
//  void *src_file   = Clang64ReadEntireFile( src, pool );
//  Clang64WriteEntireFile( dest, src_file, pool->used - start_used );
//  pool->used = start_used;
//}

#define HARMLESS_PLATFORM_TSC_TO_INTERVAL( name ) r32 name( u64 start, u64 end )
typedef HARMLESS_PLATFORM_TSC_TO_INTERVAL( hlPlatformhlTSCToInterval );

typedef struct {
  hlMem_Pool *permanent_memory;
  hlMem_Pool *transient_memory;
} hlPlatform_Memory;

typedef struct {
  hlPlatformhlTSCToInterval *hlTSCToInterval;
} hlPlatform_API;

#define HARMLESS_MAIN( name ) void name( hlPlatform_API    *api    \
                                       , hlPlatform_Memory *memory \
                                       , i64                argc   \
                                       , c8               **argv   )
typedef HARMLESS_MAIN( Main );

#define HARMLESS_UPDATE( name ) b32 name( hlPlatform_API    *api    \
                                        , hlPlatform_Memory *memory \
                                        , r32                dt     )
typedef HARMLESS_UPDATE( Update );

#if defined( PLATFORM_EXE_ENTRYPOINT ) || defined( PLATFORM_EXE_UPDATE )
#if !defined( PLATFORM_EXE_LIB )
  #error MUST DEFINE <PLATFORM_EXE_ENTRYPOINT> OR <PLATFORM_EXE_UPDATE> IF USING <PLATFORM_EXE_ENTRYPOINT>
#endif

static u64 global_tsc_constant_frequency_ms;
static HARMLESS_PLATFORM_TSC_TO_INTERVAL( hlTSCToInterval ){
  hlASSERT( global_tsc_constant_frequency_ms ); //-- RELIANT ON ZEROED DEFAULT MEMORY --//
  r32 result = (r32)(end - start) / (r32)global_tsc_constant_frequency_ms;
  //r32 result = (r32)( 1000 * (end - start) / global_tsc_constant_frequency_ms ) / 1000.0f;
  return( result );
}

#if defined( _MSC_VER )
#include <gl/gl.h>

void WinMainCRTStartup(){
  int result = WinMain( GetModuleHandle( 0 ), 0, 0, 0 );
  ExitProcess( result );
}

//static void ToggleFullscreen( HWND window ){
//  //-- COPY PASTA FROM RAYMOND CHEN --//
//  //-- <http://blogs.msdn.com/b/oldnewthing/archive/2010/04/12/9994016.aspx> --//
//  DWORD style = GetWindowLong( window, GWL_STYLE );
//  if( style & WS_OVERLAPPEDWINDOW ){
//    MONITORINFO monitor_info = { sizeof( monitor_info ) };
//    if( GetWindowPlacement( window, &global_window_prev_position )
//    &&  GetMonitorInfo( MonitorFromWindow( window, MONITOR_DEFAULTTOPRIMARY )
//                      , &monitor_info ) ){
//      SetWindowLong( window, GWL_STYLE, style & ~WS_OVERLAPPEDWINDOW );
//      SetWindowPos( window, HWND_TOP
//                  , monitor_info.rcMonitor.left
//                  , monitor_info.rcMonitor.top
//                  , monitor_info.rcMonitor.right  - monitor_info.rcMonitor.left
//                  , monitor_info.rcMonitor.bottom - monitor_info.rcMonitor.top
//                  , SWP_NOOWNERZORDER | SWP_FRAMECHANGED );
//    }
//  }else{
//    SetWindowLong( window, GWL_STYLE, style | WS_OVERLAPPEDWINDOW );
//    SetWindowPlacement( window, &global_window_prev_position );
//    SetWindowPos( window, 0, 0, 0, 0, 0
//                , SWP_NOMOVE | SWP_NOSIZE | SWP_NOZORDER | SWP_NOOWNERZORDER | SWP_FRAMECHANGED );
//  }
//}

static void Win64PushBuffer( HWND window ){
    PAINTSTRUCT paint;
    HDC device_context = BeginPaint( window, &paint );

    RECT client_rect;
    GetClientRect( window, &client_rect );
    u32 width  = client_rect.right - client_rect.left;
    u32 height = client_rect.bottom - client_rect.top;

    glViewport( 0, 0, width, height );
    
    glClearColor( 1.0f, 0.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT );

    glMatrixMode( GL_TEXTURE );
    glLoadIdentity();
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    
    glBegin( GL_TRIANGLES );
    glColor3f( 1.0f, 0.0f, 0.0f );
    glVertex2f( -1.0f, -1.0f );
    glVertex2f(  1.0f, -1.0f );
    glVertex2f(  1.0f,  1.0f );

    glColor3f( 0.0f, 1.0f, 0.0f );
    glVertex2f( -1.0f, -1.0f );
    glVertex2f(  1.0f,  1.0f );
    glVertex2f( -1.0f,  1.0f );
    glEnd();
    
    SwapBuffers( device_context );
    EndPaint( window, &paint );
}

LRESULT CALLBACK Win64MainWindowCallback( HWND window, UINT message, WPARAM w_param, LPARAM l_param ){
  LRESULT Result = 0;

  switch( message ){
  case WM_SIZE: { 
  }break;

  case WM_SETCURSOR: {
    SetCursor( 0 );
  }break;

  case WM_CLOSE: {
    DEBUG_STR( "CLOSING WINDOW CALLBACK\n" );
  }break;

  case WM_DESTROY: {
    DEBUG_STR( "DESTRYING WINDOW CALLBACK\n" );
  }break;

  case WM_ACTIVATEAPP: {
    DEBUG_STR( "APP FOCUSED\n" );
  }break;

  case WM_PAINT: {
    Win64PushBuffer( window );
    OutputDebugStringA( "CALLED PAINT\n" );
  }break;

  default: {
    Result = DefWindowProcA( window, message, w_param, l_param );
  }break;

  }
  return( Result );
;
}

static void Win64InitOpenGL( HWND window ){
  HDC window_dc = GetDC( window );

  PIXELFORMATDESCRIPTOR desired_pixel_format            = {};
                        desired_pixel_format.nSize      = sizeof( desired_pixel_format );
                        desired_pixel_format.nVersion   = 1;
                        desired_pixel_format.iPixelType = PFD_TYPE_RGBA;
                        desired_pixel_format.dwFlags    = PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW 
                                                        | PFD_DOUBLEBUFFER;
                        desired_pixel_format.cColorBits = 32; //-- "MEANT" TO EXCLUDE ALPHA --//
                        desired_pixel_format.cAlphaBits = 8;
                        desired_pixel_format.iLayerType = PFD_MAIN_PLANE;

  i32 suggested_pixel_format_index = ChoosePixelFormat( window_dc, &desired_pixel_format );
  PIXELFORMATDESCRIPTOR suggested_pixel_format;
  DescribePixelFormat( window_dc, suggested_pixel_format_index
                     , sizeof( suggested_pixel_format ), &suggested_pixel_format );
  SetPixelFormat( window_dc, suggested_pixel_format_index, &suggested_pixel_format );

  HGLRC rendering_context = wglCreateContext( window_dc );
  if( wglMakeCurrent( window_dc, rendering_context ) ){
    
  }else{
    hlASSERT( !"WAS NOT ABLE TO MAKE RC CURRENT" );
  }
  ReleaseDC( window, window_dc );
}

int CALLBACK WinMain( HINSTANCE instance, HINSTANCE prev_instance, LPSTR lp_cmd_line, i32 cmd_show ){
  Platform_API platform_api = {};
                 platform_api.hlTSCToInterval = &hlTSCToInterval;

  Platform_Memory platform_memory = {};

  //-- TODO FOR DYNAMIC RELOAD --//
  //-- GET DLL                 --//
  //-- MAKE TEMP DLL PATH      --//
  //-- LOCK FILE               --//

  HMODULE code_dll = LoadLibraryA( STR( PLATFORM_EXE_LIB ) );
  hlASSERT( code_dll );
#if defined( PLATFORM_EXE_ENTRYPOINT )
  Main *application_entrypoint = (Main *)GetProcAddress( code_dll, STR( PLATFORM_EXE_ENTRYPOINT ) );
  hlASSERT( application_entrypoint );
#endif
#if defined( PLATFORM_EXE_UPDATE )
  Update *update_function = (Update *)GetProcAddress( code_dll, STR( PLATFORM_EXE_UPDATE ) );
  hlASSERT( update_function );
#endif

#if defined( PLATFORM_EXE_ENTRYPOINT )
  application_entrypoint( &platform_api, &platform_memory, 1, lp_cmd_line );
#endif

#if defined( PLATFORM_EXE_UPDATE )
  LARGE_INTEGER frequency;
  QueryPerformanceFrequency( &frequency );
  global_tsc_constant_frequency_ms = frequency.QuadPart;

  u32 desired_scheduler_ms = 1;
  MMRESULT scheduler_return_value = timeBeginPeriod( desired_scheduler_ms );
  hlASSERT( scheduler_return_value == TIMERR_NOERROR );

  WNDCLASS window_class               = {};
           window_class.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC; 
           window_class.lpfnWndProc   = Win64MainWindowCallback;
           window_class.hInstance     = instance;
           window_class.hCursor       = LoadCursor( 0, IDC_ARROW );
           //window_class.hIcon = ;
           window_class.lpszClassName = "Test Window Class";

  if( RegisterClassA( &window_class ) ){
    HWND window = CreateWindowExA( 0
                                 , window_class.lpszClassName
                                 , "Test Window"
                                 , WS_OVERLAPPEDWINDOW | WS_VISIBLE
                                 , CW_USEDEFAULT
                                 , CW_USEDEFAULT
                                 , CW_USEDEFAULT
                                 , CW_USEDEFAULT
                                 , 0
                                 , 0
                                 , instance
                                 , 0 );

    if( window ){
      HDC window_context = GetDC( window );
      Win64InitOpenGL( window );

      i32 monitor_refreshrate = GetDeviceCaps( window_context, VREFRESH );

      u32 target_refreshrate   = monitor_refreshrate > 1 ? monitor_refreshrate : 60;
      r32 application_hz       = (r32)target_refreshrate;
      r32 target_frame_time_ms = 1000.0f / application_hz;
      
      b32 running         = true;
      u64 loop_start_tsc  = __rdtsc();
      r32 time_elapsed_ms = 0.0f;

      while( running ){
        MSG message;
        while( PeekMessageA( &message, 0, 0, 0, PM_REMOVE ) ){
          if( message.message == WM_QUIT ){
            DEBUG_STR( "APPLICATION QUIT CALL IN THE OTHER PLACE" );
          }else{
            //-- YES REALLY --//
            //-- WE DEFINED A MESSAGE CALLBACK AND THEN STIL HAVE TO DO THIS --//
            TranslateMessage( &message );
            DispatchMessageA( &message );
          }
        }

        running = update_function( &platform_api, &platform_memory, time_elapsed_ms );
        Win64PushBuffer( window );
        
        u64 loop_end_tsc = __rdtsc();
        time_elapsed_ms  = hlTSCToInterval( loop_start_tsc, loop_end_tsc );

        //char print_buffer[ 256 ];
        //_snprintf_s( print_buffer, sizeof( print_buffer ), "%.02f ms/frame\n", time_elapsed_ms );

        if( time_elapsed_ms < target_frame_time_ms ){
          i32 sleep_ms = (u32)( target_frame_time_ms - time_elapsed_ms ) - 1; //TODO BETTER ROUND
          if( sleep_ms > 0 ){
            Sleep( sleep_ms );
            time_elapsed_ms = hlTSCToInterval( loop_start_tsc, __rdtsc() );
            //hlASSERT( time_elapsed_ms < target_frame_time_ms );
            if( time_elapsed_ms < target_frame_time_ms ){
              DEBUG_STR( "SLEEP CALL IS UNRELIABLE\n" ); //--TODO FIX THIS--//
            }
          }
          while( time_elapsed_ms < target_frame_time_ms ){
            time_elapsed_ms = hlTSCToInterval( loop_start_tsc, __rdtsc() );
          }
        }else{
          DEBUG_STR( "LONG OUTPUT TIME DETECTED\n" );
        }

        loop_start_tsc = loop_end_tsc;
      }

      ReleaseDC( window, window_context );
    }else{
      hlASSERT( !"Failed to create window" );
    }

  }else{
    hlASSERT( !"Failed to register window" );
  }
#endif
  return( 0 );
}

#elif defined( __GNUC__ )

#include <X11/Xlib.h>
#include <GL/glx.h>
#include <stdio.h>//--TODO--//
#include <sys/stat.h>
i32 usleep(u32 usec);
int main( int argc, char **argv ){
  Platform_API platform_api                 = {};
               platform_api.hlTSCToInterval = &hlTSCToInterval;

  Platform_Memory platform_memory = {};

  hlMem_Pool lib_pool;
  u64        lib_pool_size = MEGABYTES( 5 );
  void      *lib_pool_mem  = hlVirtualAlloc( GIGABYTES( 6 ), lib_pool_size );
  hlInitializeMemoryPool( &lib_pool, (u8 *)lib_pool_mem, lib_pool_size );

  //--TODO PLACE IN SOME REASONABLE DIRECTORY --//
  char *lib_path              = STR( PLATFORM_EXE_LIB_PATH ) STR( PLATFORM_EXE_LIB );
  char *old_lib_path          = STR( PLATFORM_EXE_LIB_PATH ) "tmp_old_" STR( PLATFORM_EXE_LIB );
  char *current_lib_path      = STR( PLATFORM_EXE_LIB_PATH ) "tmp_new_" STR( PLATFORM_EXE_LIB );
  char *lib_compile_lock_path = STR( PLATFORM_EXE_LIB_PATH ) "lock.tmp";

  Clang64CopyEntireFile( lib_path, current_lib_path, &lib_pool );
  Clang64CopyEntireFile( lib_path, old_lib_path    , &lib_pool );

  void *code_lib = dlopen( current_lib_path, RTLD_NOW );
  hlASSERT( code_lib );
#if defined( PLATFORM_EXE_ENTRYPOINT )
  Main *application_entrypoint = (Main *)dlsym( code_lib, STR( PLATFORM_EXE_ENTRYPOINT ) );
  hlASSERT( application_entrypoint );
#endif
#if defined( PLATFORM_EXE_UPDATE )
  Update *update_function = (Update *)dlsym( code_lib, STR( PLATFORM_EXE_UPDATE ) );
  hlASSERT( update_function );
#endif

#if defined( PLATFORM_EXE_ENTRYPOINT )
  application_entrypoint( &platform_api, &platform_memory, argc, argv );
#endif

#if defined( PLATFORM_EXE_UPDATE )
  //--TODO ASSUMING CONSTANT RDTSC FREQUENCY --//
  u64 start_tsc = __rdtsc();
  usleep(1000); 
  global_tsc_constant_frequency_ms = (__rdtsc() - start_tsc);
  printf( "SET TSC FREQUENCY TO: %llu\n", global_tsc_constant_frequency_ms );

  Display *x11_display = XOpenDisplay( 0 );
  if( x11_display ){
    Window x11_root = DefaultRootWindow( x11_display );;
    GLint gl_attributes[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    XVisualInfo *x11_info = glXChooseVisual( x11_display, 0, gl_attributes );
    if( x11_info ){
      Colormap colour_map = XCreateColormap( x11_display, x11_root, x11_info->visual, AllocNone );
      XSetWindowAttributes x11_window_attributes = {};
                           x11_window_attributes.colormap   = colour_map;
                           x11_window_attributes.event_mask = ExposureMask | KeyPressMask;

      Window x11_window = XCreateWindow( x11_display, x11_root
                                       , 0, 0, 800, 600
                                       , 0, x11_info->depth, InputOutput
                                       , x11_info->visual, CWColormap | CWEventMask
                                       , &x11_window_attributes );
      XMapWindow( x11_display, x11_window ); //-- DISPLAY WINDOW --//
      XStoreName( x11_display, x11_window, "TestWindow" );

      GLXContext glx_context = glXCreateContext( x11_display, x11_info, 0, GL_TRUE );
      glXMakeCurrent( x11_display, x11_window, glx_context );
      glEnable( GL_DEPTH_TEST );

      b32 running              = true;
      u64 loop_start_tsc       = __rdtsc();
      r32 time_elapsed_ms      = 0.0f;
      r32 last_time_elapsed_ms = 0.0f;
      while( running ){

#if defined( hlDEBUG )
        //--------------------------------//
        //-- PERFORM HOT CODE RELOADING --//
        //--------------------------------//
        struct stat lock_file_stat;
        struct stat current_lib_stat;
        struct stat new_lib_stat;
        stat( current_lib_path     , &current_lib_stat );
        stat( lib_path             , &new_lib_stat     );

        struct timespec new_ts     = new_lib_stat.st_mtimespec;
        struct timespec current_ts = current_lib_stat.st_mtimespec;
        
        b32 code_compiling = !( stat( lib_compile_lock_path, &lock_file_stat ) >= 0 );

        if( !code_compiling
         || (new_ts.tv_sec > current_ts.tv_sec) 
         || ((new_ts.tv_sec == current_ts.tv_sec) && (new_ts.tv_nsec > current_ts.tv_nsec)) ){
          i32 closed = dlclose( code_lib ); //--TODO TECHNICALLY NEED TO LOOK AT application_entrypoint AS WELL--//
          hlASSERT( closed >= 0 );
          Clang64CopyEntireFile( current_lib_path, old_lib_path    , &lib_pool );
          Clang64CopyEntireFile( lib_path        , current_lib_path, &lib_pool );
          
          code_lib = dlopen( current_lib_path, RTLD_NOW );
          hlASSERT( code_lib );
          update_function = (Update *)dlsym( code_lib, STR( PLATFORM_EXE_UPDATE ) );
          hlASSERT( update_function );
        }
#endif

        //--------------------------//
        //-- HANDLE WINDOW EVENTS --//
        //--------------------------//
        //u64 x11_start_tsc = __rdtsc();
        XEvent event;
        while( XCheckWindowEvent( x11_display, x11_window, ExposureMask | KeyPressMask, &event ) ){
          if( event.type == Expose ){ //-- REDRAW DISPLAY POSSIBLE RESIZE --//
            hlDEBUG_STR( "APPLICATION RECEIVED EXPOSE COMAND\n" );
            XWindowAttributes window_attributes;
            XGetWindowAttributes( x11_display, x11_window, &window_attributes );
            glViewport( 0, 0, window_attributes.width, window_attributes.height );
            glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
            glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
            glXSwapBuffers( x11_display, x11_window );
          }
        }
        //r32 x11_time_ms = hlTSCToInterval( x11_start_tsc, __rdtsc() );
        //printf( "X11 took: %f\n", x11_time_ms );

        running = update_function( &platform_api, &platform_memory, last_time_elapsed_ms );
        
        u64 loop_end_tsc = __rdtsc();
        time_elapsed_ms  = hlTSCToInterval( loop_start_tsc, loop_end_tsc );

        r32 target_frame_time_ms = 1000.0f / 60.0f; //-- TODO --//
        //printf( "%f\n", time_elapsed_ms );
        if( time_elapsed_ms < target_frame_time_ms ){
          i32 sleep_ms = (u32)(target_frame_time_ms - time_elapsed_ms) - 1;
          if( sleep_ms > 0 ){
            usleep( 1000 * sleep_ms );
            time_elapsed_ms = hlTSCToInterval( loop_start_tsc, __rdtsc() );
            //hlASSERT( time_elapsed_ms < target_frame_time_ms );//--TODO SERIOUS LOG--//
          }

          while( time_elapsed_ms < target_frame_time_ms ){
            time_elapsed_ms = hlTSCToInterval( loop_start_tsc, __rdtsc() );
          }
        }else{
          hlDEBUG_STR( "LONG OUTPUT TIME DETECTED\n" );
        }

        loop_start_tsc       = __rdtsc();
        last_time_elapsed_ms = time_elapsed_ms;
      }
    }else{
      hlDEBUG_STR( "FAILED TO CREATE GLX VISUAL INFO\n" );
    }

  }else{
    hlDEBUG_STR( "FAILED TO CREATE X11 DISPLAY\n" );
  }
#endif

  return( 0 );
}
#else
#error UNSUPPORTED PLATFROM
#endif
#endif //-- PLATFORM_EXE_ENTRYPOINT --// 

#if defined( hlDEBUG )
#include <stdio.h>
#endif

//-- XML TYPES --//
typedef enum hl_xml_type {
  hlXML_IGNORED 
, hlXML_ARR     
, hlXML_OBJECT  
, hlXML_NUM_ARR 
, hlXML_FLT_ARR 
, hlXML_NUM     
, hlXML_FLT     
, hlXML_STR
} hl_xml_type;

c8 *hl_xml_type_to_str[] = {
  "hlXML_IGNORED"
, "hlXML_ARR"     
, "hlXML_OBJECT"  
, "hlXML_NUM_ARR" 
, "hlXML_FLT_ARR" 
, "hlXML_NUM"     
, "hlXML_FLT"     
, "hlXML_STR"
};

typedef struct hlXML_Obj {
  u64           count;
#if defined( hlDEBUG )
  u64           size;
#endif
  struct hlXML *children;
} hlXML_Obj;

typedef struct hlXML_Num_Arr {
  u64  count;
  i64 *values;
} hlXML_Num_Arr;

typedef struct hlXML_Flt_Arr {
  u64  count;
  r64 *values;
} hlXML_Flt_Arr;

typedef struct hlXML {
  struct hlXML *parent;
  hl_xml_type   type; 
  char          name[64];
  union {
    struct hlXML_Obj     obj_val;
    struct hlXML_Num_Arr num_arr_val;
    struct hlXML_Flt_Arr flt_arr_val;
    i64                  num_val;
    r64                  flt_val;
    c8                  *str;
  };
} hlXML;

//-- XML PROTO TYPES --//
enum hlXML_Tokkens {
  hlXML_TOKKEN_IGNORED
, hlXML_TOKKEN_START
, hlXML_TOKKEN_END
, hlXML_TOKKEN_ELEM_START
, hlXML_TOKKEN_ELEM_END
, hlXML_TOKKEN_VAL
};

typedef struct hlTokken_XML {
  struct hlTokken_XML *previous;
  struct hlTokken_XML *next;
  enum hlXML_Tokkens   type;
  char                *value;
  u64                  value_length;
} hlTokken_XML;

static u64 hlXMLCountChildren( hlTokken_XML *parent ){
  hlASSERT( parent->next );
  hlASSERT( parent->next->type == hlXML_TOKKEN_START || parent->next->type == hlXML_TOKKEN_ELEM_START );
  u64 result = 0;
  u64 level = 1;
  hlTokken_XML *current_tokken = parent;
  while( level ){
    current_tokken = current_tokken->next;
    if( current_tokken->type == hlXML_TOKKEN_START || current_tokken->type == hlXML_TOKKEN_ELEM_START ){
      if( level == 1 ){
        ++result;
      }
      ++level;
    }else if( current_tokken->type == hlXML_TOKKEN_END || current_tokken->type == hlXML_TOKKEN_ELEM_END ){
      --level;
    }
  }
  return( result );
}

static hlTokken_XML *hlGenerateXMLTokken( hlMem_Pool *pool, hlTokken_XML *current ){
  hlTokken_XML *result = current;
  if( result->type ){
    result = hlPUSHTYPE( pool, hlTokken_XML );
    result->type     = hlXML_TOKKEN_IGNORED;
    current->next    = result;
    result->previous = current;
  }
  return( result );
}

static hlXML *hlXMLAddChild( hlXML *parent ){
  hlASSERT( parent->type == hlXML_OBJECT || parent->type == hlXML_ARR );
  //printf( "%lu / %lu\n", parent->obj_val.count, parent->obj_val.size );
  hlASSERT( parent->obj_val.count < parent->obj_val.size );

  hlXML *result;
         result         = parent->obj_val.children + parent->obj_val.count++;
         result->parent = parent;
  if( parent->type == hlXML_ARR ){
    sprintf( result->name, "%llu", parent->obj_val.count );
  }

  return( result );
}

static hlXML hlParseXML( hlMem_Pool *tmp_pool, hlMem_Pool *pool, c8 *str, u64 str_length ){
  //printf( "PARSING: \n%s\n", str );
  hlXML result = { 0 };
  u64 tmp_pool_starting_size = tmp_pool->used;

  //-- FILL THE TEMP POOL WITH THE TOKKENISED VALUES --//
  hlTokken_XML head_tokken     = { 0 };
  hlTokken_XML *current_tokken = &head_tokken;
  char *position = str;
  position = hlGobbleWhitespaceC( position );
  while( position < str + str_length ){
    //printf( "(%lu,%lu)\n", position, str + str_length );
    //write( 1, "TAG", 3 );
    //write( 1, position, 60 );
    //write( 1, "\n", 1 );
    hlTokken_XML *previous_tokken = current_tokken;
    current_tokken = hlGenerateXMLTokken( tmp_pool, current_tokken );

    if( hlStrStartsWithC( "<?", position ) ){
      position = hlStrGobbleTilC( "?>", position );
    }else if( hlStrStartsWithC( "<!--", position ) ){
      position = hlStrGobbleTilC( "-->", position );
    }else if( hlStrStartsWithC( "<elem>", position ) ){
      position += 6;
      current_tokken->type = hlXML_TOKKEN_ELEM_START;
      //printf( "ES " );
    }else if( hlStrStartsWithC( "</elem>", position ) ){
      position += 7;
      current_tokken->type = hlXML_TOKKEN_ELEM_END;
      //printf( "EE " );
    }else if( hlStrStartsWithC( "</", position ) ){
      position = hlStrGrabTilC( tmp_pool, position + 2, ">", &current_tokken->value, &current_tokken->value_length );
      current_tokken->type = hlXML_TOKKEN_END;
      //printf( "END " );
    }else if( *position == '<' ){
      position = hlStrGrabTilC( tmp_pool, position + 1, ">", &current_tokken->value, &current_tokken->value_length );
      current_tokken->type = hlXML_TOKKEN_START;
      //printf( "START " );
    }else if( !hlIsWhitespace( *position ) && ( previous_tokken->type == hlXML_TOKKEN_START || previous_tokken->type == hlXML_TOKKEN_ELEM_START ) ){
      position = hlStrGrabTilC( tmp_pool, position, "<", &current_tokken->value, &current_tokken->value_length );
      position--;
      current_tokken->type = hlXML_TOKKEN_VAL;
      //printf( "val " );
    }else{
      write( 1, position, 100 );
      hlASSERT( !"XML PARSE FAILED" );
    }
    //write( 1, position, 100 );
    //write( 1, "\n", 1 );
    position = hlGobbleWhitespaceC( position );
    //printf( "%lu - %lu\n", (u64)position, (u64)(str + str_length) );
  }
  
  //-- FILL THE ACTUAL XML NOW --//
  hlXML *parent_xml  = 0;
  hlXML *current_xml = &result;
  current_tokken     = &head_tokken;

  while( current_tokken ){
    /*
    if( parent_xml ){
      write( 1, parent_xml->name, 64 );
      write( 1, " - ", 3 );
    }
    write( 1, current_xml->name, 64 );
    write( 1, " - ", 3 );
    write( 1, current_tokken->value, current_tokken->value_length );
    write( 1, "\n", 1 );
    */

    switch( current_tokken->type ){
    case hlXML_TOKKEN_ELEM_START:
    case hlXML_TOKKEN_START:{
      if( current_xml->type ){
        hlXML *new_xml = hlXMLAddChild( current_xml );
        parent_xml  = current_xml;
        current_xml = new_xml;
        //write( 1, current_tokken->value, current_tokken->value_length );
        //write( 1, "\n", 1 );
      }
      if( current_tokken->next->type == hlXML_TOKKEN_START ){
        current_xml->type = hlXML_OBJECT;
        u64 size = hlXMLCountChildren( current_tokken );
        current_xml->obj_val.size     = size;
        current_xml->obj_val.count    = 0;
        current_xml->obj_val.children = hlPUSHARRAY( pool, size, hlXML );
      }else if( current_tokken->next->type == hlXML_TOKKEN_ELEM_START ){
        current_xml->type = hlXML_ARR;
        u64 size = hlXMLCountChildren( current_tokken );
        current_xml->obj_val.size     = size;
        current_xml->obj_val.count    = 0;
        current_xml->obj_val.children = hlPUSHARRAY( pool, size, hlXML );
        //printf( "FOUND ARRAY OBJECT (%lu)\n", size );
      }else if( current_tokken->next->type == hlXML_TOKKEN_VAL ){
        current_xml->type = hlXML_STR;
      }else{
        hlASSERT( !"WRONG TOKKEN GIVEN" );
      }

      //printf( "COPYING (%lu)\n", current_tokken->value_length );
      hlMemCopy( current_tokken->value, current_xml->name
               , current_tokken->value_length < hlARRAYCOUNT( current_xml->name ) 
                 ? current_tokken->value_length : hlARRAYCOUNT( current_xml->name ) );
      //printf( "COPIED\n" );
    }break;
    case hlXML_TOKKEN_ELEM_END:
    case hlXML_TOKKEN_END:{
      current_xml = parent_xml;
      if( current_xml->parent != 0 ){
        parent_xml  = parent_xml->parent;
      }
    }break;
    case hlXML_TOKKEN_VAL:{
      u64 whitespace_count = hlCountWhitespaceC( current_tokken->value      );
      u64 dot_count        = hlCountCharacterC(  current_tokken->value, '.' );
      u64 a_count          = hlCountCharacterC(  current_tokken->value, 'a' );
      u64 e_count          = hlCountCharacterC(  current_tokken->value, 'e' );
      b32 numeric          = hlIsNumericC(       current_tokken->value      );
      if( whitespace_count == 0 ){
        if( numeric && dot_count < 2 ){
          if( dot_count == 1 || e_count == 1 ){
            //--FLT--//
            c8 buffer[256];
            hlASSERT( current_tokken->value_length + 1 < hlARRAYCOUNT( buffer ) );
            snprintf( buffer, current_tokken->value_length + 1, "%s", current_tokken->value );

            current_xml->flt_val = strtod( buffer, 0 );
            //printf( "SINGLE: %s -> %f\n", buffer, current_xml->flt_val );
            current_xml->type    = hlXML_FLT;
          }else{
            //--NUM--//
            c8 buffer[256];
            hlASSERT( current_tokken->value_length + 1 < hlARRAYCOUNT( buffer ) );
            snprintf( buffer, current_tokken->value_length + 1, "%s", current_tokken->value );

            current_xml->num_val = atol( buffer );   
            current_xml->type    = hlXML_NUM;
            //write( 1, current_tokken->value, current_tokken->value_length );
            //printf( " - %lu - %u\n", dot_count, numeric );
            //printf( "%ld\n", current_xml->num_val );
          }
        }else{
          //--STR--//
          current_xml->str  = hlPUSHSIZE( pool, current_tokken->value_length + 1 );
          current_xml->type = hlXML_STR;

          current_xml->str[ current_tokken->value_length - 1 ] = 0;
          hlMemCopy( current_tokken->value, current_xml->str, current_tokken->value_length );
        }
      }else{
        if( numeric ){
          //write( 1, current_tokken->value, current_tokken->value_length );
          //write( 1, "\n", 1 );
          c8 *pointer = current_tokken->value;
          if( dot_count == 0 && a_count == 0 ){
            //--NUM_ARR--//
            current_xml->num_arr_val.count  = whitespace_count + 1;
            current_xml->num_arr_val.values = hlPUSHARRAY( pool, current_xml->num_arr_val.count, i64 );
            current_xml->type               = hlXML_NUM_ARR;

            c8 *buffer = hlPUSHARRAY( tmp_pool, current_tokken->value_length + 1, c8 );
            hlMemCopy( current_tokken->value, buffer, current_tokken->value_length );
            buffer[ current_tokken->value_length ] = '\0';
            
            for( u64 arr_index = 0; arr_index < current_xml->num_arr_val.count; ++arr_index ){
              current_xml->num_arr_val.values[ arr_index ] = atol( buffer );
              buffer = hlStrGuardedNextSpaceC( buffer, '\0' );
            }
          }else{
            //--FLT_ARR--//
            current_xml->flt_arr_val.count  = whitespace_count + 1;
            current_xml->flt_arr_val.values = hlPUSHARRAY( pool, current_xml->flt_arr_val.count, r64 );
            current_xml->type               = hlXML_FLT_ARR;

            c8 *buffer = hlPUSHARRAY( tmp_pool, current_tokken->value_length + 1, c8 );
            hlMemCopy( current_tokken->value, buffer, current_tokken->value_length );
            buffer[ current_tokken->value_length ] = '\0';
            
            for( u64 arr_index = 0; arr_index < current_xml->num_arr_val.count; ++arr_index ){
              current_xml->flt_arr_val.values[ arr_index ] = atof( buffer );
              buffer = hlStrGuardedNextSpaceC( buffer, '\0' );
            }
          }
        }else{
          //--STR--//
          current_xml->str  = hlPUSHSIZE( pool, current_tokken->value_length + 1 );
          current_xml->type = hlXML_STR;

          current_xml->str[ current_tokken->value_length - 1 ] = 0;
          hlMemCopy( current_tokken->value, current_xml->str, current_tokken->value_length );
        }

      }
    }break;
    default:{
      hlASSERT( !"WRONG TOKKEN GIVEN" );
    }
    }
    current_tokken = current_tokken->next;
  }
  
  tmp_pool->used = tmp_pool_starting_size;
  return( result );
}

typedef struct {
  u64     count;
  hlXML **values;
} hlXML_Result;

typedef struct hlXMLTmpResult {
  struct hlXMLTmpResult *next;
  hlXML                 *val;
} hlXMLTmpResult;

static void hlXMLAddToTmpResult( hlMem_Pool *pool, hlXMLTmpResult *result, hlXML *new_entry ){
  hlXMLTmpResult *new_result = hlPUSHTYPE( pool, hlXMLTmpResult );

  hlXMLTmpResult *tmp_pointer = result;
  while( tmp_pointer->next ){
   tmp_pointer = tmp_pointer->next;
  }
  tmp_pointer->next = new_result;
  new_result->next = 0;
  new_result->val  = new_entry;
}

static hlXML_Result hlXMLQuery( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  hlXML_Result result = {};
  //printf( "%s\n", query );

  u64 query_len = hlStrLengthC( query );

  u64 slash_count = hlCountCharacterC( query, '/' );
  c8 *pointer = query;

  hlXMLTmpResult  base_query    = { 0, root };            
  hlXMLTmpResult *current_query = &base_query;
  for( u64 slash_index = 0; slash_index < slash_count + 1 && current_query; ++slash_index ){
    c8 *tokken = 0;
    u64 tokken_length;
    if( slash_index == slash_count ){
      tokken        = pointer;
      tokken_length = query_len - (u64)(pointer - query);
    }else{
      pointer = hlStrGrabTilC( tmp_pool, pointer, "/", &tokken, &tokken_length );
    }

    hlXMLTmpResult *new_query = 0;
    do{
      hlXML *current_focus = current_query->val;
      if( current_focus->type == hlXML_OBJECT ){
        for( u64 child_index = 0; child_index < current_focus->obj_val.count; ++child_index ){
          hlXML *child = current_focus->obj_val.children + child_index;
          if( hlStrCmpC( tokken, child->name ) ){
            if( new_query ){
              hlXMLAddToTmpResult( tmp_pool, new_query, child ); 
            }else{
              new_query = hlPUSHTYPE( tmp_pool, hlXMLTmpResult );
              new_query->next = 0;
              new_query->val  = child;
            }
          }
        }
      }else if( current_focus->type == hlXML_ARR ){
        for( u64 child_index = 0; child_index < current_focus->obj_val.count; ++child_index ){
          hlXML *child = current_focus->obj_val.children + child_index;
          for( u64 child_of_child_index = 0; child_of_child_index < child->obj_val.count; ++child_of_child_index ){
            hlXML *child_of_child = child->obj_val.children + child_of_child_index;
            if( hlStrCmpC( tokken, child_of_child->name ) ){
              if( new_query ){
                hlXMLAddToTmpResult( tmp_pool, new_query, child_of_child ); 
              }else{
                new_query = hlPUSHTYPE( tmp_pool, hlXMLTmpResult );
                new_query->next = 0;
                new_query->val  = child_of_child;
              }
            }
          }
        }

      }else{
        hlASSERT( !"FOCUSED ON WRONG TYPE OF OBJECT CHECK QUERY" );
      }
      current_query = current_query->next;
    }while( current_query );
    current_query = new_query;
  }

  u64 result_count = 0;
  if( current_query ){
    ++result_count;
    hlXMLTmpResult *query_pointer = current_query;
    while( query_pointer->next ){
      query_pointer = query_pointer->next;
      ++result_count;
    }
  }

  if( result_count ){
    result.count  = result_count;
    result.values = hlPUSHARRAY( pool, result.count, hlXML *);
    for( u64 result_index = 0; result_index < result.count; ++result_index ){
      result.values[ result_index ] = current_query->val;
      current_query = current_query->next;
    }
  }

  return( result );
}

static void hlXMLPrintQueryResult( hlXML_Result *result ){
  printf( "XML RESULT (%llu)\n", result->count );
  for( u64 i = 0; i < result->count; ++i ){
    write( 1, result->values[ i ]->name, sizeof( result->values[ i ]->name ) );
    write( 1, "\n", 1 );
  }
}

static hlXML *hlXMLQuerySingle( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query, hl_xml_type type ){
  hlXML *result = 0;
  hlXML_Result query_result = hlXMLQuery( tmp_pool, pool, root, query );
  //hlXMLPrintQueryResult( &query_result );
  hlASSERT( query_result.count == 1 );
  hlASSERT( query_result.values[ 0 ]->type == type );
  result = query_result.values[ 0 ];
  return( result );
}

static hlXML *hlXMLQueryObject( hlMem_Pool *tmp, hlMem_Pool *mem, hlXML *root, c8 *query ){
  hlXML *result = hlXMLQuerySingle( tmp, mem, root, query, hlXML_OBJECT );
  return( result );
}

static c8 *hlXMLQueryStr( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  c8 *result = 0;
  hlXML *query_result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_STR );
  result = query_result->str;
  return( result );
}

static i64 hlXMLQueryNum( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  i64 result = 0;
  hlXML *query_result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_NUM );
  result = query_result->num_val;
  return( result );
}

static r64 hlXMLQueryFlt( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  r64 result = 0;
  hlXML *query_result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_FLT );
  result = query_result->flt_val;
  return( result );
}

static hlXML *hlXMLQueryArr( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  hlXML *result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_ARR );
  hlASSERT( result->obj_val.count > 0 );
  return( result );
}

static i64 *hlXMLQueryNumVec( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query, u64 count ){
  i64 *result = 0;
  hlXML *query_result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_NUM_ARR );
  hlASSERT( query_result->num_arr_val.count == count );
  result = query_result->num_arr_val.values;
  return( result );
}

static r64 *hlXMLQueryFltVec( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query, u64 count ){
  r64 *result = 0;
  hlXML *query_result = hlXMLQuerySingle( tmp_pool, pool, root, query, hlXML_FLT_ARR );
  hlASSERT( query_result->flt_arr_val.count == count );
  result = query_result->flt_arr_val.values;
  return( result );
}

static void hlXMLPrintChildren( hlXML *parent ){
  hlASSERT( parent->type == hlXML_OBJECT );
  write( 1, parent->name, sizeof( parent->name ) );
  c8 buffer[256];
  u64 size;
  size = snprintf( buffer, hlARRAYCOUNT( buffer ), "(%s)", hl_xml_type_to_str[ parent->type ] );
  write( 1, buffer, size );

  write( 1, "\n", 1 );
  for( u64 child_count = 0; child_count < parent->obj_val.count; ++child_count ){
    hlXML *child = parent->obj_val.children + child_count;
    write( 1, "  -  ", 5 );
    write( 1, child->name, sizeof( child->name ) );
    size = snprintf( buffer, hlARRAYCOUNT( buffer ), "(%s)", hl_xml_type_to_str[ child->type ] );
    write( 1, buffer, size );
    write( 1, "\n", 1 );
  }
}

#if defined( hlDEBUG )
#include <stdio.h>
#endif

typedef enum hl_json_type {
  hlJSON_DOUBLE
, hlJSON_NUMBER
, hlJSON_STRING
, hlJSON_ARRAY
, hlJSON_OBJECT
, hlJSON_BOOL
, hlJSON_NULL
} hl_json_type;

typedef struct hlJSON_Object {
  u64 count;
#if defined( hlDEBUG )
  u64 size;
#endif
         c8   **children_names;
  struct hlJSON *children;
} hlJSON_Object;

typedef struct hlJSON_Array {
  u64 count;
#if defined( hlDEBUG )
  u64 size;
#endif
  struct hlJSON *children;
} hlJSON_Array;


typedef struct hlJSON {
  hl_json_type   type;
  union {
    struct hlJSON_Object obj_val;
    struct hlJSON_Array  arr_val;
           r64           dbl_val;
           i64           num_val;
           c8           *str_val;
           b32           bool_val;
  };
} hlJSON;

static void hlJSONGenerateObject( hlMem_Pool *pool, hlJSON *result, u64 number_children ){
  result->type                   = hlJSON_OBJECT;
  result->obj_val.count          = 0;
#if defined( hlDEBUG )
  result->obj_val.size           = number_children;
#endif
  result->obj_val.children_names = hlPUSHARRAY( pool, number_children, c8 *   );
  result->obj_val.children       = hlPUSHARRAY( pool, number_children, hlJSON );
}

static void hlJSONGenerateArray( hlMem_Pool *pool, hlJSON *result, u64 number_children ){
  result->type                   = hlJSON_ARRAY;
  result->arr_val.count          = 0;
#if defined( hlDEBUG )
  result->arr_val.size           = number_children;
#endif
  result->arr_val.children       = hlPUSHARRAY( pool, number_children, hlJSON );
}

static void hlJSONPushStringField( hlJSON *parent, c8 *name, c8 *value ){
  hlASSERT( parent->type == hlJSON_OBJECT );
  u64 count = parent->obj_val.count++;
  parent->obj_val.children_names[ count ] = name;

  parent->obj_val.children[ count ].type    = hlJSON_STRING;
  parent->obj_val.children[ count ].str_val = value;
  hlASSERT( parent->obj_val.count <= parent->obj_val.size );
}

static void hlJSONPushDoubleField( hlJSON *parent, c8 *name, r64 value ){
  hlASSERT( parent->type == hlJSON_OBJECT );
  u64 count = parent->obj_val.count++;
  parent->obj_val.children_names[ count ] = name;

  parent->obj_val.children[ count ].type    = hlJSON_DOUBLE;
  parent->obj_val.children[ count ].dbl_val = value;
  hlASSERT( parent->obj_val.count <= parent->obj_val.size );
}

static void hlJSONPushNumberField( hlJSON *parent, c8 *name, i64 value ){
  hlASSERT( parent->type == hlJSON_OBJECT );
  u64 count = parent->obj_val.count++;
  parent->obj_val.children_names[ count ] = name;

  parent->obj_val.children[ count ].type    = hlJSON_NUMBER;
  parent->obj_val.children[ count ].num_val = value;
  hlASSERT( parent->obj_val.count <= parent->obj_val.size );
}

static hlJSON *hlJSONPushObjectField( hlMem_Pool *pool, hlJSON *parent, c8 *name, u64 number_children ){
  hlASSERT( parent->type == hlJSON_OBJECT );
  u64 count = parent->obj_val.count++;
  parent->obj_val.children_names[ count ] = name;

  hlJSONGenerateObject( pool, parent->obj_val.children + count, number_children );
  hlASSERT( parent->obj_val.count <= parent->obj_val.size );

  return( parent->obj_val.children + count );
}

static hlJSON *hlJSONPushArrayField( hlMem_Pool *pool, hlJSON *parent, c8 *name, u64 number_children ){
  hlASSERT( parent->type == hlJSON_OBJECT );
  u64 count = parent->obj_val.count++;
  parent->obj_val.children_names[ count ] = name;

  hlJSONGenerateArray( pool, parent->obj_val.children + count, number_children );

  hlASSERT( parent->obj_val.count <= parent->obj_val.size );
  return( parent->obj_val.children + count );
}

static void hlJSONPushStringValue( hlJSON *parent, c8 *value ){
  hlASSERT( parent->type == hlJSON_ARRAY );
  u64 count = parent->arr_val.count++;

  parent->arr_val.children[ count ].type    = hlJSON_STRING;
  parent->arr_val.children[ count ].str_val = value;
  hlASSERT( parent->arr_val.count <= parent->arr_val.size );
}

static void hlJSONPushDoubleValue( hlJSON *parent, r64 value ){
  hlASSERT( parent->type == hlJSON_ARRAY );
  u64 count = parent->arr_val.count++;

  parent->arr_val.children[ count ].type    = hlJSON_DOUBLE;
  parent->arr_val.children[ count ].dbl_val = value;
  hlASSERT( parent->arr_val.count <= parent->arr_val.size );
}

static void hlJSONPushNumberValue( hlJSON *parent, i64 value ){
  hlASSERT( parent->type == hlJSON_ARRAY );
  u64 count = parent->arr_val.count++;

  parent->arr_val.children[ count ].type    = hlJSON_NUMBER;
  parent->arr_val.children[ count ].num_val = value;
  hlASSERT( parent->arr_val.count <= parent->arr_val.size );
}

static hlJSON *hlJSONPushObjectValue( hlMem_Pool *pool, hlJSON *parent, u64 number_children ){
  hlASSERT( parent->type == hlJSON_ARRAY );
  u64 count = parent->arr_val.count++;

  hlJSONGenerateObject( pool, parent->arr_val.children + count, number_children );
  hlASSERT( parent->arr_val.count <= parent->arr_val.size );

  return( parent->arr_val.children + count );
}

static hlJSON *hlJSONPushArrayValue( hlMem_Pool *pool, hlJSON *parent, u64 number_children ){
  hlASSERT( parent->type == hlJSON_ARRAY );
  u64 count = parent->arr_val.count++;

  hlJSONGenerateArray( pool, parent->arr_val.children + count, number_children );
  hlASSERT( parent->arr_val.count <= parent->arr_val.size );

  return( parent->arr_val.children + count );
}

static void _hlJSONPrintIdent( FILE *file,  u64 ident_level, b32 comma, b32 ident ){
  if( comma ){
    if( ident_level ){
      --ident_level;
    }
    fprintf( file, ", " );
  }
  if( ident ){
    for( u64 ident_index = 0; ident_index < ident_level; ++ident_index ){
      fprintf( file,  "  " );
    }
  }
}

static void _hlJSONPrint( FILE *file, hlJSON *object, hl_json_type parent_type, u64 ident_level, b32 first, b32 beginning_ident ){
  switch( object->type ){
    case hlJSON_OBJECT:{
      if( parent_type == hlJSON_ARRAY ){
        _hlJSONPrintIdent( file, ident_level - 1, !first, beginning_ident );
      }else{
        _hlJSONPrintIdent( file, ident_level    , !first, beginning_ident );
      }
      fprintf( file, "{\n" );

      b32 first_child = 1;
      for( u64 object_index = 0; object_index < object->obj_val.count; ++object_index ){
        c8     *child_name = object->obj_val.children_names[ object_index ];
        hlJSON *child      = object->obj_val.children + object_index;
        hlASSERT( child );

        _hlJSONPrintIdent( file, ident_level + 1, !first_child, 1 );
        fprintf( file, "\"%s\" : ", child_name );
        _hlJSONPrint( file, child, object->type, ident_level + 1, 1, !first_child );
        fprintf( file, "\n" );
        first_child = 0;
      }

      _hlJSONPrintIdent( file, ident_level, 0, 1 );
      fprintf( file, "}\n");
    }break;
    case hlJSON_ARRAY:{
      _hlJSONPrintIdent( file, ident_level, !first, !first );
      fprintf( file, "[ " );

      b32 first_child = 1;
      hlASSERT( object->arr_val.children );
      for( u64 array_index = 0; array_index < object->arr_val.count; ++array_index ){
        hlJSON *child = object->arr_val.children + array_index;
        hlASSERT( child );

        _hlJSONPrintIdent( file, 0, !first_child, !first_child );
        _hlJSONPrint( file, child, object->type, ident_level + 1, 1, !first_child );
        first_child = 0;
      }

      fprintf( file, " ]");
      if( object->type == hlJSON_OBJECT ){
        fprintf( file, "\n" );
      }
    }break;
    case hlJSON_DOUBLE:{
      fprintf( file, "%.13e", object->dbl_val );
      //hlDbl_Bits b;
      //b.dbl = object->dbl_val;
      //printf( "0x%016llx\n", b.bits );
    }break;
    case hlJSON_NUMBER:{
      fprintf( file, "%lld", object->num_val );
    }break;
    case hlJSON_STRING:{
      fprintf( file, "\"%s\"", object->str_val );
    }break;
    default:{
    }break;
  }
}

static void hlJSONPrint( FILE *file, hlJSON *object ){
  _hlJSONPrint( file, object, 0, 0, 1, 0 );
}
//-------------------------------------------------------------------------------------------//
//--                                      HASH MAP                                         --//
//-------------------------------------------------------------------------------------------//

typedef struct hlHash_Map_Item {
  u64  hash;
  void   *item;
} hlHash_Map_Item;

typedef struct hlHash_Map {
  u64 slotcount;
  u64 popcount;
  hlHash_Map_Item *slots;
} hlHash_Map;

static hlHash_Map hlHashMapGenerate( hlMem_Pool *pool, u64 slotcount ){
  hlHash_Map result = {};
             result.slotcount = slotcount;
             result.slots = hlPUSHSIZE( pool, slotcount*sizeof(hlHash_Map_Item) );
  return( result );
}

static void hlHashMapInsert( hlHash_Map *map, u64 hash, void *item ){
  hlASSERT( hash );
  hlASSERT( 3*map->popcount < map->slotcount );
  u64 preferable_index = hash % map->slotcount;
  for( u64 i = 0; i < map->slotcount; ++i ){
    u64 test_index = (preferable_index + i) % map->slotcount;
    hlASSERT( map->slots[test_index].hash != hash );
    if( map->slots[test_index].hash == 0 ){
      map->slots[test_index].hash = hash;
      map->slots[test_index].item = item;
      ++map->popcount;
      break;
    }
  }
}

static void *_hlHashMapGet( hlHash_Map *map, u64 hash ){
  hlASSERT( hash );
  void *result = 0;
  u64 preferable_index = hash % map->slotcount;
  for( u64 i = 0; i < map->slotcount; ++i ){
    u64 test_index = (preferable_index + i) % map->slotcount;
    if( map->slots[test_index].hash == 0 ){
      break;
    }else if( map->slots[test_index].hash == hash ){
      result = map->slots[test_index].item;
      break;
    }
  }
  return( result );
}
#define hlHashMapGet( map, hash, type ) (type *)_hlHashMapGet( map, hash )

#define _HARMLESS_H
#endif
